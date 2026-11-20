# Build pooled statistics for daily activities, the distributed-lag running
# model and weekly activity volumes. Run from the repository root.
library(data.table)
setDTthreads(2L)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5L) stop("Supply SESSIONS, RUNNER_DAYS, WORKOUT_METADATA, DAILY_WEATHER and DESTINATION.")
source("scripts/model_tools.R",encoding="UTF-8")
db <- DBI::dbConnect(duckdb::duckdb(home=tempdir()))
read_columns <- function(path,columns) {
  cols <- paste(DBI::dbQuoteIdentifier(db,columns),collapse=", ")
  as.data.table(DBI::dbGetQuery(db,paste("SELECT",cols,"FROM read_parquet(",DBI::dbQuoteString(db,path),")")))
}
s <- read_columns(args[1],c("user_hash","local_date","cell_10km","gps_iana_zone","day_max_c",
                            "duration_minutes","gps_distance_km"))
d <- read_columns(args[2],c("user_hash","local_date","noon_utc","south","local_weekday","day_max_c","run"))[is.finite(day_max_c)]
m <- read_columns(args[3],c("user_hash","sport","start_utc","metadata_signature_redundant"))
w <- read_columns(args[4],c("user_hash","history_anchor_date","daily_max_c"))
DBI::dbDisconnect(db,shutdown=TRUE)
d[,date:=as.IDate(local_date)]
d[,time:=as.numeric(noon_utc)]
stopifnot(!anyDuplicated(d,by=c("user_hash","date")))
knots <- as.numeric(quantile(s$day_max_c,c(.1,.5,.9),na.rm=TRUE))

# Assign each workout to a civil date in the runner's usual time zone.
locations <- s[!is.na(cell_10km),.N,by=.(user_hash,cell_10km)]
setorder(locations,user_hash,-N,cell_10km)
usual <- locations[!duplicated(user_hash)]
zones <- s[usual,on=.(user_hash,cell_10km),nomatch=0L][!is.na(gps_iana_zone),.N,by=.(user_hash,gps_iana_zone)]
setorder(zones,user_hash,-N,gps_iana_zone)
zones <- zones[!duplicated(user_hash)]
m <- m[user_hash %in% d$user_hash & metadata_signature_redundant==FALSE &
         between(start_utc,946684800,1514764800)]
m[zones,on="user_hash",zone:=i.gps_iana_zone]
stopifnot(!anyNA(m$zone))
m[,date:=as.IDate(format(as.POSIXct(as.numeric(start_utc),origin="1970-01-01",tz="UTC"),
                          "%Y-%m-%d",tz=zone[1L])),by=zone]
workouts <- m[,.(any_workouts=.N,outdoor_cycling=as.integer(any(sport %in% c("bike","bike (transport)","mountain bike")))),
              by=.(user_hash,date)]
d[, `:=`(any_workouts=0L,outdoor_cycling=0L)]
d[workouts,on=.(user_hash,date),`:=`(any_workouts=i.any_workouts,outdoor_cycling=i.outdoor_cycling)]
d[,any_sport:=as.integer(any_workouts>0)]
running <- s[,.(minutes=sum(duration_minutes,na.rm=TRUE),km=sum(gps_distance_km,na.rm=TRUE)),
             by=.(user_hash,date=as.IDate(local_date))]
d[, `:=`(minutes=0,km=0)]
d[running,on=.(user_hash,date),`:=`(minutes=i.minutes,km=i.km)]
d[,iso_week:=format(date,"%G-%V")]
d[,complete_week:=.N==7L,by=.(user_hash,iso_week)]
bins <- c("below_0","0_10","20_25","25_30","30_plus")
d[, (bins):=list(day_max_c<0,day_max_c>=0 & day_max_c<10,day_max_c>=20 & day_max_c<25,
                 day_max_c>=25 & day_max_c<30,day_max_c>=30)]
# Each runner-day contributes at most once to "runs", even if it includes several runs.
weekly <- d[complete_week==TRUE,c(list(runs=sum(run),minutes=sum(minutes),km=sum(km),
  any_workouts=sum(any_workouts),south=first(south),time=mean(time)),lapply(.SD,sum)),
  by=.(user_hash,iso_week),.SDcols=bins]

# Look up each temperature in its own runner's weather history.
w[,date:=as.IDate(history_anchor_date)]
setkey(w,user_hash,date)
relative_day <- -3:7
lag_names <- ifelse(relative_day<0,paste0("lead",abs(relative_day)),paste0("lag",relative_day))
for (i in seq_along(relative_day)) {
  date_lookup <- d$date-relative_day[i]
  set(d,j=lag_names[i],value=w[.(d$user_hash,date_lookup),daily_max_c])
}
stopifnot(max(abs(d$lag0-d$day_max_c),na.rm=TRUE)<1e-6)
d[,complete_lags:=complete.cases(.SD),.SDcols=lag_names]
lag_days <- d[complete_lags==TRUE]
# Every nonempty activity and sample-overlap cell must contain >=10 runners.
stopifnot(all(d[,uniqueN(user_hash),by=.(complete_week,complete_lags,any_sport,outdoor_cycling)]$V1>=10L))

# Remove runner intercepts and nuisance covariates. Retain only the pooled
# normal equations and CR1 covariance of the fitted model's cluster scores.
pool <- function(data, y_name, x, daily) {
  nuisance <- season_terms(data$time,data$south)
  if (daily) nuisance <- cbind(nuisance,model.matrix(~factor(local_weekday,levels=0:6),data)[,-1,drop=FALSE])
  q <- cbind(x,outcome=as.numeric(data[[y_name]]),nuisance)
  id <- match(data$user_hash,unique(data$user_hash))
  within <- q - (rowsum(q,id,reorder=FALSE)/tabulate(id))[id,,drop=FALSE]
  desired <- c(colnames(x),"outcome")
  adjusted <- qr.resid(qr(within[,colnames(nuisance),drop=FALSE]),within[,desired,drop=FALSE])
  u <- adjusted[,colnames(x),drop=FALSE]
  y <- adjusted[,"outcome"]
  n <- length(y)
  g <- max(id)
  gram <- crossprod(u)/n
  rhs <- crossprod(u,y)/n
  beta <- solve(gram,rhs)
  scores <- rowsum(u*as.numeric(y-u%*%beta),id,reorder=FALSE)
  rank <- qr(within[,c(colnames(x),colnames(nuisance)),drop=FALSE])$rank
  meat <- (g/(g-1))*((n-1)/(n-rank))*crossprod(scores)/n^2
  # Check the contributors to each published matrix/vector element, including
  # the least common weekly weather bands. Runner contributions stay in memory.
  for (j in seq_len(ncol(u))) {
    support <- rowsum(u[,j:ncol(u),drop=FALSE]*u[,j],id,reorder=FALSE)
    stopifnot(all(colSums(abs(support)>1e-12)>=10L))
    stopifnot(all(colSums(abs(scores[,j:ncol(u),drop=FALSE]*scores[,j])>1e-12)>=10L))
  }
  stopifnot(all(colSums(abs(rowsum(u*y,id,reorder=FALSE))>1e-12)>=10L))
  ij <- which(upper.tri(gram,diag=TRUE),arr.ind=TRUE)
  list(matrix=data.table(term1=colnames(x)[ij[,1]],term2=colnames(x)[ij[,2]],xx=gram[ij],meat=meat[ij]),
       vector=data.table(term=colnames(x),xy=as.numeric(rhs)))
}
fits <- list()
x <- rcs_basis(d$day_max_c,knots)
for (sport in c("any_sport","outdoor_cycling")) fits[[sport]] <- pool(d,sport,x,TRUE)
x <- do.call(cbind,lapply(lag_names,function(nm) {
  v <- rcs_basis(lag_days[[nm]],knots); colnames(v) <- paste(nm,colnames(v),sep="_"); v
}))
fits$running_lags <- pool(lag_days,"run",x,TRUE)
x <- as.matrix(weekly[,..bins])
for (v in c("runs","minutes","km","any_workouts")) fits[[paste0("weekly_",v)]] <- pool(weekly,v,x,FALSE)
dir.create(args[5],recursive=TRUE,showWarnings=FALSE)
fwrite(rbindlist(lapply(fits,"[[","matrix"),idcol="model"),file.path(args[5],"pooled_matrices.csv"))
fwrite(rbindlist(lapply(fits,"[[","vector"),idcol="model"),file.path(args[5],"pooled_vectors.csv"))
