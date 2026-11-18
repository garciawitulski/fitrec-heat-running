# Prepare whole-cohort statistics for the fixed-hour and physiological models.
# Inputs are the licensed sessions and steady-speed panels and the study's
# model configuration. The released files contain no individual scores.
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4L) stop("Usage: Rscript scripts/prepare/build_decomposition_statistics.R SESSIONS STEADY_SPEED CONFIG DESTINATION")
source("scripts/model_tools.R",encoding="UTF-8")
db <- DBI::dbConnect(duckdb::duckdb(home=tempdir()),dbdir=":memory:")
read_fields <- function(path,cols) {
  sql <- paste("SELECT",paste(DBI::dbQuoteIdentifier(db,cols),collapse=","),
               "FROM read_parquet(",DBI::dbQuoteString(db,path),")")
  as.data.table(DBI::dbGetQuery(db,sql))
}
common <- c("user_hash","canonical_start_utc","geo_mean_y")
main <- read_fields(args[1],c(common,"day_max_c","local_weekday","start_minutes_local","heat_avoided_c"))[is.finite(day_max_c)]
phys <- read_fields(args[2],c(common,"temp_c_10km","speed_lag60","hr_lag60",
                              "prior_recorded_minutes","fitness_speed_kmh_28d"))
DBI::dbDisconnect(db,shutdown=TRUE)
cfg <- jsonlite::fromJSON(args[3])
stopifnot(nrow(main)==106600L,nrow(phys)==19350L,all(complete.cases(main)),all(complete.cases(phys)))
daily_knots <- cfg$knots$rcs3
session_knots <- cfg$knots$session_temp_rcs3
speed_knots <- cfg$speed_knots
daily_basis <- rcs_basis(main$day_max_c,daily_knots)
time_controls <- cbind(season_terms(main$canonical_start_utc,main$geo_mean_y<0),
                       model.matrix(~factor(local_weekday,levels=0:6),main)[,-1,drop=FALSE])
hour <- model.matrix(~factor(floor(start_minutes_local/60)),main)[,-1,drop=FALSE]
design <- rhs <- adjusted_scores <- list()
for (label in c("gap","hour","phys")) {
  if (label=="phys") {
    X <- cbind(rcs_basis(phys$temp_c_10km,session_knots),
               rcs_basis(phys$speed_lag60,speed_knots,center=0,scale=1))
    Z <- cbind(season_terms(phys$canonical_start_utc,phys$geo_mean_y<0),
               phys$prior_recorded_minutes,phys$fitness_speed_kmh_28d)
    y <- phys$hr_lag60
    runner <- phys$user_hash
  } else {
    X <- daily_basis
    Z <- if (label=="hour") cbind(time_controls,hour) else time_controls
    y <- main$heat_avoided_c
    runner <- main$user_hash
  }
  id <- match(runner,unique(runner))
  w <- tabulate(id)
  centre <- function(a) a - (rowsum(a,id,reorder=FALSE)/w)[id,,drop=FALSE]
  nuisance <- qr(centre(Z))
  x <- qr.resid(nuisance,centre(X))
  y <- qr.resid(nuisance,centre(matrix(y,ncol=1)))
  A <- crossprod(x)
  b <- crossprod(x,y)
  e <- as.numeric(y-x%*%solve(A,b))
  s <- rowsum(x*e,id,reorder=FALSE)
  N <- nrow(x)
  G <- length(w)
  p <- ncol(x)
  factor <- G/(G-1)*(N-1)/(N-p-nuisance$rank)
  adjusted_scores[[label]] <- sqrt(factor)*(s%*%solve(A))
  # Check contributions before reducing to cohort-wide products.
  stopifnot(G>=10L,all(colSums(abs(s)>1e-12)>=10L))
  for (j in seq_len(p)) for (k in seq_len(p)) {
    stopifnot(sum(abs(rowsum(x[,j]*x[,k],id,reorder=FALSE))>1e-12)>=10L)
    stopifnot(sum(abs(rowsum(x[,j]*c(y),id,reorder=FALSE))>1e-12)>=10L)
  }
  if (label!="gap") {
    design[[label]] <- data.table(model=label,row=rep(seq_len(p),p),column=rep(seq_len(p),each=p),
                                  value=as.numeric(A/N))
    rhs[[label]] <- data.table(model=label,row=seq_len(p),value=as.numeric(b/N))
  }
}
pair_scores <- cbind(adjusted_scores$gap,adjusted_scores$hour)
stopifnot(all(outer(seq_len(ncol(pair_scores)),seq_len(ncol(pair_scores)),
  Vectorize(function(j,k) sum(abs(pair_scores[,j]*pair_scores[,k])>1e-15)>=10L))))
V1 <- crossprod(pair_scores)
V2 <- crossprod(adjusted_scores$phys)
covariances <- rbind(data.table(block="gap_and_hour",row=rep(1:4,4),column=rep(1:4,each=4),value=as.numeric(V1)),
                     data.table(block="phys",row=rep(1:4,4),column=rep(1:4,each=4),value=as.numeric(V2)))
reference <- as.numeric(rcs_basis(median(main$day_max_c),daily_knots))-colMeans(daily_basis)
constants <- c(mean_gap=mean(main$heat_avoided_c),reference_linear=unname(reference[1]),
               reference_cubic=unname(reference[2]),speed_reference=cfg$typical_segment_speed_m_s,
               setNames(session_knots,paste0("session_knot",1:3)),setNames(speed_knots,paste0("speed_knot",1:3)))
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
fwrite(rbindlist(design),file.path(args[4],"design.csv"))
fwrite(rbindlist(rhs),file.path(args[4],"outcomes.csv"))
fwrite(covariances,file.path(args[4],"coefficient_covariance.csv"))
fwrite(data.table(term=names(constants),value=as.numeric(constants)),file.path(args[4],"constants.csv"))
