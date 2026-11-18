# Build the whole-cohort statistics used by the public Table 2 analysis.
# Run from the repository root. Supply the session, steady-segment and
# runner-day Parquet files in that order. The fourth argument is the destination.
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript scripts/prepare/build_model_statistics.R SESSIONS SEGMENTS RUNNER_DAYS DESTINATION")
}
source("scripts/model_tools.R", encoding = "UTF-8")
setDTthreads(2L)
db <- DBI::dbConnect(duckdb::duckdb(home = tempdir()))
read_columns <- function(path, columns) {
  names_sql <- paste(DBI::dbQuoteIdentifier(db, columns), collapse = ", ")
  query <- paste("SELECT", names_sql, "FROM read_parquet(", DBI::dbQuoteString(db,path), ")")
  as.data.table(DBI::dbGetQuery(db, query))
}

variables <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean",
               "start_minutes_local","start_before_09","start_after_19","session_day_percentile",
               "log_distance","steady_speed","run")
sessions <- read_columns(args[1],c("session_hash","user_hash","day_max_c","canonical_start_utc",
  "geo_mean_y","local_weekday",variables[1:9]))[is.finite(day_max_c)]
segments <- read_columns(args[2],c("session_hash","speed_lag60"))
days <- read_columns(args[3],c("user_hash","year_month_local","day_max_c","local_weekday","run"))[is.finite(day_max_c)]
DBI::dbDisconnect(db, shutdown=TRUE)
stopifnot(!anyDuplicated(sessions$session_hash),!anyDuplicated(segments$session_hash))
sessions[segments, on="session_hash", steady_speed:=i.speed_lag60]

# Check the overlap of all actual session samples. The differences between
# samples must also pool at least ten distinct runners. Do not export counts.
sample_masks <- as.matrix(sessions[,lapply(.SD, is.finite),.SDcols=variables[1:10]])
strata <- apply(sample_masks, 1L, paste, collapse="")
stopifnot(all(sessions[,uniqueN(user_hash),by=.(stratum=strata)]$V1 >= 10L))
percentiles <- c(1,10,50,90,95,99)
anchors <- data.table(percentile=percentiles,
  temperature=as.numeric(quantile(sessions$day_max_c,percentiles/100)))
knots <- anchors[percentile %in% c(10,50,90),temperature]
session_terms <- cbind(rcs_basis(sessions$day_max_c,knots),
  season_terms(sessions$canonical_start_utc,sessions$geo_mean_y<0),
  model.matrix(~factor(local_weekday,levels=0:6),sessions)[,-1,drop=FALSE])
day_terms <- cbind(rcs_basis(days$day_max_c,knots),
  model.matrix(~factor(local_weekday,levels=0:6),days)[,-1,drop=FALSE])

moments <- list()
for (variable in variables) {
  if (variable=="run") {
    design <- day_terms
    y <- as.numeric(days$run)
    runner <- days$user_hash
    effect <- paste(runner,days$year_month_local)
  } else {
    design <- session_terms
    y <- as.numeric(sessions[[variable]])
    runner <- sessions$user_hash
    effect <- runner
  }
  keep <- is.finite(y) & complete.cases(design) & !is.na(effect) & !is.na(runner)
  all_terms <- cbind(design[keep,,drop=FALSE],outcome=y[keep])
  id <- match(effect[keep],unique(effect[keep]))
  runner <- runner[keep]
  n <- nrow(all_terms)
  g <- uniqueN(runner)
  group_means <- rowsum(all_terms,id,reorder=FALSE) / tabulate(id)
  within <- all_terms - group_means[id,,drop=FALSE]

  # Residualise exposure and outcome against the same nuisance predictors.
  temperature <- c("linear","cubic1")
  nuisance <- setdiff(colnames(design),temperature)
  adjusted <- qr.resid(qr(within[,nuisance,drop=FALSE]),within[,c(temperature,"outcome"),drop=FALSE])
  x <- adjusted[,temperature,drop=FALSE]
  y <- adjusted[,"outcome"]
  gram <- crossprod(x)/n
  rhs <- crossprod(x,y)/n
  beta <- solve(gram,rhs)
  residual <- as.numeric(y - x %*% beta)
  scores <- rowsum(x*residual,runner,reorder=FALSE)
  rank <- qr(within[,colnames(design),drop=FALSE])$rank
  meat <- (g/(g-1)) * ((n-1)/(n-rank)) * crossprod(scores)/n^2

  # Each published moment must have nonzero contributions from >=10 runners.
  design_moments <- cbind(x[,1]^2,x[,1]*x[,2],x[,2]^2,x[,1]*y,x[,2]*y)
  by_runner <- rowsum(design_moments,runner,reorder=FALSE)
  by_runner <- cbind(by_runner,scores[,1]^2,scores[,1]*scores[,2],scores[,2]^2)
  stopifnot(all(colSums(abs(by_runner)>1e-12)>=10L))
  moments[[variable]] <- data.table(outcome=variable,xx11=gram[1,1],xx12=gram[1,2],xx22=gram[2,2],
    xy1=rhs[1],xy2=rhs[2],meat11=meat[1,1],meat12=meat[1,2],meat22=meat[2,2])
}
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
fwrite(rbindlist(moments),file.path(args[4],"moments.csv"))
fwrite(anchors,file.path(args[4],"temperature_percentiles.csv"))
