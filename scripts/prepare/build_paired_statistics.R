# Run from the repository root with the prepared steady-speed panel, the main
# sessions file, the model configuration JSON and a destination directory.
# Only whole-cohort cross-products leave this script.
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript scripts/prepare/build_paired_statistics.R STEADY_SPEED SESSIONS CONFIG DESTINATION")
}
source("scripts/model_tools.R", encoding = "UTF-8")
db <- DBI::dbConnect(duckdb::duckdb(home = tempdir()), dbdir = ":memory:")
read_fields <- function(path, fields) {
  select <- paste(DBI::dbQuoteIdentifier(db, fields), collapse = ", ")
  sql <- paste("SELECT", select, "FROM read_parquet(", DBI::dbQuoteString(db, path), ")")
  as.data.table(DBI::dbGetQuery(db,sql))
}
d <- read_fields(args[1],c("user_hash","day_max_c","hr_lag60","session_hr_mean","speed_lag60",
  "canonical_start_utc","geo_mean_y","prior_recorded_minutes","fitness_speed_kmh_28d"))
d <- d[is.finite(day_max_c)]
temps <- read_fields(args[2],"day_max_c")[is.finite(day_max_c),day_max_c]
DBI::dbDisconnect(db, shutdown = TRUE)
cfg <- jsonlite::fromJSON(args[3])
stopifnot(nrow(d)==18626L,length(temps)==106600L,all(complete.cases(d)))
targets <- quantile(temps,c(.50,.90,.95),names=FALSE)
T <- rcs_basis(d$day_max_c,cfg$knots$rcs3)
season <- season_terms(d$canonical_start_utc,d$geo_mean_y<0)
controls <- cbind(season,training=d$prior_recorded_minutes,fitness=d$fitness_speed_kmh_28d)
speed <- rcs_basis(d$speed_lag60,cfg$speed_knots,center=0,scale=1)
runner <- match(d$user_hash,unique(d$user_hash))
runner_n <- tabulate(runner)
G <- length(runner_n)
# All three fits use exactly the same sessions and runners. The two outcomes
# must also differ for at least ten runners so their contrast cannot isolate
# a small group. No runner-specific quantity is written to disk.
stopifnot(G>=660L,G<670L,uniqueN(runner[d$hr_lag60!=d$session_hr_mean])>=10L)
within_runner <- function(x) {
  x - (rowsum(x,runner,reorder=FALSE)/runner_n)[runner,,drop=FALSE]
}
design_stats <- rhs_stats <- scores <- list()
for (model in c("fixed","free","whole")) {
  nuisance <- if (model=="fixed") cbind(controls,speed) else controls
  z <- within_runner(nuisance)
  decomp <- qr(z)
  x <- qr.resid(decomp,within_runner(T))
  y <- if (model=="whole") d$session_hr_mean else d$hr_lag60
  y <- as.numeric(qr.resid(decomp,within_runner(matrix(y,ncol=1))))
  A <- crossprod(x)
  b <- as.numeric(crossprod(x,y))
  beta <- solve(A,b)
  residual <- as.numeric(y-x %*% beta)
  score <- rowsum(x*residual,runner,reorder=FALSE)
  rank <- decomp$rank+2L
  correction <- (G/(G-1))*((nrow(d)-1)/(nrow(d)-rank))
  stopifnot(all(colSums(abs(score)>1e-10)>=10L))
  scores[[model]] <- sqrt(correction)*score
  design_stats[[model]] <- data.table(model,row=rep(1:2,2),column=rep(1:2,each=2),value=as.numeric(A))
  rhs_stats[[model]] <- data.table(model,row=1:2,value=b)
}
# The common scale factor comes after clustering. Taking cross-products of
# the six score columns retains covariance between the three fitted models.
scores <- do.call(cbind,scores)
meat <- crossprod(scores)
stopifnot(all(is.finite(meat)),nrow(meat)==6L)
contrast_basis <- rcs_basis(targets,cfg$knots$rcs3)
contrast_basis <- sweep(contrast_basis[2:3,,drop=FALSE],2,contrast_basis[1,],"-")
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
fwrite(rbindlist(design_stats),file.path(args[4],"design.csv"))
fwrite(rbindlist(rhs_stats),file.path(args[4],"outcomes.csv"))
fwrite(data.table(row=rep(1:6,6L),column=rep(1:6,each=6L),value=as.numeric(meat)),
       file.path(args[4],"joint_score_products.csv"))
fwrite(data.table(target=c("P90","P95"),linear=contrast_basis[,1],cubic=contrast_basis[,2]),
       file.path(args[4],"temperature_contrasts.csv"))
