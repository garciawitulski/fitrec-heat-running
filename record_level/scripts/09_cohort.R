# Table S1. Apply each eligibility criterion in turn and compare the resulting
# numbers of observations and runners with the published counts for each sample.
local({
  a <- read_fitrec("cohort/geolocated_workouts.parquet",
    c("session_hash","user_hash","all_eligible_HR_history"))[as.logical(all_eligible_HR_history)]
  flags <- read_fitrec("cohort/session_eligibility.parquet",
    c("session_hash","clock_eligible","event_main_keep","duration_minutes","gps_distance_km"))
  times <- read_fitrec("cohort/session_clocks.parquet",c("session_hash","gps_iana_zone"))
  gps <- read_fitrec("cohort/hr_session_audit.parquet",
    c("session_hash","is_running","gps_speed_median_m_s","hr_mean"))
  weather <- read_fitrec("cohort/session_weather.parquet",c("session_hash","temp_c_10km"))
  for (b in list(flags,times,gps,weather)) {
    stopifnot(!anyDuplicated(a$session_hash),!anyDuplicated(b$session_hash))
    a <- merge(a,b,by="session_hash",all.x=TRUE,sort=FALSE)
  }
  linked <- a[is.finite(temp_c_10km)]
  eligible <- linked[as.logical(clock_eligible) & as.logical(event_main_keep) & as.logical(is_running) &
    !is.na(gps_iana_zone) & duration_minutes>0 & gps_distance_km>0 &
    is.finite(hr_mean) & is.finite(gps_speed_median_m_s)]
  sessions <- read_fitrec("panels/sessions.parquet")
  stopifnot(setequal(sessions$session_hash,eligible$session_hash))
  main <- sessions[is.finite(day_max_c)]
  steady <- read_fitrec("panels/paired_heart_rate.parquet",
    c("session_hash","user_hash","day_max_c"))[is.finite(day_max_c)]
  days <- read_fitrec("panels/runner_days.parquet",
    c("user_hash","day_max_c"))[is.finite(day_max_c)]
  stages <- list(confirmed=a,temperature=linked,eligible=eligible,main=main,
    heart_rate=main[is.finite(hr_mean)],steady=steady,history=main[is.finite(prior_recorded_minutes)],panel=days)
  public <- data.table(stage=names(stages),paper_rows=c(115036L,113067L,112797L,106600L,
    106341L,18626L,71656L,503410L),paper_runners=c(862L,860L,850L,820L,810L,660L,801L,729L),
    lower_bound=c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE))
  public[, rows:=vapply(stages,nrow,integer(1))]
  public[, runners:=vapply(stages,function(z) uniqueN(z$user_hash),integer(1))]
  public[lower_bound==TRUE,runners:=10L*(runners%/%10L)]
  public[,status:=ifelse(rows==paper_rows & runners==paper_runners,"MATCH","DIFFERENCE")]
  fwrite(public,file.path(out,"cohort_flow.csv"))
})
