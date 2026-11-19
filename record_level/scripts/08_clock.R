# Table S16. First, measure start time relative to sunrise.
# Then compare sessions before and after a clock change in the time zone of the run.
local({
  dir <- "covariates"
  main <- read_fitrec("panels/sessions.parquet")[is.finite(day_max_c)]
  anchors <- quantile(main$day_max_c,c(.5,.9,.95),names=FALSE)
  knots <- quantile(main$day_max_c,c(.1,.5,.9),names=FALSE)
  extra <- read_fitrec(file.path(dir,"daylight_main.parquet"),
    c("session_hash","sunrise_local_min","sunset_local_min"))
  stopifnot(!anyDuplicated(main$session_hash),!anyDuplicated(extra$session_hash))
  d <- merge(main,extra,by="session_hash",all.x=TRUE,sort=FALSE)
  d[, start_minus_sunrise := start_minutes_local - sunrise_local_min]
  sun <- d[is.finite(sunrise_local_min) & is.finite(sunset_local_min)]
  sun[, start_before_sunrise := as.numeric(start_minutes_local < sunrise_local_min)]
  sun[, start_after_sunset := as.numeric(start_minutes_local >= sunset_local_min)]
  wd <- model.matrix(~factor(local_weekday,levels=0:6),sun)[,-1,drop=FALSE]
  X <- cbind(rcs_basis(sun$day_max_c,knots),season_terms(sun$canonical_start_utc,sun$geo_mean_y<0),wd)
  basis <- rcs_basis(anchors,knots)
  rows <- list()
  for (y in c("start_minus_sunrise","start_before_sunrise","start_after_sunset")) {
    fit <- fit_fe(sun[[y]],X,sun$user_hash,sun$user_hash)
    for (k in 2:3) rows[[length(rows)+1L]] <- cbind(data.table(outcome=y,
      target_label=c("P90","P95")[k-1L]),lin_contrast(fit,basis[k,]-basis[1,]))
  }
  fwrite(rbindlist(rows),file.path(out,"sunrise_starts.csv"))

  # Find the first day that the new 60-minute offset is in effect at local noon.
  calendar <- seq(as.Date("2005-10-01"),as.Date("2017-03-31"),by="day")
  noon <- as.POSIXct(paste(calendar,"12:00:00"),tz="UTC")
  changes <- list()
  for (zone in unique(d$gps_iana_zone)) {
    offset <- format(noon,format="%z",tz=zone)
    minutes <- ifelse(substr(offset,1,1) == "-", -1, 1) *
      (60*as.numeric(substr(offset,2,3)) + as.numeric(substr(offset,4,5)))
    delta <- c(0,diff(minutes))
    ok <- which(abs(delta)==60)
    if (length(ok)) changes[[zone]] <- data.table(gps_iana_zone=zone,
      change=as.IDate(calendar[ok]),forward=as.integer(delta[ok] > 0))
  }
  changes <- rbindlist(changes)
  window <- changes[, .(days_from=-21:20),by=.(gps_iana_zone,change,forward)]
  window[, date := change + days_from]
  # When transition windows overlap, use the most recent past transition.
  # If both transitions are still ahead, use the nearer one.
  window[, `:=`(future=days_from<0, distance=abs(days_from))]
  setorder(window,gps_iana_zone,date,future,distance)
  window <- unique(window,by=c("gps_iana_zone","date"))
  d[, date := as.IDate(local_date)]
  z <- merge(d,window,by=c("gps_iana_zone","date"))
  z[, after := as.numeric(days_from>=0)]
  wd <- model.matrix(~factor(local_weekday,levels=0:6),z)[,-1,drop=FALSE]
  X <- cbind(rcs_basis(z$day_max_c,knots),with(z,cbind(
    forward_jump=after*forward,back_jump=after*(1-forward),
    forward_before=days_from*forward*(1-after),forward_after=days_from*forward*after,
    back_before=days_from*(1-forward)*(1-after),back_after=days_from*(1-forward)*after)),wd)
  fe <- paste(z$user_hash,z$gps_iana_zone,z$change)
  contrasts <- list(jump_forward=c(forward_jump=1),jump_back=c(back_jump=1),
    share_of_clock_shift_followed=c(forward_jump=1,back_jump=-1)/120)
  rows <- list()
  for (y in c("start_minutes_local","start_minus_sunrise")) {
    fit <- fit_fe(z[[y]],X,fe,z$user_hash)
    for (key in names(contrasts)) rows[[length(rows)+1L]] <- cbind(
      data.table(outcome=y,quantity=key),lin_contrast(fit,contrasts[[key]]))
  }
  fwrite(rbindlist(rows),file.path(out,"clock_changes.csv"))
})
