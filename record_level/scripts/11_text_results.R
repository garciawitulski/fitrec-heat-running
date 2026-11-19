# Check the quantities printed in the abstract and results.
local({
  printed_text <- fread(file.path(repo,"data","processed","reported_text.csv"),
                        colClasses=c(printed="character"))
  check <- function(section_name,quantity_name,calculated,exact=FALSE) {
    value <- printed_text[section==section_name & quantity==quantity_name,printed]
    stopifnot(length(value)==1L,length(calculated)==1L,is.finite(calculated))
    digits <- if (grepl(".",value,fixed=TRUE)) nchar(sub(".*\\.","",value)) else 0L
    tolerance <- if (exact) 0 else .5*10^-digits
    difference <- abs(calculated-as.numeric(value))
    data.table(section=section_name,quantity=quantity_name,printed=value,R=as.numeric(calculated),
               paper=as.numeric(value),absolute_difference=difference,tolerance=tolerance,
               status=if(difference<=tolerance+1e-12) "MATCH" else "MISMATCH")
  }
  s <- read_fitrec("panels/sessions.parquet",
            c("session_hash","user_hash","day_max_c","day_min_c","hr_mean","gps_speed_median_m_s"))
  s <- s[is.finite(day_max_c)]
  subcohort <- read_fitrec("panels/steady_speed.parquet","session_hash")
  targets <- quantile(s$day_max_c,c(.50,.95),names=FALSE)
  pair <- fread(file.path(out,"heart_paired.csv"))[target_label=="P95"]
  diff <- pair[quantity=="steady_speed_free_minus_whole_session"]
  main <- fread(file.path(out,"table2.csv"))[target=="P95"]
  parts <- fread(file.path(out,"heart_decomposition.csv"))[target_label=="P95"]
  percent <- 100*parts[quantity %in% c("timing_offset_fraction_fixed_hour","speed_offset_fraction_fixed_hour"),sum(estimate)]
  bike <- fread(file.path(out,"04_sport.csv"))[outcome=="outdoor_cycling" & variant=="c3" & target_label=="P95",100*estimate]
  # Use main-cohort runs with admissible heart-rate readings to calculate artefact frequencies.
  tw <- rbindlist(lapply(1:4,function(i) read_fitrec(sprintf(
    "panels/time_weighted_hr_part%d.parquet",i),
    c("session_hash","spikes_removed","flat_readings_removed"))))
  screened <- tw[s[is.finite(hr_mean),.(session_hash)],on="session_hash"]
  stopifnot(!anyDuplicated(tw$session_hash),!anyNA(screened))
  checks <- rbindlist(list(
    check("Abstract","sessions",nrow(s),TRUE),
    check("Abstract","runners",uniqueN(s$user_hash),TRUE),
    check("Abstract","P50 temperature",targets[1]),
    check("Abstract","P95 temperature",targets[2]),
    check("Abstract","paired sessions",sum(s$session_hash %in% subcohort$session_hash),TRUE),
    check("Abstract","steady HR increase",pair[quantity=="steady_speed_free",estimate]),
    check("Abstract","whole-session HR increase",pair[quantity=="whole_session_mean",estimate]),
    check("Abstract","paired HR difference",diff$estimate),
    check("Abstract","paired HR difference lower",diff$lower),
    check("Abstract","paired HR difference upper",diff$upper),
    check("Abstract","earlier start (minutes)",-main[outcome=="start_minutes_local",estimate]),
    check("Abstract","slower than mean speed (%)",
          -100*main[outcome=="gps_speed_median_m_s",estimate]/mean(s$gps_speed_median_m_s)),
    check("Results","timing and speed offset (%)",percent),
    check("Abstract","outdoor cycling (pp)",bike),
    check("Results","temperature range near P95",
          s[abs(day_max_c-targets[2])<=1,mean(day_max_c-day_min_c,na.rm=TRUE)]),
    check("Results","temperature range near P50",
          s[abs(day_max_c-targets[1])<=1,mean(day_max_c-day_min_c,na.rm=TRUE)]),
    check("Results","sessions with isolated spikes (%)",100*mean(screened$spikes_removed>0)),
    check("Results","sessions with identical sequences (%)",100*mean(screened$flat_readings_removed>0))
  ))
  fwrite(checks,file.path(out,"check_11_text.csv"))
  print(checks[,.N,by=.(section,status)])
})
