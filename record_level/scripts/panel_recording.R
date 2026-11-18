# Table S4 and Figure S3 use the calendar-day panel and workout records prepared in
# 04_participation.R. Workouts are the observations for heart-rate recording, and
# runner-days are the observations for participation.
local({
  dir_ref <- "references/sensitivity"
  daily <- days[,.(user_hash,date,day_max_c,local_weekday,south,time,month,run,run_any)]
  info <- meta[eligible_recorded_running_history==TRUE,.(user_hash,date,has_hr_record)]
  runs <- daily[info,on=c("user_hash","date"),nomatch=0L]
  runs[,hr_recorded:=as.numeric(has_hr_record)]
  # Estimate all three curves over the full temperature grid. The table also
  # reports the recording-availability model with runner-by-month fixed effects.
  models <- c(p_run="run",p_run_any_recorded="run_any",c_hr_recorded="hr_recorded",
              c_hr_recorded_runner_month="hr_recorded")
  grid <- fread(input_file(dir_ref,"hr_recording_curves.csv"),select=c("model_id","exposure_c"))
  curve <- rows <- list()
  for(mid in names(models)) {
    d <- if(startsWith(mid,"p_")) daily else runs
    variant <- if(mid=="c_hr_recorded") "main" else "month"
    x <- rcs_basis(d$day_max_c,knots)
    colnames(x) <- c("T1","T2")
    f <- model_panel(d,models[[mid]],x,3L,variant)
    xs <- grid[model_id==mid,exposure_c]
    for(t in unique(c(xs,targets[c("P90","P95")]))) {
      delta <- as.numeric(diff(rcs_basis(c(targets["P50"],t),knots)))
      if(abs(t-targets["P50"])<1e-10) delta[] <- 0
      row <- contrast_panel(f,setNames(delta,c("T1","T2")))
      row[, `:=`(model_id=mid,exposure_c=t)]
      if(t %in% xs) curve[[length(curve)+1L]] <- row
      if(t %in% targets[c("P90","P95")]) {
        row <- copy(row)
        rows[[length(rows)+1L]] <- row[,target_label:=names(targets)[match(t,targets)]]
      }
    }
  }
  fields <- c("estimate","se","ci_low","ci_high","p_value")
  curve <- rbindlist(curve)
  rows <- rbindlist(rows)
  fwrite(curve,file.path(out,"04_recording_curves.csv"))
  fwrite(rows,file.path(out,"04_recording_contrasts.csv"))
  compare_aggregate(curve,fread(input_file(dir_ref,"hr_recording_curves.csv")),
                    c("model_id","exposure_c"),fields,"04_recording_curves")
  compare_aggregate(rows,fread(input_file(dir_ref,"hr_recording_contrasts.csv")),
                    c("model_id","target_label"),fields,"04_recording_contrasts")
  desc <- c(share_of_recorded_runs_with_heart_rate=mean(runs$hr_recorded),
            share_with_heart_rate_within_2C_of_P50=mean(runs[abs(day_max_c-targets["P50"])<=2,hr_recorded]),
            share_with_heart_rate_within_2C_of_P95=mean(runs[abs(day_max_c-targets["P95"])<=2,hr_recorded]),
            runner_days_with_a_heart_rate_run_share=mean(daily$run),runner_days_with_any_recorded_run_share=mean(daily$run_any))
  ref <- jsonlite::fromJSON(input_file(dir_ref,"hr_recording_descriptives.json"))
  stopifnot(paste0(">=", 10*floor(nrow(runs)/10))==ref$recorded_runs_in_panel_windows)
  compare_aggregate(data.table(measure=names(desc),value=unname(desc)),
                    data.table(measure=names(desc),value=as.numeric(unlist(ref[names(desc)]))),"measure","value","04_recording_shares")
})
