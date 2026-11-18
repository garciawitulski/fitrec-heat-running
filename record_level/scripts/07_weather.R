# Tables S6, S22–S23 and S26–S27. Each pair of models uses the same observations.
# The standard error of the difference accounts for the runners shared by both models.

local({
  weather_dir <- "covariates"
  cfg <- jsonlite::fromJSON(input_file("settings/session_model.json"))
  main <- read_fitrec("panels/sessions.parquet")[is.finite(day_max_c)]
  phys <- read_fitrec("panels/paired_heart_rate.parquet")
  panel <- read_fitrec("panels/runner_days.parquet")[is.finite(day_max_c)]
  panel[, date_key := format(as.IDate(local_date), "%Y-%m-%d")]
  anchors_day <- quantile(main$day_max_c, c(.5,.9,.95), names = FALSE)
  anchors_session <- quantile(phys$temp_c_10km, c(.5,.9,.95), names = FALSE)
  ids <- c("b1_heat_avoided","b2_speed","b3_log_duration","h1_hr_mean",
           "s_start_minutes","s_start_before_09","s_start_after_19")
  variables <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean",
                 "start_minutes_local","start_before_09","start_after_19")
  jobs <- data.table(id = ids, variable = variables, cohort = "main")
  jobs <- rbind(jobs, data.table(id = c("p_run","h2_pooled_session_temp","ov_segment_hr_daymax_fixed_speed"),
    variable = c("run","hr_lag60","hr_lag60"), cohort = c("panel","physio","physio")))
  estimates <- coefficients <- joint <- list()
  for (check in c("daylight","weather","leads")) {
    for (population in c("main","panel","physio")) {
      if ((check == "leads" && population == "panel") || (check == "daylight" && population == "physio")) next
      d <- switch(population, main = main, panel = panel, physio = phys)
      key <- if (population == "panel") c("user_hash","date_key") else "session_hash"
      prefix <- if (check == "weather") "weather_" else "daylight_"
      extras <- read_fitrec(file.path(weather_dir,paste0(prefix,population,".parquet")))
      if (population == "panel") extras[, date_key := as.character(date_key)]
      controls <- switch(check, daylight = c("sunrise_local_min","sunset_local_min"),
        leads = c("lead1_c","lead2_c","lead3_c"),
        weather = if (population == "physio") c("a3_s_wind_ms","a3_s_rad_wm2") else
          c("a3_precip_mm","a3_wind_ms","a3_rad_mj"))
      wanted <- c(key,controls,if (check == "leads") "at_usual_cell")
      stopifnot(!anyDuplicated(d, by = key), !anyDuplicated(extras, by = key))
      d <- merge(d, extras[, ..wanted], by = key, all.x = TRUE, sort = FALSE)
      if (check == "daylight") d[, (controls) := lapply(.SD, "/", 60), .SDcols = controls]
      if (check == "leads") d[, (controls) := lapply(.SD, "/", 5), .SDcols = controls]
      d <- d[complete.cases(d[, ..controls])]
      outcome_models <- jobs[cohort == population]
      if (check == "leads") outcome_models <- outcome_models[!id %in% ids[6:7]]
      if (check != "leads") outcome_models <- outcome_models[id != "ov_segment_hr_daymax_fixed_speed"]
      samples <- if (check == "leads") c("all with three leads","at the usual cell") else "all"
      for (j in seq_len(nrow(outcome_models))) for (s in samples) {
        model_id <- outcome_models$id[j]
        response <- outcome_models$variable[j]
        session_temp <- model_id == "h2_pooled_session_temp"
        temp <- if (session_temp) "temp_c_10km" else "day_max_c"
        a <- if (session_temp) anchors_session else anchors_day
        knots <- if (session_temp) cfg$knots$session_temp_rcs3 else cfg$knots$rcs3
        z <- d[is.finite(get(response)) & is.finite(get(temp))]
        if (s == "at the usual cell") z <- z[as.logical(at_usual_cell)]
        x <- rcs_basis(z[[temp]],knots)
        if (population != "panel")
          x <- cbind(x,season_terms(z$canonical_start_utc,z$geo_mean_y < 0))
        if (population == "physio") {
          speed <- rcs_basis(z$speed_lag60,cfg$speed_knots,center = 0,scale = 1)
          colnames(speed) <- c("speed1","speed2")
          x <- cbind(x,speed,prior = z$prior_recorded_minutes,fitness = z$fitness_speed_kmh_28d)
        } else {
          x <- cbind(x,model.matrix(~factor(local_weekday,levels = 0:6),z)[,-1,drop=FALSE])
        }
        groups <- if (population == "panel") paste(z$user_hash,z$year_month_local) else z$user_hash
        f0 <- fit_fe(as.numeric(z[[response]]),x,groups,z$user_hash)
        f1 <- fit_fe(as.numeric(z[[response]]),cbind(x,as.matrix(z[, ..controls])),groups,z$user_hash)
        stopifnot(identical(rownames(f0$score),rownames(f1$score)))
        basis <- rcs_basis(a,knots)
        for (k in 2:3) {
          contrast <- basis[k,] - basis[1L,]
          b0 <- lin_contrast(f0,contrast)
          b1 <- lin_contrast(f1,contrast)
          score <- (f1$score[,names(contrast)] - f0$score[,names(contrast)]) %*% contrast
          diff <- b1$estimate - b0$estimate
          se <- sqrt(sum(score^2))
          radius <- qt(.975,f0$df) * se
          delta <- data.table(estimate = diff, lower = diff-radius, upper = diff+radius,
                              p = 2*pt(-abs(diff/se),f0$df))
          estimates[[length(estimates)+1L]] <- cbind(data.table(check, sample=s,
            base_model=model_id,target_label=c("P90","P95")[k-1L],
            version=c("without","with","with_minus_without")),rbind(b0,b1,delta))
        }
        for (k in seq_along(controls)) {
          term <- controls[k]
          coefficients[[length(coefficients)+1L]] <- cbind(data.table(check, sample=s,
            base_model=model_id, term_no=k),lin_contrast(f1,setNames(1,term)))
        }
        if (check == "leads") {
          b <- f1$beta[controls]
          stat <- as.numeric(crossprod(b,solve(f1$vcov[controls,controls],b)))/length(b)
          joint[[length(joint)+1L]] <- data.table(sample=s,base_model=model_id,
            F_statistic=stat,p_value=pf(stat,length(b),f1$df,lower.tail=FALSE))
        }
      }
    }
  }
  fwrite(rbindlist(estimates),file.path(out,"weather_contrasts.csv"))
  fwrite(rbindlist(coefficients),file.path(out,"weather_coefficients.csv"))
  fwrite(rbindlist(joint),file.path(out,"lead_joint.csv"))
})
