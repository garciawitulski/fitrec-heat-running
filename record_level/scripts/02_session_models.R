# Table 2 and response curves. Each model compares sessions from the same runner
# with measurements of both the exposure and the outcome (Appendix B.1).
sessions <- read_fitrec("panels/sessions.parquet")
sessions <- sessions[is.finite(day_max_c)]
phys_speed <- read_fitrec("panels/steady_speed.parquet",
                          c("session_hash", "speed_lag60"))
sessions[phys_speed, on = "session_hash", steady_speed := i.speed_lag60]
stopifnot(nrow(sessions) == 106600L, sum(is.finite(sessions$hr_mean)) == 106341L,
          sum(is.finite(sessions$steady_speed)) == 18626L)

temperature_knots <- as.numeric(quantile(sessions$day_max_c, c(.10, .50, .90)))
temperature_targets <- as.numeric(quantile(sessions$day_max_c, c(.50, .90, .95)))
temperature_basis <- rcs_basis(sessions$day_max_c, temperature_knots)
target_basis <- rcs_basis(temperature_targets, temperature_knots)
season <- season_terms(sessions$canonical_start_utc, sessions$geo_mean_y < 0)
weekday <- model.matrix(~ factor(local_weekday, levels = 0:6), sessions)[, -1, drop = FALSE]
X_session <- cbind(temperature_basis, season, weekday)
outcomes <- data.table(
  outcome = c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean",
              "start_minutes_local","start_before_09","start_after_19","session_day_percentile",
              "log_distance","steady_speed","run"),
  transform = c("identity","identity","percent","identity","identity","points","points",
                "points","percent","identity","points"))
days <- read_fitrec("panels/runner_days.parquet")
days <- days[is.finite(day_max_c)]
stopifnot(nrow(days) == 503410L)
X_day <- cbind(rcs_basis(days$day_max_c, temperature_knots),
               model.matrix(~factor(local_weekday,levels=0:6),days)[,-1,drop=FALSE])

main_fits <- list()
rows <- list()
curves <- list()
limits <- quantile(sessions$day_max_c,c(.01,.99))
grid <- seq(ceiling(limits[1]*2)/2, floor(limits[2]*2)/2, by=.5)
grid <- sort(unique(c(grid,temperature_targets)))
grid_basis <- rcs_basis(grid,temperature_knots)
for (j in seq_len(nrow(outcomes))) {
  variable <- outcomes$outcome[j]
  units <- outcomes$transform[j]
  if (variable == "run") {
    id <- paste(days$user_hash, days$year_month_local)
    fit <- fit_fe(as.numeric(days$run), X_day, id, days$user_hash)
  } else {
    fit <- fit_fe(as.numeric(sessions[[variable]]), X_session, sessions$user_hash, sessions$user_hash)
  }
  main_fits[[variable]] <- fit
  for (i in 2:3) {
    contrast <- target_basis[i, ] - target_basis[1L, ]
    est <- lin_contrast(fit, contrast, units)
    idx <- names(contrast)
    est[,se_model:=sqrt(drop(crossprod(contrast,fit$vcov[idx,idx] %*% contrast)))]
    rows[[length(rows)+1L]] <- cbind(data.table(outcome=variable, target=c("P90","P95")[i-1L]),est)
  }
  curves[[variable]] <- rbindlist(lapply(seq_along(grid), function(i) {
    diff <- grid_basis[i,]-target_basis[1L,]
    idx <- names(diff)
    est <- lin_contrast(fit,diff,units)
    est[,se_model:=sqrt(drop(crossprod(diff,fit$vcov[idx,idx] %*% diff)))]
    cbind(data.table(outcome=variable,temperature=grid[i]),est)
  }))
}
table2 <- rbindlist(rows)
table2[, p_holm := NA_real_]
table2[target=="P95" & outcome %in% outcomes$outcome[1:4], p_holm := p.adjust(p,"holm")]
fwrite(table2, file.path(out,"table2.csv"))
fwrite(rbindlist(curves), file.path(out,"session_curves.csv"))
