# Compare model specifications using the design matrices and samples from script 02.
save_contrast <- function(fit, c, transform, table, model, outcome, target) {
  cbind(data.table(table=table,model=model,outcome=outcome,target=target),
        lin_contrast(fit,c,transform))
}
runner <- sessions$user_hash
ym <- paste(runner,sessions$year_month_local)
year <- factor(format(as.POSIXct(sessions$canonical_start_utc,origin="1970-01-01",tz="UTC"),
                      "%Y",tz="UTC"))
year_indicators <- model.matrix(~year)[,-1,drop=FALSE]
d95 <- target_basis[3,]-target_basis[1,]
results_alt <- list()

# Table S2. The first nested model has a common intercept and no runner fixed effects.
for (k in 1:4) {
  X <- switch(k,cbind(intercept=1,temperature_basis),
              cbind(intercept=1,X_session),X_session,cbind(X_session,year_indicators))
  fe <- if (k>=3) runner else NULL
  for (j in 1:5) {
    v <- outcomes$outcome[j]
    mod <- fit_fe(sessions[[v]],X,fe,runner)
    results_alt[[length(results_alt)+1L]] <- save_contrast(mod,d95,outcomes$transform[j],
                                                           "S2",as.character(k),v,"P95")
  }
}

# Table S25. Each row changes one component of the main model from Table 2.
specs <- c("main","daily_mean","era5","20km","two_way","runner_month","dew",
           "duration_20_180","at_least_20","hinge","four_knots","history","levels_time",
           "levels_distance","equal_runner","daily_mean_own")
n_runner <- sessions[, .N,by=user_hash][sessions,on="user_hash",N]
for (spec in specs) {
  temp <- sessions$day_max_c
  knots <- temperature_knots
  targets <- temperature_targets
  fe <- runner
  weights <- rep(1,nrow(sessions))
  controls <- cbind(season,weekday)
  dates <- NULL
  select <- rep(TRUE,nrow(sessions))
  if (spec %in% c("daily_mean","daily_mean_own")) temp <- sessions$day_mean_c
  if (spec=="daily_mean_own") targets <- as.numeric(quantile(temp,c(.5,.9,.95),na.rm=TRUE))
  if (spec=="era5") temp <- sessions$day_max_era5_c
  if (spec=="20km") temp <- sessions$day_max_c_20km
  if (spec=="two_way") dates <- sessions$local_date
  if (spec=="runner_month") { fe <- ym; controls <- weekday }
  if (spec=="dew") controls <- cbind(controls,dew=sessions$day_max_dew_c)
  if (spec=="duration_20_180") select <- between(sessions$duration_minutes,20,180)
  if (spec=="at_least_20") select <- n_runner>=20
  if (spec=="four_knots") knots <- as.numeric(quantile(temp,c(.05,.35,.65,.95)))
  if (spec=="history") {
    controls <- cbind(controls,minutes_prior=sessions$prior_recorded_minutes)
    select <- is.finite(sessions$prior_recorded_minutes)
  }
  if (spec=="equal_runner") weights <- 1/n_runner
  basis <- rcs_basis(temp,knots)
  b_targets <- rcs_basis(targets,knots)
  if (spec=="hinge") {
    basis <- cbind(hot=pmax(temp-20,0))
    b_targets <- cbind(hot=pmax(targets-20,0))
  }
  X <- cbind(basis,controls)
  vars <- if (spec=="levels_time") "duration_minutes" else if (spec=="levels_distance") "gps_distance_km" else outcomes$outcome[1:4]
  for (v in vars) {
    y <- sessions[[v]]
    # The ERA5 sensitivity analysis replaces the exposure and retains the observed thermal gap.
    if (spec=="equal_runner") weights <- 1/ave(as.numeric(is.finite(y)),runner,FUN=sum)
    u <- if (v=="log_duration") "percent" else "identity"
    mod <- fit_fe(y[select],X[select,,drop=FALSE],fe[select],runner[select],
                  weights[select],if(!is.null(dates)) dates[select])
    if (v==vars[1L]) {
      # The table reports behavioural sample sizes. Apply the disclosure rule
      # that rounds these counts down to the nearest ten.
      results_alt[[length(results_alt)+1L]] <- data.table(table="S25",model=spec,outcome="sample",
          target="sample",n_sessions=10*floor(mod$n/10),n_runners=10*floor(mod$groups/10))
    }
    for (i in 2:3) {
      results_alt[[length(results_alt)+1L]] <- save_contrast(mod,b_targets[i,]-b_targets[1,],u,
                                           "S25",spec,v,c("P90","P95")[i-1L])
    }
  }
}

# Table S28. Walking, pauses and variability in running speed.
comp <- read_fitrec("panels/composition.parquet",
                    c("session_hash","median_interval_s","walk_share_20_pct","walk_share_15_pct",
                      "pause_share_pct","speed_cv_60s"))
comp <- comp[sessions[,.(session_hash)],on="session_hash"]
stopifnot(nrow(comp)==nrow(sessions))
for (v in c("walk_share_20_pct","walk_share_15_pct","pause_share_pct","speed_cv_60s")) {
  for (spec in c("main","runner_month","sampling15")) {
    select <- if (spec=="sampling15") is.finite(comp$median_interval_s) & comp$median_interval_s<=15 else rep(TRUE,nrow(comp))
    X <- if (spec=="runner_month") cbind(temperature_basis,weekday) else X_session
    fe <- if (spec=="runner_month") ym else runner
    mod <- fit_fe(comp[[v]][select],X[select,,drop=FALSE],fe[select],runner[select])
    results_alt[[length(results_alt)+1L]] <- data.table(table="S28",model=spec,outcome=v,
          target="sample",n_sessions=10*floor(mod$n/10),n_runners=10*floor(mod$groups/10))
    for (i in 2:3) {
      result <- save_contrast(mod,target_basis[i,]-target_basis[1,],"identity",
                              "S28",spec,v,c("P90","P95")[i-1L])
      result[,mean_outcome:=mean(comp[[v]][select][mod$keep])]
      results_alt[[length(results_alt)+1L]] <- result
    }
  }
}
supp_session <- rbindlist(results_alt,fill=TRUE)
fwrite(supp_session,file.path(out,"supp_session.csv"))
