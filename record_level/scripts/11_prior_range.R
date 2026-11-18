# Supplement B.4.3. Within each runner, compare the change in heart rate from
# 20 to 28 °C after high and low heat exposure on earlier runs.
# Start with the panel before the restriction on training covariate availability,
# then use the observations and variables required by each specification.
# Sessions have equal weight. CR1 standard errors are clustered by runner, with t critical values based on the runner count minus one.
local({
  cfg <- jsonlite::fromJSON(input_file("settings/prior_heat_model.json"))
  cols <- c("user_hash","session_hash","temp_c_10km","temp_c_20km","era5_temp_c","prior_recorded_minutes",
            "minutes_ge25","fitness_speed_kmh_28d","prior_history_complete_B","tertile","family_id","geo_mean_y",
            "strict_stable_repeated_route_intersection_B","canonical_start_utc","hr_lag30","speed_lag30",
            "hr_lag60","speed_lag60","hr_lag90","speed_lag90","t_degree","t_decayed","base20","tertile20",
            "base_era5","tertile_era5")
  base <- read_fitrec("history/prior_exercise/base_cohort_check.parquet",cols)
  stopifnot(!anyDuplicated(base$session_hash))
  main <- base[is.finite(prior_recorded_minutes) & is.finite(fitness_speed_kmh_28d)]
  stopifnot(nrow(main)==19350L,uniqueN(main$user_hash)==669L)
  alternatives <- c("main","hinge","piecewise","quadratic","spline4","degree_minutes","decayed",
                    "absolute_dose","cap20","era5","lag30","lag90","complete_history","route_fe","alt_season",
                    "runner_month","interactions","no_training","speed_only")
  temperatures <- c(20,25,28)
  values <- list()
  # Use the same design matrix specification for estimation and prediction.
  # For predictions, hold covariates at their sample means.
  product_columns <- function(a,b) {
    pieces <- lapply(seq_len(ncol(a)),function(j) {
      block <- a[,j]*b
      colnames(block) <- paste(colnames(a)[j],colnames(b),sep=":")
      block
    })
    do.call(cbind,pieces)
  }
  design <- function(t,d,z,mid) {
    u <- (t-20)/8
    temperature <- switch(mid,hinge=cbind(T1=pmax(u,0)),piecewise=cbind(T1=u,T2=pmax(u,0)),
                          quadratic=cbind(T1=u,T2=u^2),spline4=rcs_basis(t,cfg$temperature_knots$F4,20,8),
                          rcs_basis(t,cfg$temperature_knots$F3,20,8))
    colnames(temperature) <- paste0("T",seq_len(ncol(temperature)))
    dose <- if(mid=="absolute_dose") cbind(positive=as.numeric(d>0),hours_log=log1p(d/60)) else
      cbind(middle=as.numeric(d==2),high=as.numeric(d==3))
    X <- cbind(temperature,dose,product_columns(temperature,dose),z)
    if(mid=="interactions") X <- cbind(X,product_columns(temperature,z),product_columns(dose,z))
    X
  }
  for(mid in alternatives) {
    d <- switch(mid,complete_history=main[prior_history_complete_B==TRUE],
                route_fe=main[strict_stable_repeated_route_intersection_B==TRUE],
                cap20=main[base20==TRUE],era5=main[base_era5==TRUE],no_training=base,main)
    exposure <- switch(mid,cap20="temp_c_20km",era5="era5_temp_c","temp_c_10km")
    dose <- switch(mid,degree_minutes="t_degree",decayed="t_decayed",absolute_dose="minutes_ge25",
                   cap20="tertile20",era5="tertile_era5","tertile")
    lag_seconds <- switch(mid,lag30="30",lag90="90","60")
    heart <- paste0("hr_lag",lag_seconds)
    speed <- paste0("speed_lag",lag_seconds)
    d <- d[is.finite(get(exposure)) & is.finite(get(dose)) & is.finite(get(heart)) & is.finite(get(speed))]
    utc <- as.POSIXct(d$canonical_start_utc,origin="1970-01-01",tz="UTC")
    z <- rcs_basis(d[[speed]],cfg$speed_knots,center=cfg$speed_knots[2],scale=1)
    colnames(z) <- c("speed1","speed2")
    if(mid=="alt_season") {
      day <- as.integer(format(utc,"%j",tz="UTC"))
      north <- ifelse(d$geo_mean_y<0,-1,1)
      z <- cbind(z,season=north*cos(2*pi*(day-196)/365.25),
                 fractional_year=as.integer(format(utc,"%Y",tz="UTC"))+(day-1)/365.25)
    } else if(!mid %in% c("speed_only","runner_month")) {
      z <- cbind(z,season_terms(utc,d$geo_mean_y<0))
    }
    if(!mid %in% c("no_training","speed_only")) {
      z <- cbind(z,minutes=d$prior_recorded_minutes,fitness=d$fitness_speed_kmh_28d)
    }
    fe <- switch(mid,runner_month=paste(d$user_hash,format(utc,"%Y-%m",tz="UTC")),
                 route_fe=paste(d$user_hash,d$family_id),d$user_hash)
    X <- design(d[[exposure]],d[[dose]],z,mid)
    fit <- fit_fe(d[[heart]],X,fe,d$user_hash)
    doses <- if(mid=="absolute_dose") c(zero=0,unlist(cfg$dose_values)) else c(low=1,middle=2,high=3)
    zbar <- matrix(colMeans(z),nrow=1,dimnames=list(NULL,colnames(z)))
    prediction <- function(t,dose) drop(design(t,dose,zbar,mid))
    for(t in temperatures) {
      levels <- lapply(doses,function(dose) prediction(t,dose))
      responses <- lapply(doses,function(dose) prediction(t,dose)-prediction(20,dose))
      vectors <- setNames(responses,paste0("R_",names(responses)))
      baseline <- names(doses)[1]
      comparisons <- if(mid=="absolute_dose") names(doses)[-1] else "high"
      for(nm in comparisons) {
        suffix <- if(mid=="absolute_dose") paste0("_",nm) else ""
        vectors[[paste0("C",suffix)]] <- responses[[nm]]-responses[[baseline]]
        vectors[[paste0("L",suffix)]] <- levels[[nm]]-levels[[baseline]]
      }
      for(nm in names(vectors)) {
        row <- contrast_panel(fit,vectors[[nm]])
        values[[length(values)+1L]] <- row[, `:=`(model_id=mid,temperature_c=t,contrast=nm)]
      }
    }
  }
  results <- rbindlist(values)
  fwrite(results,file.path(out,"11_prior_contrasts.csv"))
  source <- "references/heart_rate_and_history"
  ref <- fread(input_file(source,"recorded_heat_contrasts.csv"))
  compare_aggregate(results,ref,c("model_id","temperature_c","contrast"),
                    c("estimate","se","ci_low","ci_high","p_value"),"11_prior_contrasts")
  # The reported range covers the comparable contrasts between tertiles.
  # The continuous-dose model contributes two separate contrasts to the full results.
  alts <- results[model_id!="main" & temperature_c==28 & contrast=="C"]
  range <- data.table(quantity=c("minimum","maximum","alternatives"),
                      value=c(min(alts$estimate),max(alts$estimate),nrow(alts)),
                      model_id=c(alts$model_id[which.min(alts$estimate)],alts$model_id[which.max(alts$estimate)],""))
  fwrite(range,file.path(out,"11_prior_range.csv"))
  ref_range <- jsonlite::fromJSON(input_file(source,"recorded_heat_summary.json"))
  expected <- data.table(quantity=c("minimum","maximum","alternatives"),
                         value=c(ref_range$alternatives_C28_range,ref_range$alternatives))
  compare_aggregate(range,expected,"quantity","value","11_prior_range")
})
