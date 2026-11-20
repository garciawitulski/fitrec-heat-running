# Tables S9, S10 and S23. Sessions receive equal weight and errors are clustered by runner.
# The reference temperature is P50 in models of ambient heat history and 20 °C
# in models of prior running time in heat. Each model specifies its cohort and controls.
local({
  weighted_dir <- "history/weighted"
  reference_dir <- "references/additional"
  config <- jsonlite::fromJSON(input_file("settings/session_model.json"))
  fields <- c("estimate","se","ci_low","ci_high","p_value")
  labels <- c("P50","P90","P95")
  weather <- read_fitrec(file.path(weighted_dir,"behaviour_daily_weighted.parquet"),
                         c("user_hash","history_anchor_date","daily_max_c"))
  weather[,date:=as.IDate(history_anchor_date)]
  setkey(weather,user_hash,date)
  stopifnot(!anyDuplicated(weather,by=c("user_hash","date")))
  sessions <- read_fitrec("panels/sessions.parquet",
                          c("session_hash","year_month_local"))
  phys <- read_fitrec("panels/paired_heart_rate.parquet")
  stopifnot(!anyDuplicated(phys$session_hash),!anyDuplicated(sessions$session_hash))
  setorder(phys,user_hash,local_date,session_hash)
  temp <- rcs_basis(phys$temp_c_10km, config$knots$session_temp_rcs3)
  speed <- rcs_basis(phys$speed_lag60, config$speed_knots, center=0, scale=1)
  colnames(temp) <- c("D1","D2")
  colnames(speed) <- c("speed1","speed2")
  xphys <- cbind(temp, speed, season_terms(phys$canonical_start_utc,phys$geo_mean_y<0),
                 minutes=phys$prior_recorded_minutes, fitness_speed=phys$fitness_speed_kmh_28d)
  fphys <- fit_fe(phys$hr_lag60,xphys,phys$user_hash,phys$user_hash)
  anchors_phys <- setNames(as.numeric(quantile(phys$temp_c_10km,c(.5,.9,.95))),labels)
  # For each session, average the runner's speed-adjusted heart rate during cool sessions
  # 28 to 42 days earlier, including both endpoints.
  phys[,speed_residual:=hr_lag60-drop(speed %*% fphys$beta[colnames(speed)])]
  phys[,cool:=temp_c_10km<anchors_phys["P50"]]
  phys[,fitness_index:= {
    day <- as.integer(as.IDate(local_date))
    before <- findInterval(day-43L,day)+1L
    through <- findInterval(day-28L,day)+1L
    total <- c(0,cumsum(speed_residual*cool))
    count <- c(0,cumsum(cool))
    (total[through]-total[before])/(count[through]-count[before])
  },by=user_hash]

  hist_results <- list()
  # All Table S9 models use the same contrast. The alternatives vary the heat-history
  # terms, fixed effects or sample of sessions.
  history_rows <- function(d,controls,terms,check,groups=d$user_hash) {
    fit <- fit_fe(d$y,as.matrix(d[,..controls]),groups,d$user_hash)
    result <- list()
    for (term in terms) {
      for (label in labels[-1L]) {
        delta <- as.numeric(diff(rcs_basis(c(anchors["P50"],anchors[label]),knots)))
        row <- contrast_panel(fit,setNames(delta,paste0(term,c(":D1",":D2"))))
        result[[length(result)+1L]] <- row[, `:=`(history_term=term,quantity="modification_per_5C",target_label=label)]
      }
      row <- contrast_panel(fit,setNames(1,term))
      result[[length(result)+1L]] <- row[, `:=`(history_term=term,quantity="level_per_5C_at_current_20C",target_label="")]
    }
    value <- rbindlist(result)
    value[, `:=`(history_model=mid,check=check)]
    value
  }
  for (mid in c("EW14_H2","EW14_H1")) {
    steady <- mid=="EW14_H2"
    knots <- if(steady) config$knots$session_temp_rcs3 else config$knots$rcs3
    frame <- read_fitrec(file.path(weighted_dir,paste0("fit_design_",mid,".parquet")))
    idx <- fread(input_file(weighted_dir,paste0("design_rows_",mid,".csv")),
                 select=c("session_hash","user_hash","history_anchor_date","history14_exp7_c"))
    stopifnot(identical(frame$user_hash,idx$user_hash),!anyDuplicated(idx$session_hash),
              max(abs(20+5*frame$H-idx$history14_exp7_c))<1e-8)
    cols <- setdiff(names(frame),c("user_hash","y",grep("^year_",names(frame),value=TRUE)))
    controls <- setdiff(cols,c("H","H:D1","H:D2"))
    anchors <- setNames(as.numeric(quantile(20+10*read_fitrec(file.path(weighted_dir,
      paste0("fit_design_",if(steady) "EW14_H2" else "EW14_B1",".parquet")),"D1")$D1,c(.5,.9,.95))),labels)
    frame[, `:=`(session_hash=idx$session_hash,date=as.IDate(idx$history_anchor_date))]
    frame[sessions,on="session_hash",month:=i.year_month_local]
    stopifnot(!anyNA(frame$month))
    hist_results[[length(hist_results)+1L]] <- history_rows(frame,cols,"H","Main specification")
    hist_results[[length(hist_results)+1L]] <- history_rows(frame,cols[!startsWith(cols,"season_")],"H",
            "Runner-by-year-month effects",groups=paste(frame$user_hash,frame$month))
    # Reconstruct the four weekly means from the daily series. Require a temperature
    # on each of the preceding 28 days, including days without a run.
    weeks <- sapply(1:28,function(k) weather[.(frame$user_hash,frame$date-k),daily_max_c])
    # The original two-week values use the usual cell for each subcohort.
    # The four-week analysis uses the common daily series.
    lag_names <- sprintf("lag%02d_daily_max_c",1:14)
    lag_source <- if(steady) "h2_lag_values.parquet" else "behaviour_lag_values.parquet"
    lag_table <- read_fitrec(file.path(weighted_dir,lag_source),c("session_hash",lag_names))
    stopifnot(!anyDuplicated(lag_table$session_hash),all(frame$session_hash %in% lag_table$session_hash))
    lags <- as.matrix(lag_table[match(frame$session_hash,lag_table$session_hash),..lag_names])
    w <- 2^(-(1:14)/7)
    recovered <- drop(lags[,1:14,drop=FALSE] %*% w)/sum(w)
    stopifnot(all(is.finite(recovered)),max(abs(recovered-idx$history14_exp7_c))<1e-6)
    frame[,L1:=(lags[,1]-20)/5]
    frame[,L2_14:=(drop(lags[,2:14,drop=FALSE] %*% w[2:14])/sum(w[2:14])-20)/5]
    for(j in 1:4) set(frame,j=paste0("W",j),value=(rowMeans(weeks[,7*(j-1L)+(1:7),drop=FALSE])-20)/5)
    groups <- list(c("L1","L2_14"),paste0("W",1:4))
    multiple <- c("Lag 1 and lags 2-14","Weekly means over lags 1-28")
    single <- c("Weighted history on 14-day sample","Weighted history on 28-day sample")
    for(i in 1:2) {
      terms <- groups[[i]]
      interactions <- character()
      for(term in terms) for(v in c("D1","D2")) {
        name <- paste(term,v,sep=":")
        set(frame,j=name,value=frame[[term]]*frame[[v]])
        interactions <- c(interactions,name)
      }
      sample <- frame[complete.cases(frame[,..terms])]
      hist_results[[length(hist_results)+1L]] <- history_rows(sample,c(controls,terms,interactions),terms,multiple[i])
      hist_results[[length(hist_results)+1L]] <- history_rows(sample,cols,"H",single[i])
    }
    if(steady) {
      frame[phys,on="session_hash",fitness_index:=i.fitness_index]
      sample <- frame[is.finite(fitness_index)]
      hist_results[[length(hist_results)+1L]] <- history_rows(sample,cols,"H","Without fitness index")
      hist_results[[length(hist_results)+1L]] <- history_rows(sample,c(cols,"fitness_index"),"H","With fitness index")
    }
  }
  results <- rbindlist(hist_results)
  fwrite(results,file.path(out,"06b_history_sensitivity.csv"))
  expected <- fread(input_file(reference_dir,"history_sensitivity.csv"))
  # Within each sample, the weighted-history comparison shares the group
  # prefix of its corresponding multi-lag or weekly-means specification.
  group <- sub(" .*","",expected$check)
  lag_group <- unique(group[grepl("lag 1 and lags",expected$check)])
  week_group <- unique(group[grepl("weekly means",expected$check)])
  is_weighted <- grepl("weighted history",expected$check)
  expected[,check:=fcase(
    grepl("column",check), "Main specification",
    grepl("runner-by-year-month",check), "Runner-by-year-month effects",
    grepl("lag 1 and lags",check), "Lag 1 and lags 2-14",
    grepl("weekly means",check), "Weekly means over lags 1-28",
    is_weighted & group %in% lag_group, "Weighted history on 14-day sample",
    is_weighted & group %in% week_group, "Weighted history on 28-day sample",
    grepl("without fitness",check), "Without fitness index",
    grepl("with fitness",check), "With fitness index", default=check)]
  stopifnot(setequal(expected$check,results$check))
  compare_aggregate(results,expected,
                    c("history_model","check","history_term","quantity","target_label"),fields,"06b_history_sensitivity")

  # Table S10. Fit both models on sessions with an observed fitness index. To estimate
  # uncertainty in their difference, retain each runner's contribution to both fits.
  ok <- is.finite(phys$fitness_index)
  runner <- phys$user_hash[ok]
  x <- xphys[ok,,drop=FALSE]
  y <- phys$hr_lag60[ok]
  f0 <- fit_fe(y,x,runner,runner)
  f1 <- fit_fe(y,cbind(x,fitness_index=phys$fitness_index[ok]),runner,runner)
  terms <- c("D1","D2")
  shared <- f1$score[rownames(f0$score),terms,drop=FALSE]-f0$score[,terms,drop=FALSE]
  fd <- list(beta=f1$beta[terms]-f0$beta[terms],vcov=crossprod(shared),df=f0$df)
  fits <- list(without=f0,with=f1,with_minus_without=fd)
  rows <- list()
  for(label in labels[-1L]) {
    delta <- as.numeric(diff(rcs_basis(c(anchors_phys["P50"],anchors_phys[label]),config$knots$session_temp_rcs3)))
    for(v in names(fits)) {
      row <- contrast_panel(fits[[v]],setNames(delta,terms))
      rows[[length(rows)+1L]] <- row[, `:=`(target_label=label,version=v,term="")]
    }
  }
  rows[[length(rows)+1L]] <- contrast_panel(f1,c(fitness_index=1))[, `:=`(target_label="",version="",term="fitness_index (per beat)")]
  fitness <- rbindlist(rows)
  fwrite(fitness,file.path(out,"06b_fitness.csv"))
  compare_aggregate(fitness,fread(input_file(reference_dir,"fitness_adjustment.csv")),
                    c("target_label","version","term"),fields,"06b_fitness")

  # In the northern hemisphere, the phases are May–June, July and August–September.
  mon <- as.integer(format(as.IDate(phys$local_date),"%m"))
  mon <- ifelse(phys$geo_mean_y<0,(mon+5)%%12+1,mon)
  group <- fcase(mon %in% 5:6,"early",mon==7,"middle",mon %in% 8:9,"late",default=NA_character_)
  x <- cbind(xphys, middle=as.numeric(group=="middle"),late=as.numeric(group=="late"),
             middle1=(group=="middle")*temp[,1],middle2=(group=="middle")*temp[,2],
             late1=(group=="late")*temp[,1],late2=(group=="late")*temp[,2])
  fs <- fit_fe(phys$hr_lag60,x,phys$user_hash,phys$user_hash)
  season <- list()
  for(label in labels[-1L]) for(g in c("early","middle","late")) {
    delta <- as.numeric(diff(rcs_basis(c(anchors_phys["P50"],anchors_phys[label]),config$knots$session_temp_rcs3)))
    z <- setNames(delta,terms)
    if(g!="early") z <- c(z,setNames(delta,paste0(g,1:2)))
    row <- contrast_panel(fs,z)
    season[[length(season)+1L]] <- row[, `:=`(target_label=label,group_level=g)]
  }
  season <- rbindlist(season)
  fwrite(season,file.path(out,"06b_season.csv"))
  compare_aggregate(season,fread(input_file(reference_dir,"warm_season.csv")),
                    c("target_label","group_level"),fields,"06b_season")
  restrictions <- c("middle1","middle2","late1","late2")
  b <- fs$beta[restrictions]
  f <- drop(crossprod(b,solve(fs$vcov[restrictions,restrictions],b)))/length(b)
  test <- data.table(test="exposure_spline_by_group_interaction",F_statistic=f,p_value=pf(f,4,fs$df,lower.tail=FALSE))
  fwrite(test,file.path(out,"06b_season_joint.csv"))
  compare_aggregate(test,as.data.table(jsonlite::fromJSON(input_file(reference_dir,"warm_season_test.json"))),
                    "test",c("F_statistic","p_value"),"06b_season_joint")

  # Table S23. Define exposure tertiles across all sessions eligible for that exposure,
  # then fit the model in the steady-speed subcohort.
  prior <- read_fitrec("cohort/session_weather.parquet",
                       c("session_hash","user_hash","minutes_ge25","prior_history_complete_A","temp_c_10km"))
  eligible <- read_fitrec("cohort/geolocated_workouts.parquet",
                          c("session_hash","outcome_eligible_B"))
  prior <- prior[eligible[outcome_eligible_B==TRUE],on="session_hash",nomatch=0L]
  prior <- prior[is.finite(temp_c_10km) & prior_history_complete_A==TRUE & is.finite(minutes_ge25)]
  prior[, `:=`(q1=as.numeric(quantile(minutes_ge25,1/3)),q2=as.numeric(quantile(minutes_ge25,2/3))),by=user_hash]
  prior[,tertile:=fcase(minutes_ge25<=q1,1L,minutes_ge25>q2,3L,default=2L)]
  phys[prior,on="session_hash", `:=`(tertile=i.tertile,dose=i.minutes_ge25)]
  stopifnot(!anyNA(phys$tertile))
  middle <- as.numeric(phys$tertile==2L)
  high <- as.numeric(phys$tertile==3L)
  x <- cbind(xphys,middle=middle,high=high,M1=middle*temp[,1],M2=middle*temp[,2],
             H1=high*temp[,1],H2=high*temp[,2])
  fmain <- fit_fe(phys$hr_lag60,x,phys$user_hash,phys$user_hash)
  results <- list()
  for(t in c(20,25,28)) {
    delta <- as.numeric(diff(rcs_basis(c(20,t),config$knots$session_temp_rcs3)))
    low <- setNames(delta,terms)
    hi <- setNames(delta,c("H1","H2"))
    z <- list(R_low=low,R_middle=c(low,setNames(delta,c("M1","M2"))),R_high=c(low,hi),C=hi,L=c(high=1,hi))
    for(name in names(z)) {
      row <- contrast_panel(fmain,z[[name]])
      results[[length(results)+1L]] <- row[, `:=`(temperature_c=t,contrast=name)]
    }
  }
  recorded <- rbindlist(results)
  fwrite(recorded,file.path(out,"06b_prior_recorded_heat.csv"))
  ref <- fread(input_file("references/heart_rate_and_history/recorded_heat_contrasts.csv"))[model_id=="main"]
  compare_aggregate(recorded,ref,c("temperature_c","contrast"),fields,"06b_prior_recorded_heat")
  desc <- data.table(sessions=nrow(phys),runners=uniqueN(phys$user_hash),low_pct=100*mean(phys$tertile==1),
                     no_prior_hot_running_pct=100*mean(phys$dose==0),
                     median_positive_minutes=median(phys$dose[phys$dose>0]),
                     mde_80pct=(qnorm(.975)+qnorm(.8))*recorded[temperature_c==28 & contrast=="C",se])
  fwrite(desc,file.path(out,"06b_prior_recorded_counts.csv"))

  # The pooled-cohort sensitivity model adds calendar-year indicators.
  # Appendix B.4.3 reports the difference at 28 °C.
  year <- as.integer(format(as.POSIXct(phys$canonical_start_utc,origin="1970-01-01",tz="UTC"),"%Y"))
  yr <- model.matrix(~factor(year))[, -1L,drop=FALSE]
  fyear <- fit_fe(phys$hr_lag60,cbind(x,yr),phys$user_hash,phys$user_hash)
  delta <- as.numeric(diff(rcs_basis(c(20,28),config$knots$session_temp_rcs3)))
  year_adjusted <- contrast_panel(fyear,setNames(delta,c("H1","H2")))[, `:=`(model_id="main",temperature_c=28,contrast="C")]
  fwrite(year_adjusted,file.path(out,"06b_prior_registered.csv"))
  ref <- read_fitrec("history/prior_exercise/contrasts_full_precision.parquet")
  compare_aggregate(year_adjusted,ref[model_id=="main" & temperature_c==28 & contrast=="C"],
                    c("model_id","temperature_c","contrast"),fields,"06b_prior_registered")
  # The holdout model includes a linear temperature term above 20 °C and a linear
  # speed term. It omits training and seasonal controls.
  held <- read_fitrec("cohort/holdout_outcomes.parquet",
                      c("session_hash","user_hash","outcome_status","hr_lag60","speed_lag60"))
  held <- held[outcome_status=="OK" & is.finite(hr_lag60) & is.finite(speed_lag60)]
  held <- held[prior,on=c("user_hash","session_hash"),nomatch=0L]
  g <- pmax(held$temp_c_10km-20,0)/8
  middle <- as.numeric(held$tertile==2L)
  high <- as.numeric(held$tertile==3L)
  xx <- cbind(g=g,middle=middle,high=high,g_middle=g*middle,g_high=g*high,speed=held$speed_lag60)
  fh <- fit_fe(held$hr_lag60,xx,held$user_hash,held$user_hash)
  hold <- contrast_panel(fh,c(g_high=1))[, `:=`(contrast="g_high",sessions=fh$n)]
  fwrite(hold,file.path(out,"06b_prior_holdout.csv"))
  stored <- jsonlite::fromJSON(input_file("cohort/holdout_estimates.json"))$results$primary
  coef <- stored$coefficients$g_high
  ref <- data.table(contrast="g_high",estimate=coef$estimate,se=coef$se,ci_low=coef$ci[1],ci_high=coef$ci[2],
                    p_value=coef$p_value,sessions=stored$sessions)
  compare_aggregate(hold,ref,"contrast",c(fields,"sessions"),"06b_prior_holdout")
})
