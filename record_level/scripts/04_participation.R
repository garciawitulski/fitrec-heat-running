# Participation and other activities (Figure 6 and Tables S11–S14).
# Give each runner-day or complete ISO runner-week equal weight. Daily contrasts
# compare P90 and P95 with P50. Weekly contrasts replace one day at 10–20 °C
# with one day in another temperature band. Days without workouts count as zero.
local({
  panel_dir <- "panels"
  activity_refs <- "references/additional"
  sessions <- read_fitrec(file.path(panel_dir, "sessions.parquet"),
    c("user_hash", "local_date", "cell_10km", "gps_iana_zone", "day_max_c", "duration_minutes", "gps_distance_km"))
  days <- read_fitrec(file.path(panel_dir, "runner_days.parquet"))[is.finite(day_max_c)]
  days[, date := as.IDate(local_date)]
  days[, time := as.numeric(noon_utc)]
  stopifnot(nrow(days) == 503410L, uniqueN(days$user_hash) == 729L,
            !anyDuplicated(days, by = c("user_hash", "date")))
  knots <- as.numeric(quantile(sessions$day_max_c, c(.1,.5,.9), na.rm=TRUE))
  targets <- setNames(as.numeric(quantile(sessions$day_max_c,c(.5,.9,.95),na.rm=TRUE)),c("P50","P90","P95"))
  # Find the most frequent grid cell across eligible sessions. Break ties by the lowest cell ID.
  locations <- sessions[!is.na(cell_10km), .N, by=.(user_hash,cell_10km)]
  setorder(locations,user_hash,-N,cell_10km)
  usual <- locations[, .(cell_10km=cell_10km[1L], usual_share=N[1L]/sum(N)), by=user_hash]
  in_cell <- sessions[usual, on=.(user_hash,cell_10km), nomatch=0L]
  zones <- in_cell[!is.na(gps_iana_zone),.N,by=.(user_hash,gps_iana_zone)]
  setorder(zones,user_hash,-N,gps_iana_zone)
  zone <- zones[, .(zone=gps_iana_zone[1L]), by=user_hash]
  # Convert each UTC workout timestamp to the runner's usual time zone, accounting for clock changes.
  meta <- read_fitrec("cohort/workout_history.parquet",
    c("user_hash","sport","start_utc","has_hr_record","metadata_signature_redundant","eligible_recorded_running_history"))
  meta <- meta[user_hash %in% days$user_hash & is.finite(start_utc)]
  meta[zone,on="user_hash", zone:=i.zone]
  stopifnot(!anyNA(meta$zone))
  meta[, date:=as.IDate(format(as.POSIXct(as.numeric(start_utc), origin="1970-01-01", tz="UTC"),
                              "%Y-%m-%d",tz=zone[1L])),by=zone]
  any_runs <- unique(meta[eligible_recorded_running_history==TRUE,.(user_hash,date)])
  days[, run_any:=as.numeric(run)]
  days[any_runs,on=.(user_hash,date),run_any:=1]
  workouts <- meta[metadata_signature_redundant==FALSE & between(start_utc,946684800,1514764800)]
  sports <- list(any_sport=unique(workouts$sport), outdoor_running="run", treadmill_running="treadmill running",
    outdoor_cycling=c("bike","bike (transport)","mountain bike"), indoor_cycling="indoor cycling",
    walking_hiking=c("walk","hiking","fitness walking"),
    indoor_training=c("weight training","core stability training","circuit training","elliptical","aerobics",
                     "yoga","pilates","gymnastics","stair climing"), swimming="swimming",
    cycling_recreation=c("bike","mountain bike"), cycling_transport="bike (transport)")
  sports$indoor_all <- c("treadmill running","indoor cycling",sports$indoor_training)
  for (sport_name in names(sports)) {
    active <- unique(workouts[sport %in% sports[[sport_name]],.(user_hash,date)])
    days[, (sport_name):=0]
    days[active,on=.(user_hash,date),(sport_name):=1]
  }
  daily_totals <- sessions[,.(minutes=sum(duration_minutes,na.rm=TRUE),km=sum(gps_distance_km,na.rm=TRUE)),
                            by=.(user_hash,date=as.IDate(local_date))]
  days[daily_totals,on=.(user_hash,date),`:=`(minutes=i.minutes,km=i.km)]
  days[workouts[,.N,by=.(user_hash,date)],on=.(user_hash,date),any_workouts:=i.N]
  for (v in c("minutes","km","any_workouts")) set(days,which(is.na(days[[v]])),v,0)
  days[usual,on="user_hash",local_runner:=i.usual_share>=0.8]
  days[, iso_week:=format(date,"%G-%V")]
  days[, month:=format(date,"%Y-%m")]
  bins <- c("bin_below_0","bin_0_10","bin_20_25","bin_25_30","bin_30_plus")
  days[, (bins):=list(day_max_c<0, day_max_c>=0 & day_max_c<10, day_max_c>=20 & day_max_c<25,
                      day_max_c>=25 & day_max_c<30, day_max_c>=30)]
  weeks <- days[, c(list(days=.N,runs=sum(run),minutes=sum(minutes),km=sum(km),any_workouts=sum(any_workouts),
                 time=mean(time),south=first(south),local_runner=first(local_runner)),lapply(.SD,sum)),
                by=.(user_hash,iso_week),.SDcols=bins][days==7L]
  weeks[,month:=format(as.POSIXct(time,origin="1970-01-01",tz="UTC"),"%Y-%m")]

  # Define the controls for each of the four specifications. The monthly model
  # absorbs season and runner effects through a separate intercept for each runner and year-month.
  model_panel <- function(d,y,x,column=3L,variant="main",weekly=FALSE) {
    t <- as.POSIXct(d$time,origin="1970-01-01",tz="UTC")
    if (column>=2L && variant!="month") x <- cbind(x,season_terms(d$time,d$south))
    if (column>=2L && !weekly) {
      dow <- sapply(1:6,function(k) as.numeric(d$local_weekday==k))
      colnames(dow) <- paste0("dow",1:6)
      x <- cbind(x,dow)
    }
    if (column==4L) {
      yr <- as.integer(format(t,"%Y"))
      indicators <- sapply(2007:2016,function(k) as.numeric(yr==k))
      colnames(indicators) <- paste0("year",2007:2016)
      x <- cbind(x,indicators)
    }
    fe <- if (column<=2L) rep(1L,nrow(d)) else d$user_hash
    if (variant=="month") fe <- paste(d$user_hash,d$month)
    second <- NULL
    if (variant=="two_way") second <- if (weekly) d$iso_week else d$date
    fit_fe(d[[y]],x,fe,d$user_hash,cluster2=second)
  }
  # Days outside a runner's follow-up period contribute only to the preceding temperature history.
  source(file.path(repo,"record_level","scripts","panel_recording.R"),local=TRUE,encoding="UTF-8")
  weather <- read_fitrec("history/weighted/behaviour_daily_weighted.parquet",
                         c("user_hash","history_anchor_date","daily_max_c"))
  weather[, date:=as.IDate(history_anchor_date)]
  setkey(weather,user_hash,date)
  lag_day <- -3:7
  lag_names <- ifelse(lag_day<0,paste0("lead",abs(lag_day)),paste0("lag",lag_day))
  for (k in seq_along(lag_day)) {
    date_lookup <- days$date-lag_day[k]
    set(days,j=lag_names[k],value=weather[.(days$user_hash,date_lookup),daily_max_c])
  }
  stopifnot(max(abs(days$lag0-days$day_max_c),na.rm=TRUE)<1e-6)
  complete_days <- days[complete.cases(days[,..lag_names])]
  lag_x <- do.call(cbind,lapply(lag_names,function(nm) {
    x <- rcs_basis(complete_days[[nm]],knots)
    colnames(x) <- paste0(nm,"_D",1:2)
    x
  }))
  windows <- list(day_0=0,days_1_to_7=1:7,cumulative_days_0_to_7=0:7,leads_1_to_3_sum=-3:-1)
  lag_weights <- function(k,difference) {
    selected <- ifelse(k<0,paste0("lead",abs(k)),paste0("lag",k))
    names <- as.vector(rbind(paste0(selected,"_D1"),paste0(selected,"_D2")))
    setNames(rep(difference,length(k)),names)
  }
  labels <- c("P90","P95")
  delta <- lapply(targets[labels],function(t) drop(diff(rcs_basis(c(targets["P50"],t),knots))))
  variants <- c("c1","c2","c3","c4","month","two_way","local","any_run")
  daily_results <- list()
  profile <- list()
  for (variant in variants) {
    keep <- if (variant=="local") which(complete_days$local_runner) else seq_len(nrow(complete_days))
    column <- if (startsWith(variant,"c")) as.integer(substring(variant,2L)) else 3L
    y <- if (variant=="any_run") "run_any" else "run"
    fit <- model_panel(complete_days[keep],y,lag_x[keep,,drop=FALSE],column,variant)
    for (label in labels) for (q in names(windows)) {
      row <- contrast_panel(fit,lag_weights(windows[[q]],delta[[label]]))
      daily_results[[length(daily_results)+1L]] <- row[, `:=`(variant=variant,quantity=q,target_label=label)]
    }
    lead_terms <- names(lag_weights(-3:-1,c(1,1)))
    beta <- fit$beta[lead_terms]
    Fvalue <- drop(crossprod(beta,solve(fit$vcov[lead_terms,lead_terms],beta)))/length(beta)
    daily_results[[length(daily_results)+1L]] <- data.table(variant=variant, quantity="leads_joint_test",
       target_label="",F_statistic=Fvalue,numerator_df=length(beta),p_value=pf(Fvalue,length(beta),fit$df,lower.tail=FALSE))
    if (variant=="c3") for (label in labels) for (k in lag_day) {
      one <- contrast_panel(fit,lag_weights(k,delta[[label]]))
      setnames(one,names(one),paste0("day_",names(one)))
      if(k>=0) {
        cum <- contrast_panel(fit,lag_weights(0:k,delta[[label]]))
        setnames(cum,names(cum),paste0("cumulative_",names(cum)))
        one <- cbind(one,cum)
      }
      profile[[length(profile)+1L]] <- one[, `:=`(target_label=label,day_relative_to_hot_day=k)]
    }
  }
  daily <- rbindlist(daily_results,fill=TRUE)
  lag_profile <- rbindlist(profile,fill=TRUE)
  fwrite(daily,file.path(out,"04_lag_models.csv"))
  fwrite(lag_profile,file.path(out,"04_lag_profile.csv"))
  ref <- fread(input_file(activity_refs,"daily_model_adjustments.csv"))
  ref[,variant:=paste0("c",column)]
  fields <- c("estimate","se","ci_low","ci_high","p_value")
  compare_aggregate(daily[variant %in% paste0("c",1:4)],ref,c("variant","quantity","target_label"),
                    c(fields,"F_statistic","numerator_df"),"04_lag_ladder")
  ref <- fread(input_file(activity_refs,"daily_lag_profile.csv"))
  compare_aggregate(lag_profile,ref,c("target_label","day_relative_to_hot_day"),
                    c(paste0("day_",fields),paste0("cumulative_",fields)),"04_lag_profile")
  checks <- c(c3="column (3)", month="runner-by-year-month effects", two_way="two-way clustering, runner and date",
              local="column (3), runners at least 80% at the usual cell", any_run="column (3), any recorded run")
  df <- copy(daily[variant %in% names(checks)])
  df[,check:=unname(checks[variant])]
  compare_aggregate(df,fread(input_file(activity_refs,"daily_activity_sensitivity.csv")),c("check","quantity","target_label"),
                    c(fields,"F_statistic","numerator_df"),"04_lag_sensitivity")
  rm(lag_x,complete_days,weather,fit)
  gc()

  # Temperature bands for the weekly and daily models.
  bin_results <- list()
  for (y in c("runs","minutes","km","any_workouts","run")) {
    weekly <- y!="run"
    d <- if (weekly) weeks else days
    job <- if(y %in% c("run","any_workouts")) paste0("c",1:4) else c(paste0("c",1:4),"month","two_way","local")
    for (variant in job) {
      s <- if(variant=="local") d[local_runner==TRUE] else d
      column <- if(startsWith(variant,"c")) as.integer(substring(variant,2L)) else 3L
      x <- as.matrix(s[,..bins]); storage.mode(x) <- "double"
      fit <- model_panel(s,y,x,column,variant,weekly)
      for (band in bins) {
        row <- contrast_panel(fit,setNames(1,band))
        q <- paste0(if(weekly) "one_more_day_" else "day_in_",sub("^bin_","",band),"_instead_of_10_20")
        bin_results[[length(bin_results)+1L]] <- row[, `:=`(variant=variant,outcome=y,quantity=q)]
      }
    }
  }
  weekly_results <- rbindlist(bin_results)
  fwrite(weekly_results,file.path(out,"04_temperature_bands.csv"))
  ref <- fread(input_file(activity_refs,"weekly_model_adjustments.csv"))
  ref[,variant:=paste0("c",column)]
  compare_aggregate(weekly_results[variant %in% paste0("c",1:4)],ref,c("variant","outcome","quantity"),fields,"04_temperature_bands")
  checks["two_way"] <- "two-way clustering, runner and week"
  ws <- copy(weekly_results[variant %in% c("c3","month","two_way","local") & outcome %in% c("runs","minutes","km") &
                          grepl("25_30|30_plus",quantity)])
  ws[,check:=unname(checks[variant])]
  compare_aggregate(ws,fread(input_file(activity_refs,"weekly_activity_sensitivity.csv")),c("check","outcome","quantity"),fields,"04_weekly_sensitivity")

  # For each activity type, estimate the probability of at least one workout that day.
  dx <- rcs_basis(days$day_max_c,knots)
  colnames(dx) <- c("D1","D2")
  modalities <- list()
  for (y in names(sports)) {
    jobs <- if(y %in% c("cycling_recreation","cycling_transport")) "c3" else c(paste0("c",1:4),"month","two_way","local")
    for(variant in jobs) {
      keep <- if(variant=="local") which(days$local_runner) else seq_len(nrow(days))
      column <- if(startsWith(variant,"c")) as.integer(substring(variant,2L)) else 3L
      fit <- model_panel(days[keep],y,dx[keep,,drop=FALSE],column,variant)
      for(label in labels) {
        row <- contrast_panel(fit,setNames(delta[[label]],c("D1","D2")))
        modalities[[length(modalities)+1L]] <- row[, `:=`(variant=variant,outcome=y,target_label=label,
                                                         outcome_mean=mean(days[[y]][keep]))]
      }
    }
  }
  modalities <- rbindlist(modalities)
  fwrite(modalities,file.path(out,"04_sport.csv"))
  ref <- rbindlist(lapply(input_file(activity_refs,paste0("sport_model_adjustments_part",1:3,".csv")),fread))
  ref[,variant:=paste0("c",column)]
  compare_aggregate(modalities[variant %in% paste0("c",1:4) & !outcome %in% c("cycling_recreation","cycling_transport")],
                    ref,c("variant","outcome","target_label"),c(fields,"outcome_mean"),"04_sport")
  checks["two_way"] <- "two-way clustering, runner and date"
  sport_checks <- copy(modalities[variant %in% c("c3","month","two_way","local")])
  sport_checks[,check:=unname(checks[variant])]
  ref <- rbindlist(lapply(input_file(activity_refs,paste0("sport_sensitivity_part",1:4,".csv")),fread))
  compare_aggregate(sport_checks,ref,c("check","outcome","target_label"),c(fields,"outcome_mean"),"04_sport_sensitivity")
})
