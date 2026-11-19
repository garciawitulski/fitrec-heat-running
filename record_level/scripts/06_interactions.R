# Test how the association with heat varies by weekend status and prior ambient exposure.
# Each session has equal weight. Contrasts compare P90 and P95 with the median of
# the exposure within runners, with standard errors clustered by runner.
local({
  reference_dir <- "references/main"
  history_dir <- "references/heart_rate_and_history"
  config <- jsonlite::fromJSON(input_file("settings/session_model.json"))
  fields <- c("estimate","se","ci_low","ci_high","p_value")
  d <- read_fitrec("panels/sessions.parquet",
       c("day_max_c","canonical_start_utc","geo_mean_y","weekend","user_hash",
         "hr_mean","heat_avoided_c","gps_speed_median_m_s","duration_minutes"))[is.finite(day_max_c)]
  d[,log_duration:=log(duration_minutes)]
  pct <- c(.5,.9,.95)
  labels <- c("P50","P90","P95")
  target <- setNames(as.numeric(quantile(d$day_max_c,pct)),labels)
  knots <- unlist(config$knots$rcs3)
  x <- rcs_basis(d$day_max_c,knots)
  colnames(x) <- c("D1","D2")
  weekends <- as.numeric(d$weekend)
  X <- cbind(x, weekend=weekends, W1=weekends*x[,1], W2=weekends*x[,2],
             season_terms(d$canonical_start_utc,d$geo_mean_y<0))
  outcomes <- c(w_b1_heat_avoided="heat_avoided_c",w_b2_speed="gps_speed_median_m_s",
                w_b3_log_duration="log_duration",w_h1_hr_mean="hr_mean")
  # Read the temperatures at which the original figure displays the curves.
  grid <- fread(input_file(reference_dir,"weekend_curves.csv"),select=c("model_id","exposure_c"))
  all_weekend <- list()
  all_curves <- list()
  joint_weekend <- list()
  for (mid in names(outcomes)) {
    y <- d[[outcomes[[mid]]]]
    fit <- fit_fe(y,X,d$user_hash,d$user_hash)
    fit_terms <- c("W1","W2")
    b <- fit$beta[fit_terms]
    f <- drop(crossprod(b,solve(fit$vcov[fit_terms,fit_terms],b)))/2
    joint_weekend[[mid]] <- data.table(model_id=mid,F_statistic=f,p_value=pf(f,2,fit$df,lower.tail=FALSE))
    for(temp in unique(c(target[labels[-1L]],grid[model_id==mid,exposure_c]))) {
      delta <- as.numeric(diff(rcs_basis(c(target["P50"],temp),knots)))
      base <- setNames(delta,c("D1","D2"))
      interaction <- setNames(delta,c("W1","W2"))
      weights <- list(weekday=base,weekend=c(base,interaction),weekend_minus_weekday=interaction)
      for(stratum_name in names(weights)) {
        row <- contrast_panel(fit,weights[[stratum_name]])
        row[, `:=`(model_id=mid,stratum=stratum_name,exposure_c=temp)]
        all_curves[[length(all_curves)+1L]] <- row
        if(temp %in% target[labels[-1L]]) {
          labelled <- copy(row)
          labelled[, `:=`(target_label=names(target)[match(temp,target)],contrast=paste0("R_",stratum_name))]
          all_weekend[[length(all_weekend)+1L]] <- labelled
        }
      }
    }
  }
  weekend <- rbindlist(all_weekend)
  weekend[contrast=="R_weekend_minus_weekday",p_holm:=p.adjust(p_value,method="holm")]
  weekend_curves <- rbindlist(all_curves)
  weekend_curves <- merge(weekend_curves,unique(grid),by=c("model_id","exposure_c"))
  fwrite(weekend,file.path(out,"06_weekend.csv"))
  fwrite(weekend_curves,file.path(out,"06_weekend_curves.csv"))
  fwrite(rbindlist(joint_weekend),file.path(out,"06_weekend_joint.csv"))
  compare_aggregate(weekend,fread(input_file(reference_dir,"weekend_contrasts.csv")),
                    c("model_id","target_label","contrast"),fields,"06_weekend")
  compare_aggregate(weekend_curves,fread(input_file(reference_dir,"weekend_curves.csv")),
                    c("model_id","stratum","exposure_c"),fields,"06_weekend_curves")
  test_ref <- as.data.table(jsonlite::fromJSON(input_file(reference_dir,"joint_tests.json")))
  compare_aggregate(rbindlist(joint_weekend),
                    test_ref[model_id %in% names(outcomes) & test=="weekend_by_exposure_interaction"],
                    "model_id",c("F_statistic","p_value"),"06_weekend_joint")
  # Table S11 reports these four adjusted p-values to three decimal places.
  printed <- data.table(model_id=names(outcomes),p_holm=c(.096,1,1,.046))
  compare_aggregate(weekend[target_label=="P95" & contrast=="R_weekend_minus_weekday"],printed,
                    "model_id","p_holm","06_weekend_holm_printed",tolerance=.00050001)

  # Each design matrix has one row per session and includes the predictors for its exposure window.
  # Refit the models using the outcome y, temperature terms D1 and D2, and prior heat exposure.
  # AB denotes behaviour versus daily maximum temperature. AH denotes heart rate at a fixed
  # speed versus session temperature. EW denotes an exponentially weighted mean with a seven-day half-life.
  dirs <- c(AB="history/daily_maximum",AH="history/session_temperature",EW="history/weighted")
  model_ids <- c(as.vector(outer(c("AB14","AB07"),c("B1","B2","B3","H1","START"),paste,sep="_")),
                 "AH14","AH07",paste0("EW14_",c("B1","B2","B3","H1","START","H2")))
  measures <- list()
  for (mid in c("EW14_B1","EW14_H2")) {
    data <- read_fitrec(file.path(dirs["EW"],paste0("fit_design_",mid,".parquet")), "D1")
    measures[[mid]] <- setNames(as.numeric(quantile(20+10*data$D1,pct)),labels)
  }
  hist_grid <- fread(input_file(history_dir,"history_curves.csv"),select=c("model_id","target_c"))
  hist_contrast <- hist_level <- hist_test <- hist_curve <- list()
  for (mid in model_ids) {
    family <- substr(mid,1,2)
    h2 <- mid %in% c("AH07","AH14","EW14_H2")
    anchor <- measures[[if(h2) "EW14_H2" else "EW14_B1"]]
    knots <- unlist(if(h2) config$knots$session_temp_rcs3 else config$knots$rcs3)
    a <- read_fitrec(file.path(dirs[family],paste0("fit_design_",mid,".parquet")))
    outcome <- if("hr_lag60" %in% names(a)) "hr_lag60" else "y"
    cols <- setdiff(names(a),c("user_hash",outcome,grep("^year_",names(a),value=TRUE)))
    fit <- fit_fe(a[[outcome]],as.matrix(a[,..cols]),a$user_hash,a$user_hash)
    a[, history_c:=20+5*H]
    a[, history_within:=history_c-mean(history_c),by=user_hash]
    increment <- round(sqrt(mean(a$history_within^2))*2)/2
    hist_level[[mid]] <- contrast_panel(fit,c(H=1))[,model_id:=mid]
    hterms <- c("H:D1","H:D2")
    b <- fit$beta[hterms]
    f <- drop(crossprod(b,solve(fit$vcov[hterms,hterms],b)))/2
    hist_test[[mid]] <- data.table(model_id=mid,F_statistic=f,p_value=pf(f,2,fit$df,lower.tail=FALSE))
    grid_values <- unique(hist_grid[model_id==mid,target_c])
    for (temp in unique(c(anchor[labels[-1L]],grid_values))) {
      delta <- as.numeric(diff(rcs_basis(c(anchor["P50"],temp),knots)))
      # The temperature grid in the CSV may round the median.
      if (abs(temp-anchor["P50"]) < 1e-10) delta[] <- 0
      row <- contrast_panel(fit,setNames(delta,hterms))
      row[, `:=`(model_id=mid,target_c=temp)]
      if(temp %in% grid_values) hist_curve[[length(hist_curve)+1L]] <- row
      if(temp %in% anchor[labels[-1L]]) {
        row <- copy(row)
        row[, `:=`(target_label=names(anchor)[match(temp,anchor)],supported_increment_c=increment)]
        for (stat in c("estimate","ci_low","ci_high")) {
          val <- row[[stat]]
          display <- if (endsWith(mid,"_B3")) 100*expm1(val) else val
          scaled <- val*increment/5
          if (endsWith(mid,"_B3")) scaled <- 100*expm1(scaled)
          row[, (paste0("display_",stat)):=display]
          row[, (paste0("per_increment_",stat)):=scaled]
        }
        hist_contrast[[length(hist_contrast)+1L]] <- row
      }
    }
  }
  hist_contrast <- rbindlist(hist_contrast)
  hist_level <- rbindlist(hist_level)
  hist_test <- rbindlist(hist_test)
  hist_curve <- rbindlist(hist_curve)
  for(nm in c("hist_contrast","hist_level","hist_test","hist_curve"))
    fwrite(get(nm),file.path(out,paste0("06_",nm,".csv")))
  compare_aggregate(hist_level,fread(input_file(history_dir,"history_levels.csv")),"model_id",fields,"06_history_levels")
  compare_aggregate(hist_test,fread(input_file(history_dir,"history_joint_tests.csv")),"model_id",c("F_statistic","p_value"),"06_history_tests")
  compare_aggregate(hist_curve,fread(input_file(history_dir,"history_curves.csv")),c("model_id","target_c"),fields,"06_history_curves")
  measures_out <- c(fields,"target_c","supported_increment_c",paste0("display_",c("estimate","ci_low","ci_high")),
                    paste0("per_increment_",c("estimate","ci_low","ci_high")))
  compare_aggregate(hist_contrast,fread(input_file(history_dir,"history_contrasts.csv")),c("model_id","target_label"),
                    measures_out,"06_history_contrasts")
})
