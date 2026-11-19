# Technical checks for the competing heart-rate measurements (Tables S13, S14 and S17).
# Each fit below uses the precomputed session measurements and one observation per
# session. Exclusions for a missing outcome are specific to that outcome.
hr_contrasts <- function(fit, id, session = FALSE) {
    anchors <- if (session) hr_stargets else hr_targets
    knots <- if (session) hr_ksession else hr_kday
    rbindlist(lapply(c("P90","P95"),function(label) {
        v <- setNames(hr_spline_change(anchors[[label]], anchors[["P50"]], knots),c("T1","T2"))
        lin_contrast(fit,v)[, ':='(model_id=id,target_label=label)]
    }))
}

# Alternative physiological specifications and common-session benchmarks.
x <- hr_design(hr_phys,"temp_c_10km",speed=TRUE,training=TRUE,hour="session_start_minutes_local")
hour_fit <- fit_fe(hr_phys$hr_lag60,x,hr_phys$user_hash,hr_phys$user_hash)
x <- hr_design(hr_day,"temp_c_10km",speed=TRUE,training=TRUE,weekday=FALSE)
overlap_temp <- fit_fe(hr_day$hr_lag60,x,hr_day$user_hash,hr_day$user_hash)
overlap_net <- fit_fe(hr_day$session_hr_mean,hr_design(hr_day),hr_day$user_hash,hr_day$user_hash)
benchmark <- rbindlist(list(
    hr_contrasts(hour_fit,"h2_hour_weekday",TRUE),hr_contrasts(overlap_temp,"ov_h2_session_temp",TRUE),
    hr_contrasts(hr_fixed,"ov_segment_hr_daymax_fixed_speed"),hr_contrasts(hr_free,"ov_segment_hr_daymax_speed_free"),
    hr_contrasts(hr_whole,"ov_session_hr_daymax"),hr_contrasts(overlap_net,"ov_h1_main_spec")))
fwrite(benchmark,file.path(out,"heart_benchmark.csv"))
models_hour <- hr_models
models_hour$phys <- hour_fit
month_id <- paste(hr_main$user_hash,hr_main$year_month_local,sep="|")
month_x <- hr_design(hr_main,season=FALSE)
month_hour_x <- hr_design(hr_main,season=FALSE,hour="start_minutes_local")
models_month <- list(
    gap=fit_fe(hr_main$heat_avoided_c,month_x,month_id,hr_main$user_hash),
    hour=fit_fe(hr_main$heat_avoided_c,month_hour_x,month_id,hr_main$user_hash),
    speed=fit_fe(hr_main$gps_speed_median_m_s,month_x,month_id,hr_main$user_hash),
    net=fit_fe(hr_main$hr_mean,month_x,month_id,hr_main$user_hash),phys=hr_phfit)
for (method in c("hour","month")) {
    models <- if (method=="hour") models_hour else models_month
    results <- rbindlist(lapply(c("P90","P95"),function(label)
        hr_decomposition(models,hr_targets[[label]])[,target_label:=label]))
    fwrite(results,file.path(out,paste0("heart_",method,"_decomposition.csv")))
}

# Screened arithmetic and time-weighted means share the same raw-record cache.
# Confirm the unchanged arithmetic mean for every eligible session before linking.
tw <- rbindlist(lapply(1:4,function(k) read_fitrec(sprintf(
    "panels/time_weighted_hr_part%d.parquet",k),
    c("session_hash","hr_mean_reproduced","hr_mean_time_weighted","hr_mean_screened"))))
stopifnot(!anyDuplicated(tw$session_hash),setequal(tw$session_hash,hr_all$session_hash))
original <- hr_all$hr_mean[match(tw$session_hash,hr_all$session_hash)]
stopifnot(identical(is.na(original),is.na(tw$hr_mean_reproduced)),
          max(abs(original-tw$hr_mean_reproduced),na.rm=TRUE)<1e-9)
matching <- match(hr_main$session_hash,tw$session_hash)
weighted <- list()
for (outcome in c("hr_mean_time_weighted","hr_mean_screened")) {
    fit <- fit_fe(tw[[outcome]][matching],hr_xmain,hr_main$user_hash,hr_main$user_hash)
    weighted[[outcome]] <- hr_contrasts(fit,paste0("h1_",outcome))
}
fwrite(rbindlist(weighted),file.path(out,"heart_time_weighted.csv"))
common_y <- tw$hr_mean_time_weighted[match(hr_day$session_hash,tw$session_hash)]
ok <- is.finite(common_y)
runner <- hr_day$user_hash[ok]
weighted_models <- list(fixed=fit_fe(hr_day$hr_lag60[ok],hr_xfixed[ok,],runner,runner),
                         free=fit_fe(hr_day$hr_lag60[ok],hr_xfree[ok,],runner,runner),
                         whole=fit_fe(common_y[ok],hr_xfree[ok,],runner,runner))
s <- hr_stack(weighted_models)
paired_weighted <- list()
for (label in c("P90","P95")) {
    v <- hr_spline_change(hr_targets[[label]],hr_targets[["P50"]],hr_kday)
    effects <- function(b) {
        fixed <- sum(v*b[c("fixed.T1","fixed.T2")])
        free <- sum(v*b[c("free.T1","free.T2")])
        whole <- sum(v*b[c("whole.T1","whole.T2")])
        c(steady_fixed_speed=fixed,steady_speed_free=free,whole_session_time_weighted=whole,
          steady_speed_free_minus_time_weighted=free-whole,steady_fixed_speed_minus_time_weighted=fixed-whole)
    }
    paired_weighted[[label]] <- hr_delta(effects,s$beta,s$vcov,weighted_models$free$df)[,target_label:=label]
}
fwrite(rbindlist(paired_weighted),file.path(out,"heart_time_weighted_pairs.csv"))

# Resolution of the recorded trace: number of points and typical sampling interval.
res <- rbindlist(lapply(1:3,function(k) read_fitrec(sprintf(
    "covariates/resolution_part%d.parquet",k),c("session_hash","n_points","median_dt_s"))))
stopifnot(!anyDuplicated(res$session_hash), all(c(hr_main$session_hash,hr_phys$session_hash) %in% res$session_hash))
res[, capped := as.numeric(n_points==500L)]
match_main <- match(hr_main$session_hash,res$session_hash)
res_dt <- fit_fe(log(res$median_dt_s[match_main]),hr_xmain,hr_main$user_hash,hr_main$user_hash)
res_cap <- fit_fe(res$capped[match_main],hr_xmain,hr_main$user_hash,hr_main$user_hash)
fwrite(rbind(hr_contrasts(res_dt,"a4b_log_median_dt"),hr_contrasts(res_cap,"a4b_capped")),
       file.path(out,"heart_recording_resolution.csv"))
cutpoints <- quantile(res$median_dt_s[match(hr_phys$session_hash,res$session_hash)],c(1/3,2/3),na.rm=TRUE)
group_results <- group_tests <- list()
for (paired in c(FALSE,TRUE)) {
    d <- if (paired) hr_day else hr_phys
    original_x <- if (paired) hr_xfree else hr_xphys
    y <- if (paired) d$hr_lag60-d$session_hr_mean else d$hr_lag60
    anchors <- if (paired) hr_targets else hr_stargets
    knots <- if (paired) hr_kday else hr_ksession
    prefix <- if (paired) "a4c_paired_speed_free" else "h2_pooled_session_temp"
    lookup <- match(d$session_hash,res$session_hash)
    interval <- res$median_dt_s[lookup]
    for (group in c("dt_tercile","record")) {
        levels <- if (group=="dt_tercile") c("t1","t2","t3") else c("uncapped","capped")
        g <- if (group=="dt_tercile") {
            fifelse(interval<=cutpoints[1],"t1",fifelse(interval<=cutpoints[2],"t2","t3"))
        } else ifelse(res$capped[lookup]==1,"capped","uncapped")
        X <- original_x
        restrictions <- character()
        for (level in levels[-1L]) {
            indicator <- as.numeric(g==level)
            extra <- cbind(indicator,original_x[,c("T1","T2")]*indicator)
            colnames(extra) <- paste(level,c("level","T1","T2"),sep="_")
            X <- cbind(X,extra)
            restrictions <- c(restrictions,colnames(extra)[2:3])
        }
        id <- paste0(prefix,"_by_",group)
        fit <- fit_fe(y,X,d$user_hash,d$user_hash)
        for (label in c("P90","P95")) for (level in levels) {
            v <- setNames(hr_spline_change(anchors[[label]],anchors[["P50"]],knots),c("T1","T2"))
            if (level!=levels[1]) v <- c(v,setNames(v,paste(level,names(v),sep="_")))
            group_results[[length(group_results)+1L]] <- lin_contrast(fit,v)[,
                ':='(model_id=id,target_label=label,group_level=level)]
        }
        b <- fit$beta[restrictions]
        df1 <- length(b)
        test <- drop(crossprod(b,solve(fit$vcov[restrictions,restrictions],b)))/df1
        group_tests[[id]] <- data.table(model_id=id,F_statistic=test,p=pf(test,df1,fit$df,lower.tail=FALSE))
    }
}
fwrite(rbindlist(group_results),file.path(out,"heart_resolution_groups.csv"))
fwrite(rbindlist(group_tests),file.path(out,"heart_resolution_tests.csv"))
