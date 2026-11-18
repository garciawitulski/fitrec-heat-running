# Compare mean heart rate during steady-speed segments with the whole-session mean.
# Each session has equal weight. The models contrast P90 and P95 with P50 within
# runners, with runner-clustered errors and seasonal terms aligned by hemisphere.
# The paired models use the same sessions and controls. The fixed-speed model
# adds a spline for running speed.

hr_phys <- read_fitrec("panels/paired_heart_rate.parquet")
hr_all <- read_fitrec("panels/sessions.parquet",
    c("session_hash", "user_hash", "day_max_c", "heat_avoided_c", "gps_speed_median_m_s",
      "hr_mean", "canonical_start_utc", "geo_mean_y", "local_weekday",
      "start_minutes_local", "year_month_local", "cell_10km", "duration_minutes"))
hr_main <- hr_all[is.finite(day_max_c)]
hr_day <- hr_phys[is.finite(day_max_c)]
stopifnot(nrow(hr_phys) == 19350L, nrow(hr_main) == 106600L, nrow(hr_day) == 18626L,
          !anyDuplicated(hr_all$session_hash), !anyDuplicated(hr_phys$session_hash),
          all(is.finite(hr_day$session_hr_mean)), sum(is.finite(hr_main$hr_mean)) == 106341L)
message("The heart-rate sample includes 19,350 steady-speed sessions, of which 18,626 have daily maximum temperature.")

hr_config <- jsonlite::fromJSON(input_file("settings/session_model.json"))
hr_kday <- hr_config$knots$rcs3
hr_ksession <- hr_config$knots$session_temp_rcs3
hr_kspeed <- hr_config$speed_knots
hr_speed0 <- hr_config$typical_segment_speed_m_s
hr_targets <- quantile(hr_main$day_max_c, c(.5, .9, .95), names = FALSE)
hr_stargets <- quantile(hr_phys$temp_c_10km, c(.5, .9, .95), names = FALSE)
names(hr_targets) <- names(hr_stargets) <- c("P50", "P90", "P95")

# Consistent column names identify each coefficient in the decomposition.
hr_design <- function(d, exposure = "day_max_c", speed = FALSE, training = FALSE,
                      weekday = TRUE, season = TRUE, hour = NULL) {
    knots <- if (exposure == "day_max_c") hr_kday else hr_ksession
    X <- rcs_basis(d[[exposure]], knots)
    colnames(X) <- c("T1", "T2")
    if (season) X <- cbind(X, season_terms(d$canonical_start_utc, d$geo_mean_y < 0))
    if (weekday) {
        wd <- model.matrix(~ factor(local_weekday, levels = 0:6), data = d)[, -1L, drop = FALSE]
        X <- cbind(X, wd)
    }
    if (speed) {
        sp <- rcs_basis(d$speed_lag60, hr_kspeed, center = 0, scale = 1)
        colnames(sp) <- c("V1", "V2")
        X <- cbind(X, sp)
    }
    if (training) X <- cbind(X, prior_minutes = d$prior_recorded_minutes, fitness = d$fitness_speed_kmh_28d)
    if (!is.null(hour)) {
        # Use the most common start hour as the reference for the hour indicators.
        h <- factor(floor(d[[hour]] / 60), levels = 0:23)
        h <- relevel(h, ref = names(which.max(table(h))))
        X <- cbind(X, model.matrix(~ h)[, -1L, drop = FALSE])
    }
    X
}

# Fit the three paired models with daily maximum temperature as the exposure.
hr_xfree <- hr_design(hr_day, training = TRUE, weekday = FALSE)
hr_xfixed <- hr_design(hr_day, speed = TRUE, training = TRUE, weekday = FALSE)
hr_free <- fit_fe(hr_day$hr_lag60, hr_xfree, hr_day$user_hash, hr_day$user_hash)
hr_fixed <- fit_fe(hr_day$hr_lag60, hr_xfixed, hr_day$user_hash, hr_day$user_hash)
hr_whole <- fit_fe(hr_day$session_hr_mean, hr_xfree, hr_day$user_hash, hr_day$user_hash)

# The physiological model uses temperature during the session.
hr_xphys <- hr_design(hr_phys, "temp_c_10km", speed = TRUE, training = TRUE, weekday = FALSE)
hr_phfit <- fit_fe(hr_phys$hr_lag60, hr_xphys, hr_phys$user_hash, hr_phys$user_hash)

# Each model uses main-cohort sessions with an observed value for its outcome.
# Only the fixed-hour model includes start-hour indicators.
hr_xmain <- hr_design(hr_main)
hr_xhour <- hr_design(hr_main, hour = "start_minutes_local")
hr_gap <- fit_fe(hr_main$heat_avoided_c, hr_xmain, hr_main$user_hash, hr_main$user_hash)
hr_hour <- fit_fe(hr_main$heat_avoided_c, hr_xhour, hr_main$user_hash, hr_main$user_hash)
hr_speed <- fit_fe(hr_main$gps_speed_median_m_s, hr_xmain, hr_main$user_hash, hr_main$user_hash)
hr_valid <- is.finite(hr_main$hr_mean)
hr_net <- fit_fe(hr_main$hr_mean[hr_valid], hr_xmain[hr_valid, , drop = FALSE],
                 hr_main$user_hash[hr_valid], hr_main$user_hash[hr_valid])

# Each row collects one runner's contributions to the coefficient estimates.
# A runner absent from one sample contributes zero to the cross-model covariance.
hr_stack <- function(models) {
    ids <- sort(unique(unlist(lapply(models, function(m) rownames(m$score)), use.names = FALSE)))
    pieces <- list()
    coefficients <- numeric()
    for (label in names(models)) {
        model <- models[[label]]
        vars <- intersect(c("T1", "T2", "V1", "V2"), names(model$beta))
        keys <- paste0(label, ".", vars)
        influence <- matrix(0, nrow = length(ids), ncol = length(vars), dimnames = list(ids, keys))
        influence[rownames(model$score), ] <- model$score[, vars, drop = FALSE]
        pieces[[label]] <- influence
        coefficients <- c(coefficients, setNames(model$beta[vars], keys))
    }
    list(beta = coefficients, vcov = crossprod(do.call(cbind, pieces)))
}

# The delta method uses the derivatives of the fitted levels and their differences.
hr_delta <- function(f, beta, covariance, df = Inf) {
    point <- f(beta)
    grad <- matrix(0, length(point), length(beta))
    for (j in seq_along(beta)) {
        h <- 1e-6 * max(1, abs(beta[j]))
        plus <- minus <- beta
        plus[j] <- plus[j] + h
        minus[j] <- minus[j] - h
        grad[, j] <- (f(plus) - f(minus)) / (2 * h)
    }
    se <- sqrt(pmax(0, rowSums((grad %*% covariance) * grad)))
    critical <- qt(.975, df)
    data.table(quantity = names(point), estimate = as.numeric(point), se = se,
               lower = as.numeric(point) - critical * se, upper = as.numeric(point) + critical * se,
               p = 2 * pt(-abs(as.numeric(point) / se), df))
}

hr_spline_change <- function(end, start, knots, speed = FALSE) {
    values <- c(end, start)
    basis <- if (speed) rcs_basis(values, knots, center = 0, scale = 1) else rcs_basis(values, knots)
    as.numeric(basis[1L, ] - basis[2L, ])
}

hr_pair_models <- list(fixed = hr_fixed, free = hr_free, whole = hr_whole)
hr_pair <- hr_stack(hr_pair_models)
hr_paired_rows <- list()
for (lab in c("P90", "P95")) {
    v <- hr_spline_change(hr_targets[[lab]], hr_targets[["P50"]], hr_kday)
    paired_change <- function(b) {
        fixed <- sum(v * b[c("fixed.T1", "fixed.T2")])
        free <- sum(v * b[c("free.T1", "free.T2")])
        whole <- sum(v * b[c("whole.T1", "whole.T2")])
        c(steady_fixed_speed = fixed, steady_speed_free = free, whole_session_mean = whole,
          steady_speed_free_minus_whole_session = free - whole,
          steady_fixed_speed_minus_whole_session = fixed - whole, fixed_speed_minus_speed_free = fixed - free)
    }
    hr_paired_rows[[lab]] <- hr_delta(paired_change, hr_pair$beta, hr_pair$vcov, hr_free$df)[, target_label := lab]
}
hr_paired_rows <- rbindlist(hr_paired_rows)
fwrite(hr_paired_rows, file.path(out, "heart_paired.csv"))

# Tables S8–S9 and Figure 5. Decompose the heart-rate response at a fixed start hour
# and as timing and speed vary. Estimate the mean thermal gap on a day at P50.
hr_reference_basis <- as.numeric(rcs_basis(hr_targets[["P50"]], hr_kday)) - colMeans(hr_xmain[, c("T1","T2")])
hr_gap_mean <- mean(hr_main$heat_avoided_c)
hr_models <- list(gap = hr_gap, hour = hr_hour, speed = hr_speed, net = hr_net, phys = hr_phfit)
hr_parts <- function(b, target, mean_gap, reference_basis) {
    v <- hr_spline_change(target, hr_targets[["P50"]], hr_kday)
        b_gap <- b[c("gap.T1","gap.T2")]
        b_heat <- b[c("phys.T1","phys.T2")]
        b_pace <- b[c("phys.V1","phys.V2")]
        gap0 <- mean_gap + sum(reference_basis * b_gap)
        t0 <- hr_targets[["P50"]] + gap0
        gap_obs <- sum(v * b_gap)
        gap_fixed <- sum(v * b[c("hour.T1","hour.T2")])
        pace <- sum(v * b[c("speed.T1","speed.T2")])
        net <- sum(v * b[c("net.T1","net.T2")])
        original <- sum(b_heat * hr_spline_change(target + gap0, t0, hr_ksession))
        fixed <- sum(b_heat * hr_spline_change(target + gap0 + gap_fixed, t0, hr_ksession))
        timing <- sum(b_heat * hr_spline_change(target + gap0 + gap_obs, t0, hr_ksession))
        slower <- sum(b_pace * hr_spline_change(hr_speed0 + pace, hr_speed0, hr_kspeed, TRUE))
        predicted <- timing + slower
        c(session_temperature_on_reference_day_c = t0, speed_change_m_s = pace,
          heat_avoided_change_c = gap_obs, heat_avoided_change_fixed_hour_c = gap_fixed,
          heat_avoided_change_behavioural_c = gap_obs - gap_fixed,
          behavioural_share_of_heat_avoided_change = (gap_obs - gap_fixed) / gap_obs,
          rise_without_behavioural_change = original, rise_at_unchanged_clock_hour = fixed,
          diurnal_structure_part = original - fixed, behavioural_timing_part = fixed - timing,
          rise_with_observed_timing = timing, speed_contribution = slower,
          predicted_rise_timing_and_speed = predicted, net_heart_rate_change = net,
          total_offset_fraction_fixed_hour = 1 - net / fixed, timing_offset_fraction_fixed_hour = 1 - timing / fixed,
          speed_offset_fraction_fixed_hour = -slower / fixed, unexplained_fraction_fixed_hour = (predicted - net) / fixed,
          structural_fraction_of_registered_benchmark = 1 - fixed / original, total_offset_fraction_registered = 1 - net / original)
}
hr_decomposition <- function(models, target, mean_gap = hr_gap_mean, reference_basis = hr_reference_basis) {
    s <- hr_stack(models)
    # The published intervals include covariance between the thermal-gap and fixed-hour
    # fits. The remaining model blocks are treated as independent in this decomposition.
    block <- sub("\\..*$", "", names(s$beta))
    block[block == "hour"] <- "gap"
    covariance <- s$vcov
    covariance[outer(block, block, "!=")] <- 0
    parts <- function(b) hr_parts(b, target, mean_gap, reference_basis)
    hr_delta(parts, s$beta, covariance)
}
hr_dec <- rbindlist(lapply(c("P90", "P95"), function(lab)
    hr_decomposition(hr_models, hr_targets[[lab]])[, target_label := lab]))
fwrite(hr_dec, file.path(out, "heart_decomposition.csv"))

# Table S10. Relate predicted heart rate at the observed start times to the observed
# whole-session mean. Use joint covariance across all models, including across samples.
hr_chain <- hr_stack(c(list(gap = hr_gap, phys = hr_phfit), hr_pair_models, list(net = hr_net)))
hr_chain_rows <- list()
for (lab in c("P90", "P95")) {
    target <- hr_targets[[lab]]
    v <- hr_spline_change(target, hr_targets[["P50"]], hr_kday)
    chain_change <- function(b) {
        bg <- b[c("gap.T1", "gap.T2")]
        gap0 <- hr_gap_mean + sum(hr_reference_basis * bg)
        temperature <- hr_spline_change(target + gap0 + sum(v * bg), hr_targets[["P50"]] + gap0, hr_ksession)
        first <- sum(temperature * b[c("phys.T1", "phys.T2")])
        fixed <- sum(v * b[c("fixed.T1","fixed.T2")])
        free <- sum(v * b[c("free.T1","free.T2")])
        whole <- sum(v * b[c("whole.T1","whole.T2")])
        net <- sum(v * b[c("net.T1","net.T2")])
        c(rise_with_observed_timing = first, steady_fixed_speed = fixed, steady_speed_free = free,
          whole_session_common_sessions = whole, whole_session_main = net, link_transfer = first - fixed,
          link_speed = fixed - free, link_measure = free - whole, link_sample_and_covariates = whole - net,
          total_rise_with_timing_minus_main = first - net)
    }
    hr_chain_rows[[lab]] <- hr_delta(chain_change, hr_chain$beta, hr_chain$vcov)[, target_label := lab]
}
hr_chain_rows <- rbindlist(hr_chain_rows)
fwrite(hr_chain_rows, file.path(out, "heart_chain.csv"))
