# Reconstruct the physiological benchmark, timing and speed contributions.
# The two additional regressions use the released whole-cohort products.
local({
  source(file.path(repo,"scripts","model_tools.R"),local=TRUE,encoding="UTF-8")
  input <- file.path(repo,"data","processed","decomposition_model")
  primary <- file.path(repo,"data","processed","model_statistics")
  main <- fread(file.path(primary,"moments.csv"))
  temps <- fread(file.path(primary,"temperature_percentiles.csv"))
  design <- fread(file.path(input,"design.csv"))
  outcomes <- fread(file.path(input,"outcomes.csv"))
  covs <- fread(file.path(input,"coefficient_covariance.csv"))
  constants <- fread(file.path(input,"constants.csv"))
  a <- setNames(constants$value,constants$term)
  daily_knots <- temps[percentile %in% c(10,50,90)][order(percentile),temperature]
  session_knots <- a[paste0("session_knot",1:3)]
  speed_knots <- a[paste0("speed_knot",1:3)]
  median_day <- temps[percentile==50,temperature]
  slots <- list(gap=1:2,hour=3:4,speed=5:6,net=7:8,phys=9:12)
  b <- numeric(12)
  V <- matrix(0,12,12)
  for (label in c("gap","speed","net")) {
    outcome_name <- c(gap="heat_avoided_c",speed="gps_speed_median_m_s",net="hr_mean")[[label]]
    z <- main[outcome==outcome_name]
    A <- matrix(c(z$xx11,z$xx12,z$xx12,z$xx22),2,2)
    ix <- slots[[label]]
    b[ix] <- solve(A,c(z$xy1,z$xy2))
    meat <- matrix(c(z$meat11,z$meat12,z$meat12,z$meat22),2,2)
    V[ix,ix] <- solve(A)%*%meat%*%solve(A)
  }
  for (label in c("hour","phys")) {
    ix <- slots[[label]]
    A <- matrix(design[model==label][order(column,row),value],length(ix))
    b[ix] <- solve(A,outcomes[model==label][order(row),value])
  }
  V[1:4,1:4] <- matrix(covs[block=="gap_and_hour"][order(column,row),value],4)
  V[9:12,9:12] <- matrix(covs[block=="phys"][order(column,row),value],4)
  spline_difference <- function(after,before,knots,center=20,scale=10) {
    basis <- rcs_basis(c(after,before),knots,center=center,scale=scale)
    as.numeric(basis[1,]-basis[2,])
  }
  levels <- function(b,target) {
    daily <- spline_difference(target,median_day,daily_knots)
    gap0 <- unname(a["mean_gap"] + sum(a[c("reference_linear","reference_cubic")]*b[1:2]))
    t0 <- median_day+gap0
    observed_gap <- sum(daily*b[1:2])
    fixed_gap <- sum(daily*b[3:4])
    speed_change <- sum(daily*b[5:6])
    net <- sum(daily*b[7:8])
    benchmark <- sum(b[9:10]*spline_difference(target+gap0,t0,session_knots))
    fixed <- sum(b[9:10]*spline_difference(target+gap0+fixed_gap,t0,session_knots))
    timing <- sum(b[9:10]*spline_difference(target+gap0+observed_gap,t0,session_knots))
    speed <- sum(b[11:12]*spline_difference(a["speed_reference"]+speed_change,a["speed_reference"],
                                            speed_knots,center=0,scale=1))
    predicted <- timing+speed
    c(session_temperature_on_reference_day_c=t0,speed_change_m_s=speed_change,
      heat_avoided_change_c=observed_gap,heat_avoided_change_fixed_hour_c=fixed_gap,
      heat_avoided_change_behavioural_c=observed_gap-fixed_gap,
      behavioural_share_of_heat_avoided_change=(observed_gap-fixed_gap)/observed_gap,
      rise_without_behavioural_change=benchmark,rise_at_unchanged_clock_hour=fixed,
      diurnal_structure_part=benchmark-fixed,behavioural_timing_part=fixed-timing,
      rise_with_observed_timing=timing,speed_contribution=speed,
      predicted_rise_timing_and_speed=predicted,net_heart_rate_change=net,
      total_offset_fraction_fixed_hour=1-net/fixed,timing_offset_fraction_fixed_hour=1-timing/fixed,
      speed_offset_fraction_fixed_hour=-speed/fixed,unexplained_fraction_fixed_hour=(predicted-net)/fixed,
      structural_fraction_of_registered_benchmark=1-fixed/benchmark,total_offset_fraction_registered=1-net/benchmark,
      percentage_reduction_with_timing_and_speed=100*(1-predicted/fixed))
  }
  estimates <- list()
  for (target in c(90,95)) {
    temperature <- temps[percentile==target,temperature]
    point <- levels(b,temperature)
    gradient <- matrix(0,length(point),length(b))
    for (j in seq_along(b)) {
      step <- 1e-6*max(1,abs(b[j]))
      above <- below <- b
      above[j] <- above[j]+step
      below[j] <- below[j]-step
      gradient[,j] <- (levels(above,temperature)-levels(below,temperature))/(2*step)
    }
    se <- sqrt(pmax(0,rowSums((gradient%*%V)*gradient)))
    estimates[[as.character(target)]] <- data.table(target=paste0("P",target),quantity=names(point),
      estimate=as.numeric(point),se,lower=as.numeric(point)-qnorm(.975)*se,upper=as.numeric(point)+qnorm(.975)*se)
  }
  dest <- file.path(public_out,"models")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  fwrite(rbindlist(estimates),file.path(dest,"heart_rate_decomposition.csv"))
})
