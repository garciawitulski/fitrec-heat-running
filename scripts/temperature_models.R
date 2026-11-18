# Table 2. Solve each two-term temperature model from pooled cross-products
# after removing the fixed effects and nuisance controls. No participant data
# are needed. The root analysis script supplies repo and public_out.
local({
  library(data.table)
  source(file.path(repo,"scripts","model_tools.R"),local=TRUE,encoding="UTF-8")
  input <- file.path(repo,"data","processed","model_statistics")
  m <- fread(file.path(input,"moments.csv"))
  a <- fread(file.path(input,"temperature_percentiles.csv"))
  knots <- a[percentile %in% c(10,50,90),temperature]
  targets <- a[percentile %in% c(50,90,95),temperature]
  grid <- seq(ceiling(2*a[percentile==1,temperature])/2,
              floor(2*a[percentile==99,temperature])/2,by=.5)
  grid <- sort(unique(c(grid,targets)))
  basis <- rcs_basis(grid,knots)
  contrast <- sweep(basis,2L,as.numeric(rcs_basis(targets[1],knots)),"-")

  # The paper identifies 820 runners in the main cohort and 729 in the daily
  # panel. For heart rate and steady-speed sessions it reports lower bounds.
  # Use every possible t interval between that bound and the parent cohort.
  curve_rows <- list()
  for (i in seq_len(nrow(m))) {
    z <- m[i]
    outcome <- z$outcome
    gram <- matrix(c(z$xx11,z$xx12,z$xx12,z$xx22),2,2)
    rhs <- c(z$xy1,z$xy2)
    meat <- matrix(c(z$meat11,z$meat12,z$meat12,z$meat22),2,2)
    beta <- solve(gram,rhs)
    bread <- solve(gram)
    covariance <- bread %*% meat %*% bread
    effect <- as.numeric(contrast %*% beta)
    se <- sqrt(pmax(0,rowSums((contrast %*% covariance)*contrast)))

    runners <- switch(outcome, run=c(729,729), hr_mean=c(810,820),
                      steady_speed=c(660,669), c(820,820))
    t_wide <- qt(.975,runners[1]-1)
    t_narrow <- qt(.975,runners[2]-1)
    statistic <- abs(effect/se)
    statistic[se==0] <- NA_real_
    values <- data.table(outcome,temperature=grid,estimate=effect,se_model=se,
      lower_min=effect-t_wide*se,lower_max=effect-t_narrow*se,
      upper_min=effect+t_narrow*se,upper_max=effect+t_wide*se,
      p_min=2*pt(-statistic,runners[2]-1),p_max=2*pt(-statistic,runners[1]-1))
    cols <- c("estimate","lower_min","lower_max","upper_min","upper_max")
    if (outcome %in% c("log_duration","log_distance"))
      values[, (cols):=lapply(.SD,function(x) 100*expm1(x)),.SDcols=cols]
    if (outcome %in% c("start_before_09","start_after_19","session_day_percentile","run"))
      values[, (cols):=lapply(.SD,function(x) 100*x),.SDcols=cols]
    curve_rows[[i]] <- values
  }
  curves <- rbindlist(curve_rows)
  table2 <- curves[temperature %in% targets[2:3]]
  table2[,target:=ifelse(temperature==targets[2],"P90","P95")]
  primary <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean")
  table2[target=="P95" & outcome %in% primary,
         `:=`(p_holm_min=p.adjust(p_min,"holm"),p_holm_max=p.adjust(p_max,"holm"))]
  setcolorder(table2,c("outcome","target","temperature"))
  dir.create(file.path(public_out,"tables","manuscript"),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(public_out,"derived"),recursive=TRUE,showWarnings=FALSE)
  fwrite(table2,file.path(public_out,"tables","manuscript","table2.csv"))
  fwrite(curves,file.path(public_out,"derived","temperature_curves.csv"))

  # Validate against the printed manuscript only after finishing estimation.
  fields <- c("estimate","lower","upper","p","p_holm")
  reconstructed <- rbindlist(lapply(fields,function(field) {
    lo <- if (field=="estimate") field else paste0(field,"_min")
    hi <- if (field=="estimate") field else paste0(field,"_max")
    table2[,.(outcome,target,statistic=field,r_min=get(lo),r_max=get(hi))]
  }))
  printed <- fread(file.path(repo,"data","processed","primary_printed.csv"),
                   colClasses=c(printed="character"))[table=="T2"]
  check <- merge(printed,reconstructed,by=c("outcome","target","statistic"),all.x=TRUE)
  stopifnot(nrow(check)==81L,all(is.finite(check$r_min)),all(is.finite(check$r_max)))
  # Every admissible value must agree after rounding to the printed precision.
  check[,matches:=ifelse(upper_bound,r_max<reference,
    r_min>=reference-tolerance-1e-12 & r_max<=reference+tolerance+1e-12)]
  dir.create(file.path(public_out,"validation"),recursive=TRUE,showWarnings=FALSE)
  fwrite(check[,.(outcome,target,statistic,printed,r_min,r_max,matches)],
         file.path(public_out,"validation","table2_printed.csv"))
  if (!all(check$matches)) stop("Check table2_printed.csv for a rounding or estimation difference.")
  message("Table 2: all 81 printed values reproduced from whole-cohort model statistics.")
})
