# Tables S3, S4 and S12 use the models fitted in scripts 02 and 03.
# Compare nested adjustments and report the temperature range for each sample.
model_rows <- list()
yr <- function(t) factor(format(as.POSIXct(t,origin="1970-01-01",tz="UTC"),"%Y",tz="UTC"))
pm_year <- model.matrix(~yr(hr_phys$canonical_start_utc))[,-1,drop=FALSE]
p95 <- setNames(hr_spline_change(hr_stargets["P95"],hr_stargets["P50"],hr_ksession),c("T1","T2"))
anchor <- c("25v20"=25,"30v20"=30)

# Fit the four nested models for heart rate at a fixed speed.
for (k in 1:4) {
  X <- switch(k,cbind(intercept=1,hr_xphys[,c("T1","T2")]),cbind(intercept=1,hr_xphys),
              hr_xphys,cbind(hr_xphys,pm_year))
  fe <- if (k>=3) hr_phys$user_hash else NULL
  fit <- fit_fe(hr_phys$hr_lag60,X,fe,hr_phys$user_hash)
  model_rows[[length(model_rows)+1L]] <- save_contrast(fit,p95,"identity","S3",as.character(k),"fixed_hr","P95")
  if (k==4) phys_with_year <- fit
}

# The reference tables call the year-adjusted model "registered".
# Temperature contrasts for panels A and B. Apply Holm's correction to the family of 30 versus 20 °C tests.
fits_with_year <- list()
for (j in seq_len(nrow(outcomes))) {
  v <- outcomes$outcome[j]
  fit <- if (v=="run") main_fits$run else fit_fe(sessions[[v]],cbind(X_session,year_indicators),runner,runner)
  fits_with_year[[v]] <- fit
  for (label in c("main","registered")) {
    m <- if (label=="main") main_fits[[v]] else fit
    model_rows[[length(model_rows)+1L]] <- save_contrast(m,d95,outcomes$transform[j],"S12",label,v,"P95")
  }
  for (target in names(anchor)) {
    diff <- rcs_basis(anchor[[target]],temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
    model_rows[[length(model_rows)+1L]] <- save_contrast(fit,diff,outcomes$transform[j],"S12","registered",v,target)
  }
}
for (label in c("main","registered")) {
  fit <- if (label=="main") hr_phfit else phys_with_year
  model_rows[[length(model_rows)+1L]] <- save_contrast(fit,p95,"identity","S12",label,"fixed_hr","P95")
}
for (target in names(anchor)) {
  diff <- setNames(hr_spline_change(anchor[[target]],20,hr_ksession),c("T1","T2"))
  model_rows[[length(model_rows)+1L]] <- save_contrast(phys_with_year,diff,"identity","S12","registered","fixed_hr",target)
}
weekend <- as.numeric(sessions$weekend)
interaction <- temperature_basis*weekend
colnames(interaction) <- paste0("weekend_",colnames(temperature_basis))
for (j in 1:4) {
  v <- outcomes$outcome[j]
  fit <- fit_fe(sessions[[v]],cbind(temperature_basis,season,year_indicators,weekend,interaction),runner,runner)
  term <- colnames(interaction)
  coef <- fit$beta[term]
  wald <- as.numeric(crossprod(coef,solve(fit$vcov[term,term],coef))) / length(term)
  model_rows[[length(model_rows)+1L]] <- data.table(table="S12",model="weekend_diff",outcome=v,target="joint",
                                                    p=pf(wald,length(term),fit$df,lower.tail=FALSE))
  for (target in names(anchor)) {
    diff <- rcs_basis(anchor[[target]],temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
    names(diff) <- colnames(interaction)
    model_rows[[length(model_rows)+1L]] <- save_contrast(fit,diff,outcomes$transform[j],"S12","weekend_diff",v,target)
  }
}

# Table S12, panel C. The normal-approximation delta intervals use
# independent covariance blocks for the thermal gap, speed, heart rate and physiology.
b <- rcs_basis(20,temperature_knots)[1,]-colMeans(temperature_basis)
change <- rcs_basis(30,temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
fits_parts <- list(gap=fits_with_year$heat_avoided_c,pace=fits_with_year$gps_speed_median_m_s,
                   net=fits_with_year$hr_mean,phys=phys_with_year)
for (key in c("gap","pace","net")) {
  m <- fits_parts[[key]]
  idx <- match(c("linear","cubic1"),names(m$beta))
  names(m$beta)[idx] <- colnames(m$score)[idx] <- c("T1","T2")
  fits_parts[[key]] <- m
}
blocks <- hr_stack(fits_parts)
group <- sub("\\..*$","",names(blocks$beta))
blocks$vcov[outer(group,group,"!=")] <- 0
decomposition_terms <- function(beta) {
  gap <- beta[c("gap.T1","gap.T2")]
  gap0 <- mean(sessions$heat_avoided_c)+sum(b*gap)
  speed <- sum(change*beta[c("pace.T1","pace.T2")])
  g <- sum(change*gap)
  net <- sum(change*beta[c("net.T1","net.T2")])
  thermal <- beta[c("phys.T1","phys.T2")]
  bench <- sum(thermal*hr_spline_change(30+gap0,20+gap0,hr_ksession))
  timed <- sum(thermal*hr_spline_change(30+gap0+g,20+gap0,hr_ksession))
  slow <- sum(beta[c("phys.V1","phys.V2")]*hr_spline_change(hr_speed0+speed,hr_speed0,hr_kspeed,TRUE))
  c(ref_temp=20+gap0,gap=g,speed=speed,benchmark=bench,timing=timed,slower=slow,
    expected=timed+slow,net=net,offset_timing=1-timed/bench,offset_total=1-net/bench)
}
pc <- hr_delta(decomposition_terms,blocks$beta,blocks$vcov)
model_rows[[length(model_rows)+1L]] <- pc[,.(table="S12",model="decomposition",outcome=quantity,target="30v20",
                                              estimate,lower,upper,p)]
model_contrasts <- rbindlist(model_rows,fill=TRUE)
model_contrasts[,p_holm:=NA_real_]
model_contrasts[table=="S12" & model=="registered" & target=="30v20" & outcome %in% outcomes$outcome[1:4],p_holm:=p.adjust(p,"holm")]

# Table S4A. Compute exposure percentiles within the corresponding sample.
percentiles <- c(.01,.10,.50,.90,.95,.99)
dists <- list(day_max=sessions$day_max_c,day_mean=sessions$day_mean_c,session=sessions$temp_c_10km,
              dew=sessions$day_max_dew_c,fixed_hr=hr_phys$temp_c_10km)
support_rows <- list()
for (key in names(dists)) {
  support_rows[[length(support_rows)+1L]] <- data.table(table="S4",model="percentile",outcome=key,
                                     target=paste0("P",100*percentiles),estimate=quantile(dists[[key]],percentiles,na.rm=TRUE,names=FALSE))
}
# Table S4B. Round runner counts down for disclosure when reporting temperature support.
sample_support <- list(main=sessions,mean_hr=sessions[is.finite(hr_mean)],steady=sessions[is.finite(steady_speed)],
                       panel=days,fixed_hr=hr_phys)
for (key in names(sample_support)) {
  group <- sample_support[[key]]
  temp <- if (key=="fixed_hr") group$temp_c_10km else group$day_max_c
  targets <- if (key=="fixed_hr") hr_stargets else temperature_targets
  for (i in 2:3) {
    flags <- data.table(id=group$user_hash,ref=abs(temp-targets[1])<=2,high=abs(temp-targets[i])<=2)
    counts <- flags[,.(ref=any(ref),high=any(high)),by=id][,.(near_ref=sum(ref),near_target=sum(high),both=sum(ref&high))]
    rounded <- floor(unlist(counts)/10)*10
    support_rows[[length(support_rows)+1L]] <- data.table(table="S4",model=names(rounded),outcome=key,
                                        target=c("P90","P95")[i-1L],estimate=as.numeric(rounded))
  }
}
fwrite(rbindlist(c(list(model_contrasts),support_rows),fill=TRUE),file.path(out,"supp_registered.csv"))
