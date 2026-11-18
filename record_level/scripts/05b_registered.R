# Tables S2, S3 and S29 use the models fitted in scripts 02 and 03.
# The prespecified models include calendar-year indicators.
registered <- list()
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
  registered[[length(registered)+1L]] <- save_contrast(fit,p95,"identity","S2",as.character(k),"fixed_hr","P95")
  if (k==4) phys_registered <- fit
}

# Temperature contrasts for panels A and B. Apply Holm's correction to the family of 30 versus 20 °C tests.
fits_registered <- list()
for (j in seq_len(nrow(outcomes))) {
  v <- outcomes$outcome[j]
  fit <- if (v=="run") main_fits$run else fit_fe(sessions[[v]],cbind(X_session,year_indicators),runner,runner)
  fits_registered[[v]] <- fit
  for (label in c("main","registered")) {
    m <- if (label=="main") main_fits[[v]] else fit
    registered[[length(registered)+1L]] <- save_contrast(m,d95,outcomes$transform[j],"S3",label,v,"P95")
  }
  for (target in names(anchor)) {
    diff <- rcs_basis(anchor[[target]],temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
    registered[[length(registered)+1L]] <- save_contrast(fit,diff,outcomes$transform[j],"S3","registered",v,target)
  }
}
for (label in c("main","registered")) {
  fit <- if (label=="main") hr_phfit else phys_registered
  registered[[length(registered)+1L]] <- save_contrast(fit,p95,"identity","S3",label,"fixed_hr","P95")
}
for (target in names(anchor)) {
  diff <- setNames(hr_spline_change(anchor[[target]],20,hr_ksession),c("T1","T2"))
  registered[[length(registered)+1L]] <- save_contrast(phys_registered,diff,"identity","S3","registered","fixed_hr",target)
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
  registered[[length(registered)+1L]] <- data.table(table="S3",model="weekend_diff",outcome=v,target="joint",
                                                    p=pf(wald,length(term),fit$df,lower.tail=FALSE))
  for (target in names(anchor)) {
    diff <- rcs_basis(anchor[[target]],temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
    names(diff) <- colnames(interaction)
    registered[[length(registered)+1L]] <- save_contrast(fit,diff,outcomes$transform[j],"S3","weekend_diff",v,target)
  }
}

# Table S3, panel C. The prespecified normal-approximation delta intervals use
# independent covariance blocks for the thermal gap, speed, heart rate and physiology.
b <- rcs_basis(20,temperature_knots)[1,]-colMeans(temperature_basis)
change <- rcs_basis(30,temperature_knots)[1,]-rcs_basis(20,temperature_knots)[1,]
fits_parts <- list(gap=fits_registered$heat_avoided_c,pace=fits_registered$gps_speed_median_m_s,
                   net=fits_registered$hr_mean,phys=phys_registered)
for (key in c("gap","pace","net")) {
  m <- fits_parts[[key]]
  idx <- match(c("linear","cubic1"),names(m$beta))
  names(m$beta)[idx] <- colnames(m$score)[idx] <- c("T1","T2")
  fits_parts[[key]] <- m
}
blocks <- hr_stack(fits_parts)
group <- sub("\\..*$","",names(blocks$beta))
blocks$vcov[outer(group,group,"!=")] <- 0
planned_parts <- function(beta) {
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
pc <- hr_delta(planned_parts,blocks$beta,blocks$vcov)
registered[[length(registered)+1L]] <- pc[,.(table="S3",model="decomposition",outcome=quantity,target="30v20",
                                              estimate,lower,upper,p)]
reg <- rbindlist(registered,fill=TRUE)
reg[,p_holm:=NA_real_]
reg[table=="S3" & model=="registered" & target=="30v20" & outcome %in% outcomes$outcome[1:4],p_holm:=p.adjust(p,"holm")]

# Table S29A. Compute exposure percentiles within the corresponding sample.
percentiles <- c(.01,.10,.50,.90,.95,.99)
dists <- list(day_max=sessions$day_max_c,day_mean=sessions$day_mean_c,session=sessions$temp_c_10km,
              dew=sessions$day_max_dew_c,fixed_hr=hr_phys$temp_c_10km)
r29 <- list()
for (key in names(dists)) {
  r29[[length(r29)+1L]] <- data.table(table="S29",model="percentile",outcome=key,
                                     target=paste0("P",100*percentiles),estimate=quantile(dists[[key]],percentiles,na.rm=TRUE,names=FALSE))
}
# Table S29B. Round runner counts down for disclosure when reporting temperature support.
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
    r29[[length(r29)+1L]] <- data.table(table="S29",model=names(rounded),outcome=key,
                                        target=c("P90","P95")[i-1L],estimate=as.numeric(rounded))
  }
}
fwrite(rbindlist(c(list(reg),r29),fill=TRUE),file.path(out,"supp_registered.csv"))
