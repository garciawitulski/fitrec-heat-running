# Figure 4 and Supplementary Figure 2. The fixed-hour model adds indicators for
# local start hour. The uncertainty in the difference is based on the joint
# covariance of the two models, clustered by runner.
local({
  daily <- read_fitrec("panels/sessions.parquet",
    c("day_max_c","heat_avoided_c","user_hash","start_minutes_local","canonical_start_utc",
      "geo_mean_y","local_weekday","day_min_c"))[is.finite(day_max_c)]
  knots <- quantile(daily$day_max_c,c(.1,.5,.9),names=FALSE)
  targets <- quantile(daily$day_max_c,c(.5,.9,.95),names=FALSE)
  window <- quantile(daily$day_max_c,c(.01,.99),names=FALSE)
  grid <- sort(unique(c(seq(ceiling(2*window[1])/2,floor(2*window[2])/2,by=.5),targets)))
  wd <- model.matrix(~factor(local_weekday,levels=0:6),daily)[,-1,drop=FALSE]
  x <- cbind(rcs_basis(daily$day_max_c,knots),season_terms(daily$canonical_start_utc,daily$geo_mean_y<0),wd)
  hour <- model.matrix(~factor(floor(start_minutes_local/60)),daily)[,-1,drop=FALSE]
  free <- fit_fe(daily$heat_avoided_c,x,daily$user_hash,daily$user_hash)
  fixed <- fit_fe(daily$heat_avoided_c,cbind(x,hour),daily$user_hash,daily$user_hash)
  basis <- rcs_basis(grid,knots)
  ref <- rcs_basis(targets[1],knots)
  fixed_hours <- rbindlist(lapply(seq_along(grid),function(i) {
    contrast <- basis[i,]-ref[1,]
    a <- lin_contrast(free,contrast)
    b <- lin_contrast(fixed,contrast)
    influence <- (free$score[,names(contrast)]-fixed$score[,names(contrast)]) %*% contrast
    se <- sqrt(sum(influence^2))
    delta <- a$estimate-b$estimate
    gap <- qt(.975,free$df)*se
    timing <- data.table(estimate=delta,lower=delta-gap,upper=delta+gap,
      p=if(se==0) NA_real_ else 2*pt(-abs(delta/se),free$df))
    cbind(data.table(temperature=grid[i],quantity=c("heat_avoided_change_c","heat_avoided_change_fixed_hour_c",
      "heat_avoided_change_behavioural_c")),rbind(a,b,timing))
  }))
  fwrite(fixed_hours,file.path(out,"fixed_hour_curves.csv"))
  workdir <- "references/main"
  what <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean")
  paper <- fread(input_file(workdir,"curves.csv"))[outcome %in% what]
  calc <- fread(file.path(out,"session_curves.csv"))[outcome %in% what]
  # The source file stores log duration. The figure reports percentage changes.
  for (item in c("estimate","ci_low","ci_high"))
    paper[outcome=="log_duration",(item):=100*expm1(get(item))]
  paper[,temperature:=round(exposure_c,7)]
  calc[,temperature:=round(temperature,7)]
  setnames(paper,c("ci_low","ci_high","p_value"),c("lower","upper","p"))
  measures <- c("estimate","lower","upper","p")
  checks <- list(compare_aggregate(calc,paper,c("outcome","temperature"),measures,"09c_fig4"))
  paper_fixed <- fread(input_file(workdir,"fixed_hour_curves.csv"))
  setnames(paper_fixed,c("ci_low","ci_high","p_value"),c("lower","upper","p"))
  paper_fixed[,temperature:=round(exposure_c,7)]
  fixed_hours[,temperature:=round(temperature,7)]
  checks[[2]] <- compare_aggregate(fixed_hours,paper_fixed,c("quantity","temperature"),measures,"09c_fixed_hour")
  secondary <- fread(file.path(out,"session_curves.csv"))[outcome %in% c("start_minutes_local","session_day_percentile")]
  source_secondary <- fread(input_file(workdir,"secondary_outcome_curves.csv"))
  for (item in c("estimate","ci_low","ci_high"))
    source_secondary[outcome=="session_day_percentile",(item):=100*get(item)]
  setnames(source_secondary,c("ci_low","ci_high","p_value"),c("lower","upper","p"))
  secondary[,temperature:=round(temperature,7)]
  source_secondary[,temperature:=round(exposure_c,7)]
  compare_aggregate(secondary,source_secondary,c("outcome","temperature"),measures,"09c_secondary_hours")

  ylabels <- c("Session minus maximum temperature (°C)","Median running speed (m/s)",
               "Run duration (%)","Mean heart rate (beats/min)")
  cols <- c("#1F5974","#2A715A","#8A6334","#A14346")
  pdf(file.path(out,"figure4_primary_outcomes.pdf"),width=10,height=7,useDingbats=FALSE)
  tryCatch({
    par(mfrow=c(2,2),mar=c(4,4.3,2,1))
    counts <- hist(daily$day_max_c,breaks=seq(floor(min(daily$day_max_c)/2)*2,
      ceiling(max(daily$day_max_c)/2)*2,by=2),plot=FALSE)
    for (k in seq_along(what)) {
      z <- calc[outcome==what[k]][order(temperature)]
      plot(z$temperature,z$estimate,type="n",ylim=range(z$lower,z$upper),
        xlab="Daily maximum temperature (°C)",ylab=ylabels[k],
        main=paste0(letters[k],". ",c("Heat avoided","Speed","Duration","Heart rate")[k]))
      u <- par("usr")
      rect(head(counts$breaks,-1),u[3],tail(counts$breaks,-1),
        u[3]+.13*diff(u[3:4])*counts$counts/max(counts$counts),col="grey89",border=NA)
      abline(h=0,lty=3)
      polygon(c(z$temperature,rev(z$temperature)),c(z$lower,rev(z$upper)),
        border=NA,col=adjustcolor(cols[k],.23))
      lines(z$temperature,z$estimate,col=cols[k],lwd=2)
      if (k==1) with(fixed_hours[quantity=="heat_avoided_change_fixed_hour_c"],
        lines(temperature,estimate,lty=2))
      t <- z[temperature %in% round(targets[2:3],7)]
      with(t,segments(temperature,lower,temperature,upper,col=cols[k]))
      with(t,points(temperature,estimate,pch=21,bg="white",col=cols[k]))
    }
  },finally=dev.off())
  pdf(file.path(out,"figure_S2_fixed_hour.pdf"),width=11,height=4.1,useDingbats=FALSE)
  tryCatch({
    par(mfrow=c(1,3),mar=c(4.1,4.3,2,1))
    parts <- list(fixed_hours,secondary[outcome=="start_minutes_local"],secondary[outcome=="session_day_percentile"])
    for (j in 1:3) {
      a <- parts[[j]]
      plot(range(a$temperature),range(a$lower,a$upper),type="n",xlab="Daily maximum temperature (°C)",
        ylab=c("Session minus maximum temperature (°C)","Local start time (min)","Session temperature rank (points)")[j],
        main=c("a. Heat avoided","b. Start time","c. Within-day rank")[j])
      abline(h=0,lty=3)
      linesets <- if(j==1) split(a,a$quantity) else list(a)
      for (k in seq_along(linesets)) {
        z <- linesets[[k]][order(temperature)]
        polygon(c(z$temperature,rev(z$temperature)),c(z$lower,rev(z$upper)),border=NA,col=adjustcolor(cols[k],.16))
        lines(z$temperature,z$estimate,col=cols[k],lty=k,lwd=1.7)
      }
      if(j==1) legend("bottomleft",c("Timing component","Observed","Fixed hour"),col=cols[1:3],lty=1:3,bty="n",cex=.75)
    }
  },finally=dev.off())
})
