# Figure 6 and Figures S3–S5. Scripts 04 and 06 provide the estimates and intervals.
# The published figures supply the plotting grid and the rules for displaying estimates.
local({
  blue <- "#1D596D"; orange <- "#A55029"
  shaded <- function(x,lower,upper,colour) {
    polygon(c(x,rev(x)),c(lower,rev(upper)),border=NA,col=adjustcolor(colour,.16))
  }
  pdf_png <- function(name,width,height,draw) {
    for (ext in c("pdf","png")) {
      path <- file.path(out,paste0(name,".",ext))
      if(ext=="pdf") grDevices::cairo_pdf(path,width=width,height=height)
      else grDevices::png(path,width=width,height=height,units="in",res=320,type="cairo")
      tryCatch(draw(),finally=dev.off())
    }
  }
  daily <- fread(file.path(out,"04_lag_profile.csv"))[target_label=="P95"]
  setorder(daily,day_relative_to_hot_day)
  cumulative <- daily[day_relative_to_hot_day>=0]
  sport <- fread(file.path(out,"04_sport.csv"))[target_label=="P95" & variant=="c3"]
  groups <- c(any_sport="Any sport",outdoor_running="Outdoor running",treadmill_running="Treadmill running",
              outdoor_cycling="Outdoor cycling",cycling_recreation="    Recreation",cycling_transport="    Transport",
              indoor_cycling="Indoor cycling",walking_hiking="Walking and hiking",indoor_training="Indoor training",
              swimming="Swimming",indoor_all="All indoor sessions")
  sport <- sport[match(names(groups),outcome)]
  stopifnot(!anyNA(sport$estimate),nrow(daily)==11,nrow(cumulative)==8)
  pdf_png("figure_6_participation",10.5,4.65,function() {
    par(mfrow=c(1,2),mar=c(4.5,4.3,3,0.8),las=1,mgp=c(2.6,.6,0),cex.axis=.86,bty="l")
    limits <- 100*range(daily$day_ci_low,daily$day_ci_high,cumulative$cumulative_ci_low,cumulative$cumulative_ci_high)
    plot(NA,xlim=c(-3.2,7.2),ylim=limits,xlab="Days from the hot day",ylab="",
         main="a   Days with a recorded run",xaxt="n",cex.main=1)
    mtext("A day counts if at least one recorded run has heart rate data.",side=3,line=.2,cex=.72)
    mtext("Change in runner-days per 100 windows",side=2,line=2.8,las=0,cex=.84)
    axis(1,at=-3:7)
    abline(h=0,col="grey60",lty=3)
    abline(v=0,col="grey75",lty=3)
    shaded(cumulative$day_relative_to_hot_day,100*cumulative$cumulative_ci_low,100*cumulative$cumulative_ci_high,blue)
    lines(cumulative$day_relative_to_hot_day,100*cumulative$cumulative_estimate,col=blue,lwd=2)
    with(daily,segments(day_relative_to_hot_day,100*day_ci_low,day_relative_to_hot_day,100*day_ci_high,col="grey38"))
    points(daily$day_relative_to_hot_day,100*daily$day_estimate,pch=21,col="grey25",bg="white",cex=.8)
    legend("topright",legend=c("Daily (percentage points)","Cumulative from day 0"),pch=c(21,NA),pt.bg="white",lty=c(NA,1),
           col=c("grey25",blue),bty="n",lwd=c(1,2),cex=.78)
    par(mar=c(4.5,8,3,.9))
    pos <- rev(seq_len(nrow(sport)))
    plot(NA,xlim=100*range(sport$ci_low,sport$ci_high),ylim=c(.5,length(pos)+.5),
         xlab="Change (percentage points)",ylab="",yaxt="n",main="b   Workouts on a hot day",cex.main=1)
    axis(2,at=pos,labels=unname(groups),tick=FALSE,cex.axis=.82)
    abline(v=0,col="grey60",lty=3)
    segments(100*sport$ci_low,pos,100*sport$ci_high,pos,col=blue,lwd=1.6)
    points(100*sport$estimate,pos,pch=19,col=blue,cex=.77)
  })

  # Figure S3. Panel a uses runner-days as its denominator. Panel b uses all
  # workouts recorded during follow-up.
  recording <- fread(file.path(out,"04_recording_curves.csv"))
  recording_at <- fread(file.path(out,"04_recording_contrasts.csv"))
  setorder(recording,model_id,exposure_c)
  recorded_models <- c("p_run","p_run_any_recorded","c_hr_recorded")
  med <- recording[abs(estimate)<1e-12,unique(exposure_c)]
  at <- recording_at[model_id=="p_run"]
  top_ticks <- c(med,at$exposure_c[match(c("P90","P95"),at$target_label)])
  pdf_png("figure_S3_participation",8.4,4.35,function() {
    par(mfrow=c(1,2),mar=c(4.3,4.1,3.6,.85),mgp=c(2.5,.55,0),cex.axis=.83,las=1,bty="l")
    for(j in 1:2) {
      ids <- if(j==1L) recorded_models[1:2] else recorded_models[3]
      d <- recording[model_id %in% ids]
      plot(NA,xlim=range(d$exposure_c),ylim=100*range(d$ci_low,d$ci_high),
           xlab="Daily maximum (°C)",ylab="")
      title(main=c("a   Recorded a run on the day","b   Heart rate available in a recorded run")[j],line=2,cex.main=.96)
      mtext("Change (percentage points)",2,line=2.7,las=0,cex=.87)
      axis(3,at=top_ticks,labels=c("P50","P90","P95"),cex.axis=.65,gap.axis=-1)
      abline(h=0,col="grey60",lty=3)
      for(k in seq_along(ids)) {
        one <- d[model_id==ids[k]]
        col <- c(blue,orange)[k]
        shaded(one$exposure_c,100*one$ci_low,100*one$ci_high,col)
        lines(one$exposure_c,100*one$estimate,col=col,lwd=1.8,lty=k)
      }
      if(j==1L) legend("bottom",c("Run with heart rate","Any recorded run"),
                       col=c(blue,orange),lwd=2,lty=1:2,cex=.79,bty="n")
    }
  })

  curves <- fread(file.path(out,"06_weekend_curves.csv"))
  setorder(curves,model_id,stratum,exposure_c)
  pdf_png("figure_S5_weekend",11,3.8,function() {
    par(mfrow=c(1,3),mar=c(4.3,4.1,3,.7),mgp=c(2.4,.6,0),cex.axis=.85,bty="l",las=1)
    titles <- c("a   Heat avoided","b   Whole-session heart rate","c   Weekend minus weekday")
    for(j in 1:3) {
      mid <- if(j==1L) "w_b1_heat_avoided" else "w_h1_hr_mean"
      groups <- if(j==3L) "weekend_minus_weekday" else c("weekday","weekend")
      d <- curves[model_id==mid & stratum %in% groups]
      plot(NA,xlim=range(d$exposure_c),ylim=range(d$ci_low,d$ci_high),xlab="Daily maximum (°C)",ylab="",
           main=titles[j],cex.main=.94)
      mtext(if(j==1L) "Temperature difference (°C)" else "Heart-rate difference (beats/min)",2,line=2.8,las=0,cex=.82)
      abline(h=0,col="grey60",lty=3)
      for(k in seq_along(groups)) {
        one <- d[stratum==groups[k]]
        colour <- c(blue,orange)[k]
        shaded(one$exposure_c,one$ci_low,one$ci_high,colour)
        lines(one$exposure_c,one$estimate,col=colour,lwd=1.8,lty=k)
      }
      if(j==1L) legend("bottomleft",c("Weekday","Weekend"),col=c(blue,orange),lty=1:2,lwd=2,bty="n",cex=.8)
    }
  })

  # Percentages in 1 °C bands among sessions with data for both exposure windows.
  # The denominator includes sessions outside the P1–P99 range.
  hist <- list()
  targets <- list()
  for(id in c("daily_max","session")) {
    dir <- if(id=="daily_max") "history/daily_maximum" else "history/session_temperature"
    x <- if(id=="daily_max") "day_max_c" else "temp_c_10km"
    tab <- read_fitrec(file.path(dir,"sessions_with_history.parquet"),
                       c("user_hash","history14_c","history7_c",x))
    tab <- tab[is.finite(history14_c) & is.finite(history7_c) & is.finite(get(x))]
    q <- quantile(tab[[x]],c(.01,.5,.9,.95,.99),names=FALSE)
    targets[[id]] <- q[2:4]
    n <- nrow(tab)
    obs <- tab[,.(n_sessions=.N,runners=uniqueN(user_hash)),by=.(bin_left_c=floor(get(x)))]
    band <- data.table(bin_left_c=seq(floor(q[1]),ceiling(q[5])-1))
    band <- obs[band,on="bin_left_c"]
    band[,percent:=fifelse(!is.na(runners) & runners>=10,100*n_sessions/n,NA_real_)]
    # Export a percentage only if it meets the reporting threshold.
    hist[[id]] <- band[,.(exposure_id=id,bin_left_c,percent)]
  }
  histogram <- rbindlist(hist)
  fwrite(histogram,file.path(out,"figure_S4_exposure_distribution.csv"))
  ref <- fread(input_file("references/history/distributions.csv"))
  compare_aggregate(histogram,ref,c("exposure_id","bin_left_c"),"percent","10b_history_histogram")
  curves <- fread(file.path(out,"06_hist_curve.csv"))
  display <- fread(input_file("references/heart_rate_and_history/history_curves.csv"),
                   select=c("model_id","target_c","display_allowed"))
  curves <- curves[display[display_allowed==TRUE],on=c("model_id","target_c"),nomatch=0L]
  setorder(curves,model_id,target_c)
  types <- c("B1","B2","B3","H1","START","H2")
  titles <- c("Heat avoided","Median speed","Session duration","Whole-session heart rate","Local start time","Heart rate at fixed speed")
  units <- c("°C","m/s","percent","beats/min","minutes","beats/min")
  pdf_png("figure_S4_prior_heat",11,7.1,function() {
    par(mfrow=c(2,3),mar=c(3.8,4.5,3.8,1),mgp=c(2.4,.55,0),cex.axis=.79,bty="l",las=1)
    for(j in seq_along(types)) {
      type <- types[j]
      id <- if(type=="H2") "session" else "daily_max"
      models <- c(if(type=="H2") "AH14" else paste0("AB14_",type),paste0("EW14_",type))
      d <- copy(curves[model_id %in% models])
      if(type=="B3") for(k in c("ci_low","estimate","ci_high")) set(d,j=k,value=100*expm1(d[[k]]))
      y <- range(c(d$ci_low,d$ci_high))
      width <- diff(y)
      bottom <- y[1]-.24*width
      plot(NA,xlim=range(d$target_c),ylim=c(bottom,y[2]),xlab=if(id=="session") "Session temperature (°C)" else "Daily maximum (°C)",ylab="")
      title(main=paste0(letters[j],"   ",titles[j]),line=2.1,cex.main=.94)
      mtext(paste0("Modification (",units[j],")"),side=2,line=3.05,las=0,cex=.8)
      bins <- histogram[exposure_id==id & !is.na(percent)]
      rect(bins$bin_left_c,bottom,bins$bin_left_c+1,bottom+.18*width*bins$percent/max(bins$percent),
           col="grey86",border=NA)
      axis(3,at=targets[[id]],labels=c("P50","P90","P95"),cex.axis=.69,tck=-.018,padj=.6)
      abline(h=0,col="grey65",lty=3)
      abline(v=targets[[id]][1L],col="grey65",lty=3)
      for(k in 2:1) {
        one <- d[model_id==models[k]]
        col <- c(blue,orange)[k]
        if(k==2L) shaded(one$target_c,one$ci_low,one$ci_high,col)
        for(v in c("ci_low","ci_high")) lines(one$target_c,one[[v]],col=adjustcolor(col,.62),lty=k,lwd=.65)
        lines(one$target_c,one$estimate,col=col,lty=k,lwd=1.8)
      }
      if(j==1L) legend("bottomleft",legend=c("Equal-day","Weighted"),lty=1:2,lwd=2,col=c(blue,orange),cex=.65,
                       bg="white",box.lty=0)
    }
  })
})
