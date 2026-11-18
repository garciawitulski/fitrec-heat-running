# Draw the manuscript and supplementary figures from their aggregate data.
# The study design and map are provided as rendered images.
# Some confidence bands use the published lower bound for the runner count.
# The figure source guide identifies them and describes their precision.
local({
  source_data <- file.path(repo,"data","published","figures")
  manuscript <- file.path(public_out,"figures","manuscript")
  supplement <- file.path(public_out,"figures","supplementary")
  dir.create(manuscript,recursive=TRUE,showWarnings=FALSE)
  dir.create(supplement,recursive=TRUE,showWarnings=FALSE)
  read_series <- function(name) fread(file.path(source_data,paste0(name,".csv")))
  colours <- c("#1F5974","#B65A34","#397A68","#785380","#383B3D")
  # Open the device for one figure and close it even if a panel fails.
  draw_pdf <- function(path,width,height,draw) {
    grDevices::pdf(path,width=width,height=height,useDingbats=FALSE)
    tryCatch({par(mgp=c(2.4,.6,0),bty="l",cex.axis=.86); draw()},finally=grDevices::dev.off())
  }
  ribbon <- function(x,low,high,col) {
    polygon(c(x,rev(x)),c(low,rev(high)),border=NA,col=adjustcolor(col,.15))
  }
  forest <- function(x,labels,title,xlab) {
    at <- rev(seq_len(nrow(x)))
    plot(NA,xlim=range(c(0,x$ci_low,x$ci_high)),ylim=c(.5,nrow(x)+.5),
         yaxt="n",ylab="",main=title,xlab=xlab,cex.main=.92)
    axis(2,at=at,labels=labels,las=1,tick=FALSE,cex.axis=.8)
    abline(v=0,lty=3,col="grey55")
    segments(x$ci_low,at,x$ci_high,at,col=colours[1],lwd=1.5)
    points(x$estimate,at,pch=19,col=colours[1],cex=.8)
  }
  curves <- function(x,ids,labels,title,ylab,xlab="Daily maximum (degrees C)",key="topleft") {
    x <- x[series %in% ids]
    stopifnot(nrow(x)>0L)
    plot(range(x$temperature),range(x$ci_low,x$ci_high),type="n",main=title,ylab=ylab,xlab=xlab,cex.main=.92)
    abline(h=0,lty=3,col="grey60")
    for (j in seq_along(ids)) {
      z <- x[series==ids[j]][order(temperature)]
      ribbon(z$temperature,z$ci_low,z$ci_high,colours[j])
      lines(z$temperature,z$estimate,col=colours[j],lty=j,lwd=1.7)
    }
    if(length(ids)>1L) legend(key,labels,col=colours[seq_along(ids)],lty=seq_along(ids),bty="n",cex=.70)
  }
  stopifnot(file.copy(file.path(source_data,"fig1_design.png"),file.path(manuscript,"figure_1_design.png"),overwrite=TRUE),
            file.copy(file.path(source_data,"supp_fig1_map.png"),file.path(supplement,"figure_S1_map.png"),overwrite=TRUE))

  paired <- read_series("figure2_paired")
  draw_pdf(file.path(manuscript,"figure_2_paired.pdf"),8.3,6.0,function() {
    par(mfrow=c(2,1),mar=c(3.7,13,2.1,1))
    keys <- c("steady_fixed_speed","steady_speed_free","whole_session_mean",
              "steady_speed_free_minus_whole_session","steady_fixed_speed_minus_whole_session")
    labels <- c("Steady, fixed speed","Steady, speed free","Whole session","Speed free minus whole session","Fixed speed minus whole session")
    for (j in 1:2) {
      i <- if(j==1) 1:3 else 4:5
      z <- paired[quantity %in% keys[i]]
      at <- rev(seq_along(i))
      plot(NA,xlim=range(c(0,z$ci_low,z$ci_high)),ylim=c(.4,length(i)+.6),yaxt="n",ylab="",
           main=c("a  Heart rate in the same runs","b  Within-run differences")[j],xlab="Heart-rate change (beats/min)")
      axis(2,at=at,labels=labels[i],las=1,tick=FALSE)
      abline(v=0,lty=3,col="grey60")
      for(k in 1:2) {
        v <- z[target_label==c("P90","P95")[k]][match(keys[i],quantity)]
        y <- at+c(.12,-.12)[k]
        segments(v$ci_low,y,v$ci_high,y,col=colours[k],lwd=1.5)
        points(v$estimate,y,pch=c(1,19)[k],col=colours[k])
      }
      legend("bottomright",c("P90","P95"),pch=c(1,19),col=colours[1:2],bty="n",horiz=TRUE,cex=.78)
    }
  })

  hours <- read_series("figure3_hourly")
  bins <- read_series("figure3_daily_temperature")
  days <- read_series("figure3_daily_participation")
  bands <- c("le10","10to20","20to25","25to30","gt30")
  labels <- c("<=10","10-20","20-25","25-30",">30")
  draw_pdf(file.path(manuscript,"figure_3_descriptive.pdf"),10.1,7,function() {
    par(mfrow=c(2,2),mar=c(4.0,4.2,2.2,1))
    for (k in 1:2) {
      v <- c("temperature","share")[k]
      plot(0:23,rep(0,24),type="n",ylim=range(hours[[v]]),xlab="Hour of local day",
           main=c("a  Local temperature","b  Run start time")[k],
           ylab=c("Hourly temperature (degrees C)","Share of runs starting")[k])
      for(j in seq_along(bands)) {
        z <- hours[band==bands[j]][order(hour)]
        lines(z$hour,z[[v]],col=c("#173D70","#476F9F","#538173","#CB6C46","#8A3030")[j],lwd=1.8)
      }
      if(k==1) legend("topleft",labels,title="Daily maximum (degrees C)",
        col=c("#173D70","#476F9F","#538173","#CB6C46","#8A3030"),lty=1,bty="n",cex=.67)
    }
    plot(bins$mean_daily_max_c,bins$mean_session_temp_c,type="b",pch=16,main="c  Temperature during the run",
         xlab="Mean daily maximum (degrees C)",ylab="Mean session temperature (degrees C)")
    abline(0,1,lty=3,col="grey60")
    plot(range(days$bin_low_c)+1,range(days$run_share),type="n",main="d  Recorded runs",
         xlab="Daily maximum (degrees C)",ylab="Share of days with a run")
    for(k in 0:1) {
      z <- days[weekend==as.logical(k)][order(bin_low_c)]
      lines(z$bin_low_c+1,z$run_share,type="b",pch=16,col=colours[k+1L])
    }
    legend("topright",c("Weekdays","Weekends"),col=colours[1:2],lty=1,bty="n",cex=.78)
  })
  primary <- read_series("figure4_primary")
  fixed <- read_series("figureS2_curves")
  ids <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean")
  titles <- c("Heat avoided","Median speed","Duration","Whole-session mean heart rate")
  units <- c("degrees C","m/s","%","beats/min")
  draw_pdf(file.path(manuscript,"figure_4_primary.pdf"),10,7,function() {
    par(mfrow=c(2,2),mar=c(4.0,4.2,2.2,1))
    for(j in 1:4) {
      z <- primary[series==ids[j]][order(temperature)]
      plot(z$temperature,z$estimate,type="n",ylim=range(z$ci_low,z$ci_high),xlab="Daily maximum (degrees C)",
           ylab=paste0("Change (",units[j],")"),main=paste(letters[j],titles[j]))
      u <- par("usr")
      rect(bins$bin_low_c,u[3],bins$bin_high_c,u[3]+.12*diff(u[3:4])*bins$sessions/max(bins$sessions),col="grey88",border=NA)
      abline(h=0,lty=3,col="grey60")
      ribbon(z$temperature,z$ci_low,z$ci_high,colours[j]); lines(z$temperature,z$estimate,col=colours[j],lwd=2)
      if(j==1) with(fixed[series=="heat_avoided_change_fixed_hour_c"][order(temperature)],lines(temperature,estimate,lty=2))
      pts <- z[target %in% c("P90","P95")]
      segments(pts$temperature,pts$ci_low,pts$temperature,pts$ci_high,col=colours[j])
      points(pts$temperature,pts$estimate,pch=21,bg="white",col=colours[j])
    }
  })
  stages <- read_series("figure5_stages")
  heat <- read_series("figure5_stages_curve")
  checks <- read_series("figure5_checks")
  stage_labels <- c("Initial benchmark","Fixed hour","Observed timing","Timing and speed","Measured heart rate")
  draw_pdf(file.path(manuscript,"figure_5_decomposition.pdf"),11.8,8.4,function() {
    layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(1,1.25))
    par(mar=c(6.7,4,2.3,1))
    level <- stages$estimate
    lower <- c(0,level[1],0,level[2],level[3],0,level[4],0)
    upper <- level[c(1,2,2,3,4,4,5,5)]
    positions <- seq_along(upper)
    plot(NA,xlim=c(.4,8.6),ylim=range(0,stages$ci_low,1.15*stages$ci_high),
      xaxt="n",xlab="",ylab="Heart-rate change (beats/min)",main="a  Accounting for the P95 contrast")
    abline(h=0,col="grey70")
    segments(positions[-8]+.28,upper[-8],positions[-1]-.28,upper[-8],col="grey65",lty=3)
    rect(positions-.28,pmin(lower,upper),positions+.28,pmax(lower,upper),
         col=c(colours[1],"grey75",colours[2],"grey75","grey75",colours[3],"white",colours[5]))
    absolute <- c(1,3,6,8)
    interval <- stages[c(1,2,4,5)]
    segments(absolute,interval$ci_low,absolute,interval$ci_high)
    text(absolute,interval$ci_high,formatC(interval$estimate,digits=2,format="f"),pos=3,cex=.7)
    steps <- c(2,4,5,7)
    text(steps,lower[steps],formatC(upper[steps]-lower[steps],digits=2,format="f"),pos=3,cex=.7)
    labels <- c("Original benchmark","Wider daily range","Fixed-hour benchmark","Earlier timing",
                "Slower speed","Timing and speed","Unexplained","Observed")
    text(positions,par("usr")[3],labels,srt=38,adj=1,xpd=NA,cex=.70)
    par(mar=c(4.0,4.2,2.3,1))
    curves(heat,stages$series,stage_labels,"b  Days warmer than the median","Change (beats/min)")
    par(mar=c(3.7,16,2.3,1))
    forest(checks,checks$label,"c  Alternative measurements and specifications","P95 versus P50 (beats/min)")
  })
  draw_pdf(file.path(supplement,"figure_S2_fixed_hour.pdf"),11.4,4.1,function() {
    par(mfrow=c(1,3),mar=c(4.0,4.2,2.2,1))
    curves(fixed,c("heat_avoided_change_c","heat_avoided_change_fixed_hour_c","heat_avoided_change_behavioural_c"),
           c("Observed","Fixed hour","Timing component"),"a  Heat avoided","Change (degrees C)",key="bottomleft")
    curves(fixed,"start_time_minutes","","b  Local start time","Change (min)")
    curves(fixed,"within_day_temperature_rank","","c  Within-day temperature rank","Change (percentile points)")
  })
  recording <- read_series("figureS3_recording")
  draw_pdf(file.path(supplement,"figure_S3_recording.pdf"),8.6,4.3,function() {
    par(mfrow=c(1,2),mar=c(4.0,4.2,2.8,1))
    curves(recording,c("run_with_heart_rate","any_recorded_run"),c("Run with heart rate","Any recorded run"),
          "a  Recorded a run on the day","Change (percentage points)",key="bottom")
    curves(recording,"heart_rate_recorded_during_workout","","b  Heart rate recorded during a workout","Change (percentage points)")
  })
  history <- read_series("figureS4_history")
  hist <- read_series("figureS4_distribution")
  ticks <- read_series("figureS4_temperatures")
  draw_pdf(file.path(supplement,"figure_S4_history.pdf"),11.2,7.4,function() {
    par(mfrow=c(2,3),mar=c(4,4.1,4.0,1))
    outcomes <- c(ids,"start_time_minutes","fixed_speed_heart_rate")
    for(j in seq_along(outcomes)) {
      id <- if(j==6) "session" else "daily_max"
      linesets <- paste0(c("equal_day_14d_","weighted_14d_"),outcomes[j])
      d <- history[series %in% linesets]
      y <- range(d$ci_low,d$ci_high)
      base <- y[1]-.22*diff(y)
      plot(range(d$temperature),c(base,y[2]),type="n",xlab=if(j==6) "Session temperature (degrees C)" else "Daily maximum (degrees C)",
           ylab=paste0("Modification (",c(units,"min","beats/min")[j],")"))
      title(main=paste(letters[j],c(titles,"Local start time","Heart rate at fixed speed")[j]),line=2.3,cex.main=.88)
      axis(3,at=ticks[exposure_id==id,temperature_c],labels=paste0("P",ticks[exposure_id==id,percentile]),cex.axis=.7)
      h <- hist[exposure_id==id];rect(h$bin_left_c,base,h$bin_right_c,base+.17*diff(y)*h$percent/max(h$percent),col="grey86",border=NA)
      abline(h=0,col="grey60",lty=3)
      for(k in 2:1) {
        z <- d[series==linesets[k]][order(temperature)]
        if(k==2) ribbon(z$temperature,z$ci_low,z$ci_high,colours[k])
        lines(z$temperature,z$ci_low,col=colours[k],lty=k,lwd=.55)
        lines(z$temperature,z$ci_high,col=colours[k],lty=k,lwd=.55)
        lines(z$temperature,z$estimate,col=colours[k],lty=k,lwd=1.7)
      }
      if(j==1) legend("bottomleft",c("Equal-day","Weighted"),col=colours[1:2],lty=1:2,bty="n",cex=.73)
    }
  })
  weekend <- read_series("figureS5_weekends")
  draw_pdf(file.path(supplement,"figure_S5_weekend.pdf"),10.8,3.9,function() {
    par(mfrow=c(1,3),mar=c(4,4.3,2.4,1))
    for(j in 1:3) {
      id <- if(j==1) "heat_avoided_c" else "hr_mean"
      d <- copy(weekend[series==id]);d[,series:=stratum]
      curves(d,if(j==3) "weekend_minus_weekday" else c("weekday","weekend"),c("Weekday","Weekend"),
             c("a  Heat avoided","b  Whole-session heart rate","c  Weekend minus weekday")[j],
             c("Change (degrees C)","Change (beats/min)","Difference (beats/min)")[j],key="bottomleft")
    }
  })
  sensitivity <- read_series("figureS6_sensitivity")
  draw_pdf(file.path(supplement,"figure_S6_sensitivity.pdf"),11.1,8,function() {
    par(mfrow=c(2,2),mar=c(4.1,12.5,2.2,1))
    labels <- c("Main specification","Daily mean temperature","ERA5 reanalysis","Weather cell within 20 km",
      "Two-way clustering","Runner by year-month effects","Maximum dew point","Sessions of 20-180 min",
      "At least 20 sessions","Hinge at 20 degrees C","Four-knot spline","Prior 14-day running time")
    for(j in 1:4) forest(sensitivity[outcome==ids[j]],labels,paste(letters[j],titles[j]),paste0("P95 vs P50 (",units[j],")"))
  })
})
