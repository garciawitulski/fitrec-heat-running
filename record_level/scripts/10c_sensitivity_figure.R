# Figure S6. The 12 rows in each panel show the main model and the
# eleven alternative specifications in the supplement.
local({
  est <- fread(file.path(out,"supp_session.csv"))[table=="S27" & target=="P95"]
  specs <- c("main","daily_mean","era5","20km","two_way","runner_month","dew",
             "duration_20_180","at_least_20","hinge","four_knots","history")
  labels <- c("Main specification","Daily mean temperature","ERA5 0.25° reanalysis",
              "Weather cell within 20 km","Clustering by runner and date","Runner × year-month effects",
              "Maximum dew point","Sessions of 20–180 minutes","Runners with at least 20 sessions",
              "Hinge at 20 °C","Four-knot spline","Prior 14-day running time")
  vars <- c("heat_avoided_c","gps_speed_median_m_s","log_duration","hr_mean")
  titles <- c("a   Heat avoided","b   Median speed","c   Duration","d   Whole-session mean heart rate")
  units <- c("Temperature difference (°C)","Speed difference (m/s)","Change in duration (%)",
             "Heart-rate difference (beats/min)")
  digits <- c(2,3,1,2)
  y <- rev(seq_along(specs))
  for (ext in c("pdf","png","svg")) {
    path <- file.path(out,paste0("figure_S6_sensitivity.",ext))
    switch(ext,pdf=cairo_pdf(path,width=13.5,height=7.6),
           png=png(path,width=13.5,height=7.6,units="in",res=320,type="cairo"),
           svg=svg(path,width=13.5,height=7.6))
    tryCatch({
      layout(matrix(1:10,nrow=2,byrow=TRUE),widths=c(2.5,2.4,1.85,2.4,1.85))
      par(cex=.95,oma=c(.2,.25,.25,.2),family="sans",fg="#263438",col.axis="#263438",col.lab="#263438")
      for (j in seq_along(vars)) {
        rows <- est[outcome==vars[j]][match(specs,model)]
        stopifnot(nrow(rows)==length(specs),all(is.finite(rows$estimate)))
        par(mar=c(3.3,.15,2,.15))
        if (j %% 2L==1L) {
          plot.new()
          plot.window(c(0,1),c(.4,12.6))
          text(.98,y,labels,adj=1,cex=.77,col="#263438")
        }
        limits <- range(c(0,rows$lower,rows$upper))
        tick <- pretty(limits,n=4)
        plot.new()
        plot.window(range(tick),c(.4,12.6))
        abline(v=0,col="#A8ADAD",lwd=.7)
        abline(v=rows$estimate[1],lty=2,col="#437184",lwd=.8)
        segments(rows$lower,y,rows$upper,y,col="#3A5966",lwd=1.6)
        points(rows$estimate,y,pch=c(18,rep(16,length(y)-1)),cex=c(1.5,rep(.8,length(y)-1)),
               col=c("#134D60",rep("#3A5966",length(y)-1)))
        axis(1,at=tick,cex.axis=.75,mgp=c(2,.4,0),tck=-.015)
        title(xlab=units[j],line=1.65,cex.lab=.84)
        mtext(titles[j],3,line=.8,adj=0,cex=.88,font=2)
        plot.new()
        plot.window(c(0,1),c(.4,12.6))
        fmt <- paste0("%.",digits[j],"f (%.",digits[j],"f, %.",digits[j],"f)")
        text(.5,y,sprintf(fmt,rows$estimate,rows$lower,rows$upper),cex=.74,col="#263438")
        mtext("Estimate (95% CI)",3,line=.8,cex=.72)
      }
    }, finally=dev.off())
  }
})
