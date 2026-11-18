# Heat, heart rate and place (Table S9). The usual cell is identified before
# excluding sessions with missing daily weather. Ties go to the smaller cell ID.
cell_counts <- hr_all[!is.na(cell_10km),.N,by=.(user_hash,cell_10km)]
setorder(cell_counts,user_hash,-N,cell_10km)
usual_cell <- cell_counts[!duplicated(user_hash)]
is_usual <- hr_main$cell_10km==usual_cell$cell_10km[match(hr_main$user_hash,usual_cell$user_hash)]
is_usual[is.na(is_usual)] <- FALSE
edge <- read_fitrec("panels/sessions_edge.parquet",names(hr_all))[is.finite(day_max_c)]
edge_x <- hr_design(edge)
edge_xhour <- hr_design(edge,hour="start_minutes_local")
edge_hour <- fit_fe(edge$heat_avoided_c,edge_xhour,edge$user_hash,edge$user_hash)
ids_spatial <- c("b1_heat_avoided","b2_speed","b3_log_duration","h1_hr_mean","s_start_minutes")
fields_spatial <- c("heat_avoided_c","gps_speed_median_m_s","duration_minutes","hr_mean","start_minutes_local")
spatial_fits <- edge_fits <- month_estimates <- list()
for (i in seq_along(ids_spatial)) {
    id <- ids_spatial[i]
    y <- hr_main[[fields_spatial[i]]]
    yedge <- edge[[fields_spatial[i]]]
    if (id=="b3_log_duration") {
        y <- log(y)
        yedge <- log(yedge)
    }
    cell_id <- paste(hr_main$user_hash,hr_main$cell_10km)
    a <- fit_fe(y,hr_xmain,cell_id,hr_main$user_hash)
    b <- fit_fe(y[is_usual],hr_xmain[is_usual,,drop=FALSE],hr_main$user_hash[is_usual],hr_main$user_hash[is_usual])
    spatial_fits[[id]] <- rbind(hr_contrasts(a,paste0(id,"_runner_cell")),hr_contrasts(b,paste0(id,"_usual_cell")))
    edge_fits[[id]] <- fit_fe(yedge,edge_x,edge$user_hash,edge$user_hash)
    m <- if (id=="b1_heat_avoided") models_month$gap else if (id=="b2_speed") models_month$speed else if (id=="h1_hr_mean") models_month$net else fit_fe(y,month_x,month_id,hr_main$user_hash)
    month_estimates[[id]] <- hr_contrasts(m,paste0(id,"_runner_month"))
}
spatial_results <- rbindlist(spatial_fits)
fwrite(spatial_results,file.path(out,"heart_place_checks.csv"))
fwrite(rbindlist(month_estimates),file.path(out,"heart_month_checks.csv"))
fwrite(rbindlist(lapply(ids_spatial,function(id) hr_contrasts(edge_fits[[id]],id))),file.path(out,"heart_edge_checks.csv"))
models_edge <- list(gap=edge_fits$b1_heat_avoided,hour=edge_hour,speed=edge_fits$b2_speed,
                    net=edge_fits$h1_hr_mean,phys=hr_phfit)
edge_basis <- as.numeric(rcs_basis(hr_targets[["P50"]],hr_kday))-colMeans(edge_x[,c("T1","T2")])
edge_parts <- rbindlist(lapply(c("P90","P95"),function(label)
    hr_decomposition(models_edge,hr_targets[[label]],mean(edge$heat_avoided_c),edge_basis)[,target_label:=label]))
fwrite(edge_parts,file.path(out,"heart_edge_decomposition.csv"))

# Complete registered temperature grid; the three percentile anchors join the
# 0.5-degree sequence. Fractions at the median have zero denominators.
env <- hr_config$plot_envelope
grid <- sort(unique(c(seq(ceiling(env[1]*2)/2,floor(env[2]*2)/2,by=.5),hr_targets[c("P90","P95")])))
curve <- rbindlist(lapply(grid,function(value)
    hr_decomposition(hr_models,value)[,exposure_c:=value]))
fwrite(curve,file.path(out,"heart_decomposition_curves.csv"))

# Figure 2. Matched levels and paired differences; separate marks for P90/P95.
blue <- "#21688A"
amber <- "#BF5735"
fig_paired <- function(keys,labels,title) {
    a <- hr_paired_rows[quantity %in% keys]
    plot(NA,xlim=range(c(0,a$lower,a$upper)) + c(-.15,.15),ylim=c(.45,length(keys)+.55),
         yaxt="n",ylab="",xlab="Change in heart rate (beats/min)",bty="n",main=title,adj=0)
    axis(2,at=rev(seq_along(keys)),labels=labels,las=1,tick=FALSE)
    abline(v=0,col="grey65",lty=2)
    for (lab in c("P90","P95")) {
        d <- a[target_label==lab]
        y <- length(keys)+1L-match(d$quantity,keys)+if (lab=="P90") .10 else -.10
        segments(d$lower,y,d$upper,y,col=blue,lwd=1.5)
        points(d$estimate,y,pch=if (lab=="P90") 21 else 19,col=blue,bg="white",cex=1.2)
    }
    legend("bottomright",legend=c("P90","P95"),pch=c(21,19),col=blue,pt.bg="white",bty="n",horiz=TRUE,cex=.85)
}
pdf(file.path(out,"figure_2_heart_rate.pdf"),width=9.2,height=6,useDingbats=FALSE)
par(mfrow=c(2,1),mar=c(3.3,13,2.1,1),mgp=c(1.7,.5,0),family="sans",cex.axis=.89)
fig_paired(c("steady_fixed_speed","steady_speed_free","whole_session_mean"),
           c("Steady, fixed speed","Steady, speed free","Whole session"),"A   Heart rate in the same runs")
fig_paired(c("steady_speed_free_minus_whole_session","steady_fixed_speed_minus_whole_session"),
           c("Speed free minus whole session","Fixed speed minus whole session"),"B   Paired differences")
invisible(dev.off())

# Figure 5. Fixed-hour benchmark, contrasts along the temperature range and checks.
levels_fig <- c("rise_without_behavioural_change","rise_at_unchanged_clock_hour","rise_with_observed_timing",
                "predicted_rise_timing_and_speed","net_heart_rate_change")
labels_fig <- c("Initial benchmark","Fixed hour","Observed timing","Timing and speed","Measured heart rate")
parts_fig <- hr_dec[target_label=="P95"][match(levels_fig,quantity)]
cold <- "#0B4E63"
colors <- c("#70A4B5",cold,"#258D88","#9D7154",amber)
pdf(file.path(out,"figure_5_decomposition.pdf"),width=11.7,height=8.1,useDingbats=FALSE)
layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(1,1.3))
par(mar=c(5.4,4.4,2.1,.7),mgp=c(2.4,.55,0),family="sans",cex.axis=.86)
limit <- range(c(0,parts_fig$lower,parts_fig$upper))
plot(NA,xlim=c(.5,5.55),ylim=limit+c(-.20,.15),xaxt="n",xlab="",ylab="Heart-rate change (beats/min)",
     main="A   Heart-rate accounting at P95",adj=0,bty="n")
abline(h=0,col="grey70")
for (j in seq_len(5L)) {
    rect(j-.23,0,j+.23,parts_fig$estimate[j],col=colors[j],border=NA)
    segments(j,parts_fig$lower[j],j,parts_fig$upper[j])
    if (j<5) rect(j+.31,parts_fig$estimate[j+1],j+.68,parts_fig$estimate[j],col="grey78",border=NA)
}
text(1:5,par("usr")[3],labels_fig,srt=28,adj=1,xpd=NA,cex=.76)
par(mar=c(4.4,4.1,2.1,1))
hot_curve <- curve[exposure_c>hr_targets[["P50"]] & quantity %in% levels_fig]
plot(NA,xlim=range(c(hr_targets[["P50"]],hot_curve$exposure_c)),ylim=range(hot_curve$lower,hot_curve$upper),
     xlab="Daily maximum temperature (degrees C)",ylab="Heart-rate change (beats/min)",
     main="B   Hotter than the median day",adj=0,bty="n")
for (j in seq_along(levels_fig)) {
    z <- hot_curve[quantity==levels_fig[j]][order(exposure_c)]
    x <- c(hr_targets[["P50"]],z$exposure_c)
    if (j %in% c(2,5)) polygon(c(x,rev(x)),c(0,z$lower,rev(c(0,z$upper))),border=NA,
                                 col=adjustcolor(colors[j],alpha.f=.16))
    lines(x,c(0,z$estimate),col=colors[j],lwd=1.6,lty=j)
}
legend("topleft",labels_fig,col=colors,lty=1:5,lwd=1.6,bty="n",cex=.68)
par(mar=c(3.8,14.2,2.2,1))
q <- "rise_at_unchanged_clock_hour"
alt_hour <- fread(file.path(out,"heart_hour_decomposition.csv"))[quantity==q & target_label=="P95"]
alt_month <- fread(file.path(out,"heart_month_decomposition.csv"))[quantity==q & target_label=="P95"]
net_month <- hr_contrasts(models_month$net,"month")[target_label=="P95"]
locations <- spatial_results[target_label=="P95" & model_id %in% c("h1_hr_mean_runner_cell","h1_hr_mean_usual_cell")]
qq <- c("steady_fixed_speed","steady_speed_free","whole_session_mean","steady_speed_free_minus_whole_session","steady_fixed_speed_minus_whole_session")
series <- c("estimate","lower","upper")
z <- rbindlist(lapply(list(parts_fig[2],alt_hour,alt_month,parts_fig[5],net_month,locations,
                           hr_paired_rows[target_label=="P95"][match(qq,quantity)]),function(d) d[,..series]))
lab <- c("Benchmark","Benchmark: hour and weekday","Benchmark: runner-month",
         "Whole session","Whole session: runner-month","Whole session: runner-cell","Whole session: usual cell",
         "Steady, fixed speed","Steady, speed free","Whole session, same runs",
         "Speed free minus whole","Fixed speed minus whole")
ypos <- rev(seq_along(lab))
plot(NA,xlim=range(z$lower,z$upper)+c(-.15,.15),ylim=c(.4,length(lab)+.5),yaxt="n",
     xlab="P95 versus P50, beats/min",ylab="",main="C   Alternative measurements and specifications",adj=0,bty="n")
abline(v=0,lty=2,col="grey65")
axis(2,at=ypos,labels=lab,las=1,tick=FALSE,cex.axis=.82)
segments(z$lower,ypos,z$upper,ypos,col=cold,lwd=1.5)
points(z$estimate,ypos,pch=19,col=cold)
invisible(dev.off())
