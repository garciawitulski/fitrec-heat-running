# Figure 3. Pool sessions across runners to calculate unadjusted averages.
local({
  s <- read_fitrec("panels/sessions.parquet",c("user_hash","cell_10km",
    "day_start_utc","day_end_utc","day_max_c","temp_c_10km","heat_avoided_c","start_minutes_local"))[is.finite(day_max_c)]
  days <- read_fitrec("panels/runner_days.parquet",
    c("user_hash","day_max_c","run","weekend"))[is.finite(day_max_c)]
  bands <- c("le10","10to20","20to25","25to30","gt30")
  s[,band:=cut(day_max_c,c(-Inf,10,20,25,30,Inf),labels=bands)]
  s[,hour:=pmin(23L,pmax(0L,as.integer(start_minutes_local%/%60)))]
  tab <- s[,.(sessions=.N,runners=uniqueN(user_hash),mean_daily_max_c=mean(day_max_c),
    mean_session_temp_c=mean(temp_c_10km),mean_heat_avoided_c=mean(heat_avoided_c)),by=band]
  stopifnot(all(tab$sessions>=100L & tab$runners>=10L))
  h <- s[,.N,by=.(band,hour)]
  h <- merge(CJ(band=factor(bands,levels=bands),hour=0:23),h,all.x=TRUE)
  h[is.na(N),N:=0L]
  h[,share:=N/sum(N),by=band]
  edges <- seq(-6,36,2)
  s[,bin:=cut(day_max_c,edges,labels=FALSE)]
  bins <- s[!is.na(bin),.(sessions=.N,runners=uniqueN(user_hash),mean_daily_max_c=mean(day_max_c),
    mean_session_temp_c=mean(temp_c_10km),mean_heat_avoided_c=mean(heat_avoided_c)),by=bin][sessions>=100L & runners>=10L]
  bins[,bin_low_c:=edges[bin]]
  days[,bin:=cut(day_max_c,edges,labels=FALSE)]
  participation <- days[!is.na(bin),.(runner_days=.N,runners=uniqueN(user_hash),run_share=mean(run)),
    by=.(bin,weekend)][runner_days>=100L & runners>=10L]
  participation[,bin_low_c:=edges[bin]]

  # The local weather archive stores the grid cell, time-window start, series position and temperature
  # in four one-dimensional NumPy arrays within a ZIP file.
  # Read the original hourly ERA5-Land values with base R.
  cache <- input_file("cohort/hourly_era5_land.npz")
  read_vector <- function(name,type) {
    con <- unz(cache,paste0(name,".npy"),open="rb")
    on.exit(close(con))
    header_start <- readBin(con,"raw",8)
    stopifnot(identical(header_start[1:6],as.raw(c(147,78,85,77,80,89))),as.integer(header_start[7])==1)
    length_header <- readBin(con,integer(),1,size=2,signed=FALSE,endian="little")
    header <- rawToChar(readBin(con,"raw",length_header))
    n <- as.integer(sub(".*'shape': *\\(([0-9]+),.*","\\1",header))
    dtype <- sub(".*'descr': *'([^']+)'.*","\\1",header)
    if (type=="key") {
      stopifnot(grepl("^<U[0-9]+$",dtype))
      width <- as.integer(sub("^<U","",dtype))
      key <- matrix(readBin(con,integer(),n*width,size=4,endian="little"),nrow=width)
      return(apply(key,2,function(chars) intToUtf8(chars[chars>0L])))
    }
    if (type=="offset") {
      stopifnot(dtype=="<i8")
      words <- matrix(readBin(con,integer(),2L*n,size=4,endian="little"),nrow=2)
      stopifnot(all(words[2,]==0L),all(words[1,]>=0L))
      return(as.numeric(words[1,]))
    }
    stopifnot(dtype=="<f4")
    readBin(con,numeric(),n,size=4,endian="little")
  }
  windows <- data.table(cell_10km=read_vector("keys","key"),start=read_vector("firsts","offset"))
  offsets <- read_vector("offsets","offset")
  temp <- read_vector("t","temperature")
  stopifnot(tail(offsets,1)==length(temp),length(offsets)==nrow(windows)+1L)
  windows[, `:=`(pos=head(offsets,-1),last=start+3600*(diff(offsets)-1),archive_order=.I)]
  counts <- s[,.N,by=.(cell_10km,band,day_start_utc,day_end_utc)]
  counts[,id:=.I]
  counts[, first_hour:=ceiling(as.numeric(day_start_utc)/3600)*3600]
  q <- counts[,.(hour=seq_len(min(24L,ceiling((as.numeric(day_end_utc)-first_hour)/3600)))-1L),
    by=.(id,cell_10km,band,first_hour,N)]
  q[,utc:=first_hour+hour*3600]
  q[,qid:=.I]
  matches <- windows[q,on=.(cell_10km,start<=utc,last>=utc),allow.cartesian=TRUE,
    .(qid=i.qid,archive_order,pos=x.pos+(i.utc-x.start)/3600+1)]
  stopifnot(!anyNA(matches$pos))
  setorder(matches,qid,-archive_order)
  matches <- unique(matches,by="qid")
  q[matches,on="qid",temperature:=temp[i.pos]-273.15]
  profile <- q[,.(temperature=weighted.mean(temperature,N)),by=.(band,hour)]

  # Check every displayed cell at the precision of its published aggregate.
  ref <- jsonlite::fromJSON(input_file("references/descriptive/daily_and_hourly_temperature.json"))
  check <- function(block,label,r,paper,digits) data.table(block,label,R=r,paper=paper,
    status=ifelse(abs(r-paper)<=0.5*10^(-digits)+1e-10,"MATCH","DIFFERENCE"))
  checks <- list()
  for (b in bands) {
    z <- ref$bands[[b]]
    for (what in c("sessions","mean_daily_max_c","mean_session_temp_c","mean_heat_avoided_c"))
      checks[[length(checks)+1L]] <- check("temperature_band",paste(b,what),tab[band==b][[what]],z[[what]],if(what=="sessions") 0 else 2)
    checks[[length(checks)+1L]] <- check("start_hour",paste(b,0:23),h[band==b][order(hour),share],z$start_hour_share,4)
    checks[[length(checks)+1L]] <- check("hourly_temperature",paste(b,0:23),profile[band==b][order(hour),temperature],z$mean_hourly_temp_c_by_local_hour,2)
  }
  ref_bins <- as.data.table(ref$session_temperature_vs_daily_maximum)
  setorder(bins,bin_low_c)
  setorder(ref_bins,bin_low_c)
  stopifnot(identical(as.numeric(bins$bin_low_c),as.numeric(ref_bins$bin_low_c)))
  for (what in c("sessions","mean_daily_max_c","mean_session_temp_c","mean_heat_avoided_c"))
    checks[[length(checks)+1L]] <- check("temperature_bin",paste(bins$bin_low_c,what),bins[[what]],ref_bins[[what]],if(what=="sessions") 0 else 2)
  ref_days <- as.data.table(ref$participation_by_bin_and_weekend)
  setorder(participation,bin_low_c,weekend)
  setorder(ref_days,bin_low_c,weekend)
  stopifnot(all(participation$bin_low_c==ref_days$bin_low_c),all(participation$weekend==ref_days$weekend))
  for (what in c("runner_days","run_share")) checks[[length(checks)+1L]] <- check("participation_bin",
    paste(participation$bin_low_c,participation$weekend,what),participation[[what]],ref_days[[what]],if(what=="run_share") 4 else 0)
  fwrite(rbindlist(checks),file.path(out,"figure3_verification.csv"))
  # Plot the estimates from this run in the order used in the paper.
  grDevices::pdf(file.path(out,"figure3_descriptive.pdf"),width=10,height=7)
  tryCatch({
    par(mfrow=c(2,2),mar=c(4,4,1.5,1))
    cols <- c("#173D70","#476F9F","#538173","#CB6C46","#8A3030")
    for (v in c("temperature","share")) {
      frame <- if (v=="temperature") profile else h
      mat <- as.matrix(dcast(frame,hour~band,value.var=v)[,..bands])
      matplot(0:23,mat,type="l",lty=1,col=cols,lwd=2,xlab="Hour of local day",
        ylab=if(v=="temperature") "Hourly temperature (°C)" else "Share of runs starting")
      if(v=="temperature") legend("topleft",c("<=10","10–20","20–25","25–30",">30"),
        col=cols,lty=1,title="Daily maximum (°C)",cex=.66,bty="n")
    }
    plot(bins$mean_daily_max_c,bins$mean_session_temp_c,type="b",pch=16,
      xlab="Daily maximum (°C)",ylab="Mean session temperature (°C)")
    abline(0,1,lty=2,col="grey50")
    plot(range(edges),range(participation$run_share),type="n",xlab="Daily maximum (°C)",ylab="Share of days with a run")
    for (w in 0:1) with(participation[weekend==w],lines(bin_low_c+1,run_share,type="b",col=c("#205782","#A94328")[w+1],pch=16))
    legend("topright",c("Weekdays","Weekends"),col=c("#205782","#A94328"),lty=1,bty="n")
  },finally=grDevices::dev.off())
})
