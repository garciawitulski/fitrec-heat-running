# Figure 6 and the main activity-volume estimates. After removing runner
# intercepts and seasonal/calendar controls, each model solves A*beta = b
# using one pooled cross-product matrix A and one outcome cross-product b.
local({
  library(data.table)
  source(file.path(repo,"scripts","model_tools.R"),encoding="UTF-8",local=TRUE)
  input <- file.path(repo,"data","processed","activity_model")
  matrices <- fread(file.path(input,"pooled_matrices.csv"))
  vectors <- fread(file.path(input,"pooled_vectors.csv"))
  a <- fread(file.path(repo,"data","processed","model_statistics","temperature_percentiles.csv"))
  knots <- a[percentile %in% c(10,50,90),temperature]
  basis <- rcs_basis(a[percentile %in% c(50,90,95),temperature],knots)
  targets <- c("P90","P95")
  rows <- list()
  for (id in unique(vectors$model)) {
    m <- matrices[model==id]
    v <- vectors[model==id]
    terms <- v$term
    gram <- meat <- matrix(0,length(terms),length(terms),dimnames=list(terms,terms))
    ij <- cbind(match(m$term1,terms),match(m$term2,terms))
    gram[ij] <- gram[ij[,2:1,drop=FALSE]] <- m$xx
    meat[ij] <- meat[ij[,2:1,drop=FALSE]] <- m$meat
    beta <- solve(gram,v$xy)
    bread <- solve(gram)
    cov <- bread %*% meat %*% bread
    contrasts <- list()
    if (startsWith(id,"weekly_")) {
      for (term in terms) {
        h <- setNames(numeric(length(terms)),terms)
        h[term] <- 1
        contrasts[[paste("band",term,sep="/")]] <- h
      }
    } else for (i in 1:2) {
      delta <- basis[i+1,]-basis[1,]
      if (id!="running_lags") {
        contrasts[[paste(targets[i],"day",sep="/")]] <- delta
      } else {
        # Each day is counted once if the runner recorded at least one run with
        # heart rate data. The sum from day 0 to j gives the change in runner-days
        # with a run per 100 windows of j + 1 days (eight days when j = 7).
        day <- -3:7
        pieces <- list()
        for (j in day) pieces[[paste0("day_",j)]] <- j
        for (j in 0:7) pieces[[paste0("cumulative_0_",j)]] <- 0:j
        pieces[["sum_1_7"]] <- 1:7
        pieces[["sum_minus3_minus1"]] <- -3:-1
        for (name in names(pieces)) {
          h <- setNames(numeric(length(terms)),terms)
          for (k in pieces[[name]]) {
            key <- if (k<0) paste0("lead",abs(k)) else paste0("lag",k)
            h[paste(key,names(delta),sep="_")] <- delta
          }
          contrasts[[paste(targets[i],name,sep="/")]] <- h
        }
      }
    }
    H <- do.call(rbind,contrasts)
    multiplier <- if (startsWith(id,"weekly_")) 1 else 100
    est <- multiplier*as.numeric(H %*% beta)
    se <- multiplier*sqrt(pmax(0,rowSums((H %*% cov)*H)))
    bounds <- if (id %in% c("any_sport","outdoor_cycling")) c(729,729) else c(720,729)
    wide <- qt(.975,bounds[1]-1); narrow <- qt(.975,bounds[2]-1)
    labels <- tstrsplit(rownames(H),"/",fixed=TRUE)
    rows[[id]] <- data.table(model=id,target=labels[[1]],contrast=labels[[2]],estimate=est,se=se,
      lower_min=est-wide*se,lower_max=est-narrow*se,upper_min=est+narrow*se,upper_max=est+wide*se,
      p_min=2*pt(-abs(est/se),bounds[2]-1),p_max=2*pt(-abs(est/se),bounds[1]-1))
  }
  all_models <- rbindlist(rows)
  tables <- file.path(public_out,"tables","supplementary")
  figures <- file.path(public_out,"figures","manuscript")
  dir.create(tables,recursive=TRUE,showWarnings=FALSE)
  dir.create(figures,recursive=TRUE,showWarnings=FALSE)
  fwrite(all_models,file.path(tables,"activity_models.csv"))
  # Test eight headline contrasts against the paper and its supplement.
  published <- fread(file.path(input,"reported_activity.csv"))
  setnames(published,c("estimate","lower","upper"),c("paper_estimate","paper_lower","paper_upper"))
  check <- merge(published,all_models,by=c("model","target","contrast"),all.x=TRUE)
  check[, tolerance:=.5*10^(-digits)+1e-12]
  check[, matches:=abs(estimate-paper_estimate)<=tolerance &
    lower_min>=paper_lower-tolerance & lower_max<=paper_lower+tolerance &
    upper_min>=paper_upper-tolerance & upper_max<=paper_upper+tolerance]
  dir.create(file.path(public_out,"validation"),recursive=TRUE,showWarnings=FALSE)
  fwrite(check,file.path(public_out,"validation","activity_printed.csv"))
  stopifnot(nrow(check)==8L,all(check$matches))
  sports <- all_models[model %in% c("any_sport","outdoor_cycling") & target=="P95",
    .(outcome=model,estimate,lower=lower_min,upper=upper_max,reference="Refitted from pooled model statistics")]
  reported <- fread(file.path(input,"published_sports.csv"))
  sports <- rbind(sports,reported)
  sport_order <- c("any_sport","outdoor_running","treadmill_running","outdoor_cycling","cycling_recreation",
    "cycling_transport","indoor_cycling","walking_hiking","indoor_training","swimming","indoor_all")
  row_index <- match(sport_order,sports$outcome)
  stopifnot(!anyNA(row_index))
  sports <- sports[row_index]
  fwrite(sports,file.path(tables,"figure6_activity_types.csv"))

  # Figure 6a uses the outer confidence limits allowed by the reported runner counts.
  profile <- all_models[model=="running_lags" & target=="P95" & startsWith(contrast,"day_")]
  cumulative <- all_models[model=="running_lags" & target=="P95" & startsWith(contrast,"cumulative_")]
  pdf(file.path(figures,"figure_6_participation.pdf"),width=12,height=5.5,useDingbats=FALSE)
  tryCatch({
    par(mfrow=c(1,2),mar=c(4.5,4.5,2.5,1))
    plot(-3:7,profile$estimate,ylim=range(profile$lower_min,profile$upper_max,cumulative$lower_min,cumulative$upper_max),
      type="n",xlab="Days relative to the hot day",ylab="Change in runner-days per 100 windows",
      main="a. Days with a recorded run")
    mtext("A day counts if at least one recorded run has heart rate data.",side=3,line=.2,cex=.7)
    abline(h=0,lty=3)
    polygon(c(0:7,7:0),c(cumulative$lower_min,rev(cumulative$upper_max)),border=NA,col="#C6DCE7")
    lines(0:7,cumulative$estimate,lwd=2,col="#2B5D79")
    segments(-3:7,profile$lower_min,-3:7,profile$upper_max)
    points(-3:7,profile$estimate,pch=21,bg="white")
    legend("topright",c("Daily (percentage points)","Cumulative from day 0"),lty=c(NA,1),pch=c(1,NA),bty="n",cex=.8)
    par(mar=c(4.5,8.5,2.5,1))
    pos <- rev(seq_len(nrow(sports)))
    plot(range(sports$lower,sports$upper),range(pos)+c(-.5,.5),type="n",yaxt="n",ylab="",
         xlab="Change in probability (percentage points)",main="b. Workouts at P95 versus P50")
    mtext("Filled points show refits. Open points show published estimates.",side=3,line=.2,cex=.7)
    abline(v=0,lty=3)
    axis(2,at=pos,labels=gsub("_"," ",sports$outcome),las=1,tick=FALSE,cex.axis=.8)
    segments(sports$lower,pos,sports$upper,pos)
    points(sports$estimate,pos,pch=21,bg=ifelse(sports$outcome %in% c("any_sport","outdoor_cycling"),"#2B5D79","white"))
  },finally=dev.off())
})
