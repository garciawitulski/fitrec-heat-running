# Compare the daylight, clock-change and weather models with the aggregate
# reference values for these analyses.
local({
  ref_dir <- "references/additional"
  checks <- list()
  compare <- function(name, calculated, reference, by, fields) {
    stopifnot(!anyDuplicated(calculated,by=by),!anyDuplicated(reference,by=by))
    joined <- merge(calculated,reference,by=by,all=TRUE,suffixes=c("_r","_paper"))
    stopifnot(nrow(joined)==nrow(reference),nrow(joined)==nrow(calculated))
    rbindlist(lapply(fields,function(field) {
      computed <- joined[[paste0(field,"_r")]]
      printed <- joined[[paste0(field,"_paper")]]
      # Some aggregates were released to three or four significant figures.
      mantissa <- sub("[eE].*","",printed)
      exponent <- ifelse(grepl("[eE]",printed),as.numeric(sub(".*[eE]","",printed)),0)
      digits <- ifelse(grepl(".",mantissa,fixed=TRUE),nchar(sub(".*\\.","",mantissa)),0)
      precision <- 10^(exponent-digits)
      tolerance <- pmax(1e-8,precision/2)
      paper <- as.numeric(printed)
      cbind(joined[,..by], data.table(source=name,field=field,R=computed,paper=paper,
        tolerance=tolerance,status=ifelse(is.finite(computed) & is.finite(paper) &
          abs(computed-paper)<=tolerance,"MATCH","DIFFERENCE")))
    }))
  }
  fields <- c("estimate","lower","upper","p")
  for (record in c("contrasts","coefficients")) {
    calc <- fread(file.path(out,paste0("weather_",record,".csv")))
    refs <- list()
    for (prefix in c("daylight_sessions","daylight_participation","future_weather","weather_sessions","weather_participation","weather_steady_speed")) {
      path <- input_file(ref_dir,paste0(prefix,"_",record,".csv"))
      z <- fread(path,colClasses="character")
      z[, check:=if (startsWith(prefix,"daylight_")) "daylight" else if (startsWith(prefix,"future_")) "leads" else "weather"]
      if (!"sample" %in% names(z)) z[,sample:="all"]
      if (record=="coefficients") {
        z[,base_model:=sub("_(daylight|leads|a3_weather)$","",model_id)]
        z[base_model=="h2",base_model:="h2_pooled_session_temp"]
        z[,term_no:=fcase(grepl("^sunrise|^lead1|^a3_precip|^a3_s_wind",term),1L,
          grepl("^sunset|^lead2|^a3_wind|^a3_s_rad",term),2L,
          grepl("^lead3|^a3_rad",term),3L)]
        stopifnot(!anyNA(z$term_no))
      }
      setnames(z,c("ci_low","ci_high","p_value"),c("lower","upper","p"))
      refs[[prefix]] <- z
    }
    reference <- rbindlist(refs,fill=TRUE)
    by <- c("check","sample","base_model",if (record=="contrasts") c("version","target_label") else "term_no")
    checks[[record]] <- compare(record,calc,reference[,c(by,fields),with=FALSE],by,fields)
  }
  z <- fread(input_file(ref_dir,"future_weather_joint_tests.csv"),colClasses="character")
  z[,base_model:=sub("_leads$","",model_id)]
  checks[["joint"]] <- compare("lead_joint",fread(file.path(out,"lead_joint.csv")),z,
    c("sample","base_model"),c("F_statistic","p_value"))
  for (name in c("sunrise_starts","clock_changes")) {
    file <- if (name=="sunrise_starts") "sunrise_relative_starts" else "clock_change"
    z <- fread(input_file(ref_dir,paste0(file,".csv")),colClasses="character")
    setnames(z,c("ci_low","ci_high","p_value"),c("lower","upper","p"))
    by <- c("outcome",if (name=="sunrise_starts") "target_label" else "quantity")
    checks[[name]] <- compare(name,fread(file.path(out,paste0(name,".csv"))),z,by,fields)
  }
  fwrite(rbindlist(checks,fill=TRUE),file.path(out,"weather_verification.csv"))
})
