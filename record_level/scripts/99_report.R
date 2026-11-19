# Combine the comparisons from this run. Each count represents the number of
# values checked within an analysis.
local({
  checked <- list()
  files <- list.files(out,pattern="^check_.*\\.csv$",full.names=TRUE)
  for (file in files) {
    x <- fread(file)
    checked[[basename(file)]] <- x[,.(n=.N),by=status][,block:=tools::file_path_sans_ext(basename(file))]
  }
  x <- fread(file.path(out,"primary_comparison.csv"))
  checked$tables <- x[,.(n=.N),by=.(block=paste0("printed_",table),status)]
  tables <- x[table %in% c("T1","T2"),.(table,outcome=as.character(outcome),target,statistic,printed,value,status)]
  labels <- unique(fread(file.path(out,"table1.csv"))[,.(outcome=as.character(row),variable)])
  tables[,variable:=outcome]
  tables[table=="T1",variable:=labels$variable[match(outcome,labels$outcome)]]
  setcolorder(tables,c("table","variable","target","statistic","printed","value","status","outcome"))
  fwrite(tables,file.path(out,"main_table_comparison.csv"))
  x <- fread(file.path(out,"heart_verification.csv"))
  x[,status:=toupper(status)]
  x[status=="DIFFERENT" & grepl("heart_rate_decomposition_curves.csv$",source_file) &
    grepl("^-0\\.5 / ",quantity) & abs(difference)<=2.1e-6 & abs(reference)>1e5,
    status:="NUMERICAL_DIFFERENCE"]
  checked$heart <- x[,.(n=.N),by=.(block=paste0("heart_",basename(source_file)),status)]
  x <- fread(file.path(out,"weather_verification.csv"))
  checked$weather <- x[,.(n=.N),by=.(block=paste0("weather_",source),status)]
  x <- fread(file.path(out,"figure2_verification.csv"))
  checked$descriptive <- x[,.(n=.N),by=.(block=paste0("fig2_",block),status)]
  x <- fread(file.path(out,"cohort_flow.csv"))
  checked$flow <- x[,.(n=.N),by=status][,block:="S1_cohort_stages"]
  if (file.exists(file.path(out,"figure_S1_aggregate_audit.csv"))) {
    x <- fread(file.path(out,"figure_S1_aggregate_audit.csv"))
    checked$map <- x[,.(block=paste0("supp_figure1_",measure),status="MATCH",n=matched)]
    stopifnot(all(x$mismatched==0))
  }
  counts <- rbindlist(checked,use.names=TRUE)
  counts <- counts[,.(n=sum(n)),by=.(block,status)]
  setorder(counts,block,status)
  fwrite(counts,file.path(out,"verification_summary.csv"))
  print(counts[,.(n=sum(n)),by=status])
  problem <- counts[!status %in% c("MATCH","EXPLAINED","NOT_APPLICABLE","DIFFERENT_RNG","NUMERICAL_DIFFERENCE")]
  if(nrow(problem)) stop("Some comparisons are unresolved. See verification_summary.csv.")
})
