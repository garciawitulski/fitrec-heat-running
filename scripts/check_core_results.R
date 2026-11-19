# Compare the paired models with Table S13 and the decomposition with Figure 3.
# For paired confidence limits, use the range allowed by the reported runner count.
local({
  figure_data <- file.path(repo,"data","published","figures")
  paired <- fread(file.path(public_out,"models","paired_heart_rate.csv"))
  paired_rows <- data.table(row=c(6L,7L,8L,37L,38L,39L),
    quantity=c("steady_fixed_speed","steady_speed_free","whole_session_mean",
               "steady_speed_free_minus_whole_session","steady_fixed_speed_minus_whole_session",
               "fixed_speed_minus_speed_free"))
  cells <- fread(file.path(repo,"data","published","table_cells.csv"))
  reported <- merge(paired_rows,cells[table_id=="S13" & column %in% 2:3,.(row,column,text)],by="row")
  reported[,target:=ifelse(column==2L,"P90","P95")]
  numbers <- strsplit(gsub("[(),]"," ",chartr("\u2212","-",reported$text)),"[[:space:]]+")
  values <- t(vapply(numbers,as.numeric,numeric(3)))
  reported[,c("printed_estimate","printed_lower","printed_upper"):=as.data.table(values)]
  z <- merge(paired,reported,by=c("quantity","target"),all=TRUE)
  stopifnot(nrow(z)==12L,!anyDuplicated(z,by=c("quantity","target")),!anyNA(z))
  z[,matches:=round(estimate,2)==printed_estimate &
    round(lower_conservative,2)<=printed_lower & printed_lower<=round(lower_at_max_df,2) &
    round(upper_at_max_df,2)<=printed_upper & printed_upper<=round(upper_conservative,2)]
  decomp <- fread(file.path(public_out,"models","heart_rate_decomposition.csv"))
  levels <- fread(file.path(figure_data,"figure3_stages.csv"))
  setnames(levels,c("series","estimate"),c("quantity","reference"))
  stages <- merge(decomp[target=="P95"],levels,by="quantity",all.y=TRUE)
  stopifnot(nrow(stages)==5L,!anyNA(stages))
  stages[,matches:=abs(estimate-reference)<1e-7 & abs(lower-ci_low)<1e-6 & abs(upper-ci_high)<1e-6]
  # The article reports the paired difference to two decimal places and
  # the reduction associated with timing and speed to a whole percentage.
  headline <- rbindlist(list(
    z[target=="P95" & quantity=="steady_speed_free_minus_whole_session",
      .(quantity,estimate,reference=0.81,matches=round(estimate,2)==0.81 &
        round(lower_conservative,2)==0.29 & round(lower_at_max_df,2)==0.29 &
        round(upper_conservative,2)==1.33 & round(upper_at_max_df,2)==1.33)],
    decomp[target=="P95" & quantity=="percentage_reduction_with_timing_and_speed",
      .(quantity,estimate,reference=30,matches=round(estimate)==30)]))
  stopifnot(nrow(headline)==2L)
  dest <- file.path(public_out,"validation")
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  fwrite(z,file.path(dest,"heart_rate_pairs.csv"))
  fwrite(stages,file.path(dest,"heart_rate_stages.csv"))
  fwrite(headline,file.path(dest,"heart_rate_headlines.csv"))
  if (!all(c(z$matches,stages$matches,headline$matches))) stop("The heart-rate comparisons did not reproduce the manuscript.")
})
