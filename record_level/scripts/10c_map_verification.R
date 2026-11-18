# Figure S1. Compare both map statistics within the approved rectangles in memory.
# Export the total number of matches for each statistic.
local({
  sessions <- read_fitrec("panels/sessions.parquet",
                          c("geo_mean_x","geo_mean_y","day_max_c","user_hash"))[is.finite(day_max_c)]
  regions <- fread(input_file("geography/aggregate_regions.csv"),
                   select=c("squares","sessions","percent_hot_days"))
  map_summary <- jsonlite::fromJSON(input_file("references/descriptive/map_totals.json"))
  stopifnot(nrow(sessions)==map_summary$sessions,nrow(regions)==map_summary$released_regions,
            sum(is.finite(regions$percent_hot_days))==map_summary$regions_with_hot_share)
  x <- pmin(179.999999,pmax(-180,sessions$geo_mean_x))
  y <- pmin(89.999999,pmax(-90,sessions$geo_mean_y))
  count_ok <- logical(nrow(regions))
  percent_ok <- rep(NA,nrow(regions))
  assigned <- logical(nrow(sessions))
  for(i in seq_len(nrow(regions))) {
    squares <- as.matrix(jsonlite::fromJSON(regions$squares[i]))
    inside <- rep(FALSE,nrow(sessions))
    for(j in seq_len(nrow(squares))) {
      size <- squares[j,3]
      inside <- inside | (x>=squares[j,1] & x<squares[j,1]+size & y>=squares[j,2] & y<squares[j,2]+size)
    }
    stopifnot(!any(assigned & inside),sum(inside)>=100L,uniqueN(sessions$user_hash[inside])>=10L)
    assigned <- assigned | inside
    count_ok[i] <- sum(inside)==regions$sessions[i]
    if(is.finite(regions$percent_hot_days[i])) {
      hot <- sessions$day_max_c[inside]>=25
      ids <- sessions$user_hash[inside]
      stopifnot(all(vapply(split(ids,hot),uniqueN,integer(1))>=10L))
      percent_ok[i] <- round(100*mean(hot),1)==regions$percent_hot_days[i]
    }
  }
  stopifnot(round(100*mean(assigned),1)==map_summary$percent_sessions_in_released_regions)
  report <- data.table(measure=c("session counts","released hot-day percentages"),
                       matched=c(sum(count_ok),sum(percent_ok,na.rm=TRUE)),
                       mismatched=c(sum(!count_ok),sum(!percent_ok,na.rm=TRUE)))
  fwrite(report,file.path(out,"figure_S1_aggregate_audit.csv"))
  print(report)
  stopifnot(all(report$mismatched==0))
})
