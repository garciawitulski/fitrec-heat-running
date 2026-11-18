# Table 1. Each row describes runners, sessions or runner-days.
s <- read_fitrec("panels/sessions.parquet")
s <- s[is.finite(day_max_c)]
d <- read_fitrec("panels/runner_days.parquet")
d <- d[is.finite(day_max_c)]
steady <- read_fitrec("panels/steady_speed.parquet","session_hash")
s[, year := as.integer(format(as.POSIXct(canonical_start_utc,origin="1970-01-01",tz="UTC"),"%Y",tz="UTC"))]
persons <- s[,.(sessions=.N, gender=first(gender),
                span=as.numeric(max(local_date)-min(local_date)),
                hot=any(day_max_c>=30)), by=user_hash]
followup <- d[, .(days=.N),by=user_hash]

# Identify the four countries reported individually in the paper.
shp <- input_file("geography/ne_10m_admin_0_countries.shp")
countries <- sf::st_read(shp,quiet=TRUE)
countries <- countries[countries$ADM0_A3 %in% c("DNK","POL","ESP","USA"),]
# Match session coordinates to country polygons in the longitude–latitude plane.
sf::st_crs(countries) <- NA
points <- sf::st_as_sf(as.data.frame(s[,.(geo_mean_x,geo_mean_y)]),
                      coords=c("geo_mean_x","geo_mean_y"),crs=NA_character_)
hit <- sf::st_intersects(points,countries)
iso <- rep("other",nrow(s))
found <- lengths(hit)>0
iso[found] <- countries$ADM0_A3[vapply(hit[found],function(v) v[1L],integer(1))]

count <- function(row, label, n, denominator) {
  data.table(row=row,variable=label,statistic=c("n","percent"),value=c(n,100*n/denominator))
}
summ <- function(row,label,x,which=c("median","q1","q3")) {
  x <- x[is.finite(x)]
  values <- c(mean=mean(x),sd=sd(x),median=median(x),q1=unname(quantile(x,.25)),
              q3=unname(quantile(x,.75)))
  data.table(row=row,variable=label,statistic=which,value=unname(values[which]))
}
five <- c("mean","sd","median","q1","q3")
np <- nrow(persons)
ns <- nrow(s)
table1_rows <- list(
  count(1,"Male",sum(persons$gender=="male"),np),
  count(2,"Female or gender not reported",sum(persons$gender!="male"),np),
  summ(3,"Sessions per runner",persons$sessions),
  summ(4,"Recording span, days",persons$span),
  count(5,"At least 20 sessions",sum(persons$sessions>=20L),np),
  count(6,"At least one session at or above 30 C",sum(persons$hot),np),
  summ(7,"Duration, min",s$duration_minutes,five),
  summ(8,"Distance, km",s$gps_distance_km),
  summ(9,"Median speed, m/s",s$gps_speed_median_m_s,five),
  summ(10,"Mean heart rate, beats/min",s$hr_mean,five),
  count(11,"Eligible for steady-speed subcohort",sum(s$session_hash %in% steady$session_hash),ns),
  count(12,"Southern Hemisphere",sum(s$geo_mean_y<0),ns),
  summ(13,"Local start time (minutes after midnight)",s$start_minutes_local),
  count(14,"Before 09:00",sum(s$start_before_09),ns),
  count(15,"09:00 to 18:59",sum(!s$start_before_09 & !s$start_after_19),ns),
  count(16,"At or after 19:00",sum(s$start_after_19),ns),
  count(17,"Weekend",sum(s$weekend),ns),
  count(18,"2006-2010",sum(s$year<=2010),ns),
  count(19,"2011",sum(s$year==2011),ns),
  count(20,"2012",sum(s$year==2012),ns),
  count(21,"2013",sum(s$year==2013),ns),
  count(22,"2014",sum(s$year==2014),ns),
  count(23,"2015-2016",sum(s$year>=2015),ns),
  count(24,"Denmark",sum(iso=="DNK"),ns),
  count(25,"Poland",sum(iso=="POL"),ns),
  count(26,"Spain",sum(iso=="ESP"),ns),
  count(27,"United States",sum(iso=="USA"),ns),
  count(28,"Other countries or unassigned",sum(iso=="other"),ns),
  summ(29,"Daily maximum temperature, C",s$day_max_c,five),
  summ(30,"Daily maximum dew point, C",s$day_max_dew_c),
  summ(31,"Session temperature, C",s$temp_c_10km,five),
  summ(32,"Session temperature percentile",100*s$session_day_percentile),
  summ(33,"Heat avoided, C",s$heat_avoided_c,five),
  count(34,"Daily maximum at or below 10 C",sum(s$day_max_c<=10),ns),
  count(35,"Daily maximum at or above 25 C",sum(s$day_max_c>=25),ns),
  count(36,"Daily maximum at or above 30 C",sum(s$day_max_c>=30),ns),
  summ(37,"Days per runner",followup$days),
  count(38,"Days with recorded run",sum(d$run),nrow(d)))
table1 <- rbindlist(table1_rows)
stopifnot(nrow(table1)==101L)
fwrite(table1,file.path(out,"table1.csv"))
