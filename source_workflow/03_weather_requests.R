# Construct the first set of hourly weather requests from the retained runs.
# Run at the repository root after 02_source_counts.py --prepare-weather-input.
library(data.table)

build_weather_plan <- function(runs, mask_path) {
  stopifnot(nrow(runs) == 115036L, uniqueN(runs$runner_index) == 862L,
            all(is.finite(runs$median_latitude)), all(is.finite(runs$median_longitude)),
            all(runs$end_utc_seconds > runs$start_utc_seconds))
  mask <- ncdf4::nc_open(mask_path)
  on.exit(ncdf4::nc_close(mask))
  land <- ncdf4::ncvar_get(mask, "lsm")
  mask_lat <- round(mask$dim$latitude$vals * 10)
  mask_lon <- (round(mask$dim$longitude$vals * 10) + 1800) %% 3600 - 1800

  # Dates include the hour before the start and the hour after the end when
  # interpolation requires them. Each run lasts at most one day.
  runs <- copy(runs)
  runs[, first_day := floor(start_utc_seconds / 86400)]
  runs[, last_day := floor(ceiling(end_utc_seconds / 3600) / 24)]
  runs[, lat_land := floor(median_latitude * 10 + 0.5)]
  runs[, lon_land := (floor(median_longitude * 10 + 0.5) + 1800) %% 3600 - 1800]
  runs[, lat_era5 := round(median_latitude / 0.25) * 0.25]
  runs[, lon_era5 := (round(median_longitude / 0.25) * 0.25 + 180) %% 360 - 180]
  blocks <- unique(runs[, .(lat_land, lon_land, lat_era5, lon_era5, first_day, last_day)])
  days <- blocks[, .(day = seq.int(first_day, last_day)),
                 by = .(lat_land, lon_land, lat_era5, lon_era5, first_day, last_day)]

  cells <- unique(days[, .(lat_land, lon_land)])
  cells[, fraction := land[cbind(match(lon_land, mask_lon), match(lat_land, mask_lat))]]
  centres <- cells[!is.na(fraction) & fraction > 0.5,
                   .(lat_land, lon_land, latitude = lat_land / 10, longitude = lon_land / 10)]
  coast <- cells[is.na(fraction) | fraction <= 0.5,
                 .(dlat = rep(-1:1, each = 3), dlon = rep(-1:1, 3)), by = .(lat_land, lon_land)]
  coast <- coast[(dlat != 0 | dlon != 0) & abs(lat_land + dlat) <= 900]
  coast[, latitude := (lat_land + dlat) / 10]
  coast[, longitude := ((lon_land + dlon + 1800) %% 3600 - 1800) / 10]
  coast[, fraction := land[cbind(match(round(longitude * 10), mask_lon),
                                  match(round(latitude * 10), mask_lat))]]
  centres <- rbind(centres, coast[!is.na(fraction) & fraction > 0.5,
                                   .(lat_land, lon_land, latitude, longitude)])

  land_days <- merge(centres, unique(days[, .(lat_land, lon_land, day)]),
                     by = c("lat_land", "lon_land"), allow.cartesian = TRUE)
  needed <- unique(rbind(
    land_days[, .(dataset = "reanalysis-era5-land-timeseries", latitude, longitude, day)],
    days[, .(dataset = "reanalysis-era5-single-levels-timeseries", latitude = lat_era5,
             longitude = lon_era5, day)]
  ))
  setorder(needed, dataset, latitude, longitude, day)
  needed[, period := cumsum(c(TRUE, diff(day) > 365)), by = .(dataset, latitude, longitude)]
  plan <- needed[, .(start_date = as.IDate(min(day), origin = "1970-01-01"),
                     end_date = as.IDate(max(day), origin = "1970-01-01")),
                 by = .(dataset, latitude, longitude, period)]
  plan[, period := NULL]
  totals <- plan[, .(grid_cells = uniqueN(.SD, by = c("latitude", "longitude")),
                     request_windows = .N), by = dataset]
  if (!identical(totals$grid_cells, c(6003L, 3321L)) ||
      !identical(totals$request_windows, c(6903L, 3896L))) {
    stop("Initial weather request counts differ from the published workflow: ",
         paste(totals$grid_cells, totals$request_windows, collapse = " / "))
  }
  list(requests = plan, totals = totals)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) <= 1L)
  local_dir <- file.path("data", "processed", "local")
  runs_path <- file.path(local_dir, "preweather_sessions.csv")
  if (!file.exists(runs_path)) stop("Run 02_source_counts.py --prepare-weather-input first.")
  mask_name <- "lsm_1279l4_0.1x0.1.grb_v4_unpack.nc"
  mask_path <- if (length(args)) args[1L] else file.path("data", "raw", mask_name)
  if (!file.exists(mask_path)) {
    dir.create(dirname(mask_path), recursive = TRUE, showWarnings = FALSE)
    mask_url <- paste0("https://confluence.ecmwf.int/download/attachments/140385202/",
                       mask_name, "?version=1&modificationDate=1591983422208&api=v2")
    mask_temp <- tempfile(tmpdir = dirname(mask_path))
    utils::download.file(mask_url, mask_temp, mode = "wb", method = "libcurl")
  } else {
    mask_temp <- mask_path
  }
  mask_hash <- "9df65151d1608990851e2aa8d0e1e4ca063b5181b0568493c42676d2572b202c"
  if (unname(tools::sha256sum(mask_temp)) != mask_hash) {
    if (mask_temp != mask_path) unlink(mask_temp)
    stop("The ECMWF invariant land-sea mask differs from the study version.")
  }
  if (mask_temp != mask_path && !file.rename(mask_temp, mask_path)) stop("Could not save the verified mask.")
  result <- build_weather_plan(fread(runs_path), mask_path)
  # This plan contains dates and coordinates. Keep it in the ignored input directory.
  fwrite(result$requests, file.path(local_dir, "initial_weather_requests.csv"))
  dir.create(file.path("outputs", "validation"), recursive = TRUE, showWarnings = FALSE)
  fwrite(result$totals, file.path("outputs", "validation", "initial_weather_requests.csv"))
  print(result$totals)
}
