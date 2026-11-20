# Rebuild the FitRec cohort from the original records

These scripts reconstruct the first row of Table S1 from the two original FitRec archives, before linking the sessions to weather. They also prepare the initial requests for hourly weather from Copernicus. Run the commands below from the top-level directory of the repository.

## Get the original files and reconstruct the cohort

You need R 4.6 or later, Python 3.10 or later, and about 3.1 GB for the two compressed archives. The source reader streams their contents. Run the following once from a terminal:

```sh
python -m pip install numpy pyproj
Rscript --vanilla source_workflow/01_get_fitrec.R
python -B source_workflow/02_source_counts.py --cohort
```

The R script downloads both archives from UC San Diego and validates their sizes and SHA-256 hashes. If you already have the original files, pass their directory as an argument to the R script and use `--raw-dir DIRECTORY` with the Python script. Use these records in accordance with the provider's academic-use terms. Source and attribution: [FitRec dataset](https://cseweb.ucsd.edu/~jmcauley/datasets/fitrec.html), [Ni, Muhlstein and McAuley (2019)](https://doi.org/10.1145/3308558.3313643).

The Python reader handles the provider's serialized Python dictionaries, matches the workout and metadata records, and checks the elapsed-time clock against UTC at each reading. For simultaneous records from the same runner it checks sensor traces, routes and time spans, then applies the session and location criteria described in Appendix A.1. It writes one aggregate CSV to `outputs/validation/source_counts.csv`. Expect 117,781 running workouts with confirmed clocks from 863 runners, and 115,036 eligible geolocated running sessions from 862 runners before weather linkage. Use `--output FILE` to change the aggregate output path. Omitting `--cohort` runs just the source linkage and clock checks.

## Generate the initial weather requests

To retain the event times and representative locations needed for weather linkage, run the preparation option when reconstructing the cohort. This also runs all the cohort checks:

```sh
python -B source_workflow/02_source_counts.py --prepare-weather-input
Rscript --vanilla -e "install.packages(c('data.table','ncdf4'), repos='https://cloud.r-project.org')"
Rscript --vanilla source_workflow/03_weather_requests.R
```

The R script downloads and verifies the ECMWF global land–sea mask. To use an existing copy, pass the path to the original NetCDF mask as the script's argument. On the 0.1° ERA5-Land grid, a land-dominant central cell contributes its own weather requests. A marine central cell contributes candidate land cells in the surrounding 3 × 3 neighbourhood. The ERA5 sensitivity calculation uses the nearest 0.25° cell. Dates cover all UTC hours needed to bracket each run. Dates at most 365 days apart share a request window. [ECMWF ERA5-Land documentation](https://confluence.ecmwf.int/pages/viewpage.action?pageId=140385202)

The script checks 6,903 initial ERA5-Land windows across 6,003 cells and 3,896 ERA5 windows across 3,321 cells. The aggregate counts go to `outputs/validation/initial_weather_requests.csv`. The session file and detailed weather plan contain dates and locations. Both stay under the Git-ignored directory `data/processed/local/`.

To obtain hourly temperatures and dew points, configure your own Climate Data Store account and accept the terms of each dataset. The following commands check the first cached weather response or request it directly from Copernicus. Increase or omit `--limit` to process more windows. Existing responses are validated before reuse:

```sh
python -m pip install netCDF4 cdsapi
python -B source_workflow/04_fetch_weather.py era5-land --limit 1
python -B source_workflow/04_fetch_weather.py era5-land --limit 1 --download
python -B source_workflow/04_fetch_weather.py era5 --limit 1 --download
```

The check verifies the variables, units, coordinates, hourly axis and finite temperature values, including responses in NetCDF or ZIP format. Downloaded files stay in `data/processed/local/weather/`. [CDS API setup](https://cds.climate.copernicus.eu/how-to-api), [ERA5-Land hourly time series](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries), [ERA5 hourly time series](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries).

## What these steps reproduce

The original FitRec files reproduce all 13 checks on linkage and cohort selection. Each of the 6,903 ERA5-Land and 3,896 ERA5 requests matches the dataset, variables, coordinates and date range used in the initial weather retrieval.

The scripts in this directory stop after the initial weather requests. For the paper, I made another 343 ERA5-Land requests to complete hourly coverage, linked the weather to individual runs and constructed the steady-speed segments, training histories and daily follow-up. To refit the models from individual records, you need those prepared panels and, for analyses of the development and confirmation samples, the original runner assignments. [Guide to the record-level analysis](../record_level/README.md).
