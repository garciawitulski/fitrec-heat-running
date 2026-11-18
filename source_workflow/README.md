# Rebuild the FitRec cohort from the original records

Run these commands from the repository root. This workflow starts with the two unmodified UC San Diego archives. It reproduces the FitRec source counts and the first row of supplementary Table S1, before weather linkage. A subsequent step prepares the initial Copernicus hourly weather requests.

## Get the original files and reconstruct the cohort

You need R 4.6 or later, Python 3.10 or later, and about 3.1 GB for the two compressed archives. The source reader streams their contents. Run the following once from a terminal:

```sh
python -m pip install numpy pyproj
Rscript --vanilla source_workflow/01_get_fitrec.R
python -B source_workflow/02_source_counts.py --cohort
```

The R script downloads both archives from UC San Diego and validates their sizes and SHA-256 hashes. If you already have the original files, pass their directory as an argument to the R script and use `--raw-dir DIRECTORY` with the Python script. Use these records in accordance with the provider's academic-use terms. Source and attribution: [FitRec dataset](https://cseweb.ucsd.edu/~jmcauley/datasets/fitrec.html), [Ni, Muhlstein and McAuley (2019)](https://doi.org/10.1145/3308558.3313643).

The Python reader handles the provider's serialized Python dictionaries, matches the workout and metadata records, and checks the elapsed-time clock against UTC at each reading. For simultaneous records from the same runner it checks sensor traces, routes and time spans, then applies the prespecified event-selection and geographic criteria. It writes one aggregate CSV to `outputs/validation/source_counts.csv`. Expect 117,781 running workouts with confirmed clocks from 863 runners, and 115,036 eligible geolocated running sessions from 862 runners before weather linkage. Use `--output FILE` to change the aggregate output path. Omitting `--cohort` runs just the source linkage and clock checks.

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

## Coverage

Running the steps above on the original FitRec release reproduced all 13 source and cohort checks. The resulting 6,903 ERA5-Land and 3,896 ERA5 request definitions matched the archived first-stage requests individually, including the dataset, variables, coordinates and date ranges. The response reader was tested against three original Copernicus responses from each product. No new hourly download was needed for these checks.

This source workflow ends with the initial weather acquisition. The published analysis also used 343 follow-up ERA5-Land requests after inspecting hourly availability, constructed exposures and physiological outcomes, and assembled training histories and daily follow-up. These later record-level data preparation steps are outside this workflow. Reproducing the development and confirmation samples requires the runner assignments used in the study. The repository's main README describes the analyses available from the included aggregate inputs and the optional models that use prepared records.
