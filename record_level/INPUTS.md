# Inputs for the individual-record models

FitRec links workout metadata to readings from wearable devices. The session panels use local time to assign exposures and calculate whole-session heart rate and heart rate during steady-speed segments. The daily panel assigns each runner the weather at their usual grid cell on each day of follow-up, including days without a run. Daily maximum temperatures require hourly coverage across the entire local calendar day. [FitRec source and data terms](https://cseweb.ucsd.edu/~jmcauley/datasets/fitrec.html) · [Original source workflow](../source_workflow/README.md)

The [inventory](input_checksums.csv) lists the input names used in the code and their SHA-256 checksums. Place the files under `data/record_level/`, or set `FITREC_INPUT_DIR` to another root. For example, `panels/sessions.parquet` is the prepared session file under that directory. The licensed records must be obtained separately.

If you already have the files, you can point to them in place. Make a CSV with the columns `input` and `path`, then set `FITREC_INPUT_MAP` to its location. Put one input on each row. Relative paths in the second column start from the directory containing your CSV. An input omitted from the CSV uses its inventory path under `FITREC_INPUT_DIR`. Keep the boundary shapefile and its companion files together under the same base name.

```csv
input,path
panels/sessions.parquet,linked/running_sessions.parquet
panels/runner_days.parquet,linked/daily_followup.parquet
```

## Data and model inputs

| Input group | Use in the analysis |
|---|---|
| `panels/` | One row per run in `sessions.parquet` and one per runner-day in `runner_days.parquet`. The steady-speed file contains sessions with eligible segments. The paired heart-rate panel contains both heart-rate outcomes and their common adjustment variables. The other files describe screened heart-rate readings, run composition and coverage at the ends of the weather windows. |
| `covariates/` | Daylight, weather, recording resolution and clock-change analyses. Covariates cover the main session cohort, runner-day panel or steady-speed cohort, according to the file name. |
| `cohort/` | Workout and clock eligibility, geographic linkage, session weather, hourly ERA5-Land series, prior workouts and the held-out exercise-history sample. |
| `history/` | Design matrices, exercise histories, lagged exposures and row indices for the prior-heat models. The four directories distinguish the daily maximum, session temperature, weighted history and prior exercise. |
| `geography/` | Adaptive regions for the aggregate map and Natural Earth country boundaries at 1:10 million scale, version 5.1.1, for Table 1. |
| `settings/` | Spline knots, reference speed and other fixed model parameters for the session and prior-heat analyses. |
| `references/` | Aggregate estimates, grids of temperatures at which curves are evaluated, published map totals, and contrasts for checking model results before rounding. The paper's printed values are distributed separately in the repository's `data/processed/` directory. |

ERA5-Land supplies the main hourly temperature series. A sensitivity analysis uses ERA5. The individual-record models require the linked panels and measures listed above, including the assignments for analyses in the held-out sample. [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries) · [ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries) · [Natural Earth](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/)

## Model and field identifiers

The identifiers in the design matrices also serve as keys when matching the aggregate reference estimates. In the prior-heat models, `AB` denotes behaviour and whole-session heart rate in relation to the daily maximum, `AH` denotes fixed-speed heart rate in relation to session temperature and `EW` denotes exponentially weighted exposure history with a seven-day half-life. The digits `07` and `14` give the length of the exposure window in days. The outcomes are heat avoided (`B1`), median speed (`B2`), log duration (`B3`), whole-session heart rate (`H1`), fixed-speed heart rate (`H2`) and local start time (`START`).

The weather columns whose names begin with `a3_` contain precipitation, wind or radiation. Units appear in each column name. The `a4b_` and `a4c_` identifiers in the heart-rate reference tables distinguish the sampling-resolution and paired subgroup models. The code retains these names where it joins to the supplied files.
