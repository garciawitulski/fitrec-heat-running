# Behavioural and heart-rate responses to heat in recreational runners: evidence from wearable data

**Author** Christian García-Witulski

Download or clone this repository to reproduce the main temperature contrasts, the paired heart-rate comparison and selected activity models. The scripts also rebuild the two manuscript tables, the 29 supplementary tables and the figure set.

## Run the analysis

Open a terminal in the downloaded repository and run:

```sh
Rscript --vanilla run.R --install
```

You need R. The `--install` option installs the one required package, `data.table`, in `.R-library/` inside the clone if it is missing. If `data.table` is already installed, `Rscript --vanilla run.R` works offline. Results are written to `outputs/reproduced/`. All inputs for that run are included in the download.

## What the scripts calculate

| Analysis | Inputs and calculation | Script |
|---|---|---|
| Table 2 and the main temperature curves | Fit the temperature coefficients of eleven within-runner models from whole-cohort sums of products. Calculate the 22 P90 and P95 contrasts, clustered standard errors and Holm-adjusted p values. | [Table 2 models](scripts/temperature_models.R) |
| Steady-speed and whole-session heart rate | Fit three regressions on the same 18,626 sessions. Estimate the difference between responses with its paired standard error. | [Paired heart-rate analysis](scripts/paired_heart_rate.R) |
| Timing, speed and the heart-rate response | Fit the fixed-start-hour and physiological models, then combine them with the main models to calculate the heart-rate response at fixed and observed timing and speed. | [Decomposition](scripts/heart_rate_decomposition.R) |
| Participation and other sports | Refit seven models for any sport, outdoor cycling, lagged daily running and weekly activity. The Figure 6 activity file labels the two refitted sport series and the nine series drawn from reported contrasts. | [Activity models](scripts/activity_models.R) |
| All manuscript and supplementary tables | Render the released cells at their reported precision, including counts disclosed as lower bounds. Each table is available as HTML and CSV. | [Table rendering](scripts/tables.R) |
| Figures | Draw Figure 6 from the activity analysis and redraw nine other numerical figures from the released aggregate series. The study diagram and the geographic map are included as rendered images. | [Figure rendering](scripts/figures.R) |

The regressions use small matrices of sums over the cohort. These retain the information needed to fit the registered temperature coefficients after adjustment for runner effects and other covariates. Cluster score products supply the uncertainty estimates. For cohorts whose exact runner count is protected, the scripts calculate the range of confidence limits and p values consistent with the disclosed counts. All 81 printed values in Table 2 agree across those ranges. Changing the cohort, adding predictors or resampling individual runners requires the underlying records. [Data and methods](data/README.md).

## Find an output

Within `outputs/reproduced/`:

- **`tables/manuscript/table2.csv`** contains the recalculated estimates. **`validation/table2_printed.csv`** compares them with the printed values. **`tables/manuscript/table_2.html`** displays Table 2 at publication precision.
- **`models/paired_heart_rate.csv`** and **`models/heart_rate_decomposition.csv`** contain the paired estimates and physiological decomposition.
- **`tables/supplementary/activity_models.csv`** gives the activity estimates. **`tables/supplementary/table_S1.html`** through **`table_S29.html`** are the supplementary displays.
- **`figures/manuscript/`** and **`figures/supplementary/`** contain the figures. **`validation/`** records the numerical checks.

## Working from individual records

The source workflow starts from the two original FitRec files available from UC San Diego. It reconstructs the cohort before weather linkage and prepares the first weather requests. To refit the remaining models at the individual-record level, use the linked session and runner-day panels. Reproducing the development and confirmation samples also requires their original assignments. [Original data and source checks](source_workflow/README.md) · [Record-level analyses](record_level/README.md).
