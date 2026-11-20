# Data for the replication

The estimates are calculated from the model statistics in `processed/`. The tables and figures also use the reported values in `published/`.

| Directory | Contents and purpose |
|---|---|
| [processed/model_statistics](processed/model_statistics/README.md) | Whole-cohort cross-products for the eleven Table 2 models, clustered score products and the temperature percentiles used for comparison. |
| [processed/paired_model](processed/paired_model/README.md) | Cross-products for three models fit on the same sessions, including covariance between their scores. |
| [processed/decomposition_model](processed/decomposition_model/README.md) | Cross-products for the fixed-hour temperature gap and physiological model, covariance and cohort-level reference values. |
| [processed/activity_model](processed/activity_model/README.md) | Pooled matrices for the participation and activity models, together with the comparison values used for the published figure. |
| [published](published/README.md) | The reported table entries, titles and notes, and the data series plotted in the figures. |

The model files contain sums across the analysis cohort. Individual runs, runner identifiers, routes and timestamps are excluded. Where the paper gives a lower bound on the number of runners, the scripts calculate the range of confidence limits and p values consistent with that bound and the size of the parent cohort. The scripts in `scripts/prepare/` document how a researcher with access to the prepared records can construct the model inputs.

In the Table S12 reference data, `registered` identifies session models with additional calendar-year indicators. The participation model uses runner-by-year-month fixed effects in both specifications.

The original wearable records come from FitRec at UC San Diego, with academic-use and redistribution conditions set by the provider. Weather exposures were linked using Copernicus ERA5-Land and ERA5. [FitRec](https://cseweb.ucsd.edu/~jmcauley/datasets/fitrec.html) · [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries) · [ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries).
