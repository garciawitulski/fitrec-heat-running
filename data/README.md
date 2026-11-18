# Data included in the replication

The default analysis uses two kinds of input, both in the repository.

| Directory | Contents and purpose |
|---|---|
| [processed/model_statistics](processed/model_statistics/README.md) | Whole-cohort cross-products for the eleven Table 2 models, clustered score products and the temperature anchors. |
| [processed/paired_model](processed/paired_model/README.md) | Cross-products for three models fit on the same sessions, including covariance between their scores. |
| [processed/decomposition_model](processed/decomposition_model/README.md) | Cross-products for the fixed-hour temperature gap and physiological model, covariance and cohort-level reference values. |
| [processed/activity_model](processed/activity_model/README.md) | Pooled matrices for the participation and activity models, together with the comparison values used for the published figure. |
| [published](published/README.md) | Printed cells, table metadata and aggregate figure series. These supply the manuscript displays and checks on the independently calculated contrasts. |

The joint products combine observations across the analysis cohort. They contain no session rows, runner identifiers, routes or timestamps. Inference uses the published range of runner counts when the exact number is protected. The scripts in `scripts/prepare/` document how a researcher with access to the prepared records can construct the model inputs.

The original wearable records come from FitRec at UC San Diego, with academic-use and redistribution conditions set by the provider. Weather exposures were linked using Copernicus ERA5-Land and ERA5. [FitRec](https://cseweb.ucsd.edu/~jmcauley/datasets/fitrec.html) · [ERA5-Land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries) · [ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries).
