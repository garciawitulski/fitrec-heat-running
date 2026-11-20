# Behavioural and heart-rate responses to heat in recreational runners: evidence from wearable data

Christian García-Witulski

I used FitRec records to study how recreational runners changed their training on hotter days and whether the average heart rate for a whole run captured the rise during steady running. This repository contains the R scripts and aggregate data for the paper.

## Run the analysis

Download or clone the repository, open a terminal in its top-level directory and run

```sh
Rscript --vanilla run.R --install
```

The analysis uses R and `data.table`. If the package is missing, the command installs it in `.R-library/`. If you already have it, run `Rscript --vanilla run.R`. All input files are included, so this second command also works offline.

## Results

The scripts write their results to `outputs/reproduced/`. Four analyses can be run directly from the aggregate data.

| Analysis | What the script calculates |
|---|---|
| [Temperature and running](scripts/temperature_models.R) | Contrasts at the 90th and 95th percentiles of daily maximum temperature relative to the median for all 11 outcomes in Table 2. Standard errors are clustered by runner. The four primary comparisons at the 95th percentile use Holm-adjusted p values. |
| [Heart rate within the same run](scripts/paired_heart_rate.R) | The difference between the heart rate response during steady running and the whole-session response, using three models fitted to the same 18,626 runs. |
| [Timing and speed](scripts/heart_rate_decomposition.R) | The predicted heart rate response when runners change their start time and speed, compared with the response at a fixed start time and speed. |
| [Participation and other activities](scripts/activity_models.R) | Seven models of daily activity, running around an isolated hot day and weekly activity. Figure 6b identifies the two activity estimates fitted here and the nine taken from the paper. |

For Table 2, `tables/manuscript/table2.csv` contains the calculated estimates and `validation/table2_printed.csv` checks them against all 81 values reported in the paper. The [table script](scripts/tables.R) uses the reported values to format the two main tables and 29 supplementary tables as HTML and CSV. The [figure script](scripts/figures.R) redraws nine figures from aggregate data and copies the study diagram and map. The activity analysis draws Figure 6. [Guide to output files](outputs/README.md).

## Individual records

Changing the sample or predictors, or resampling individual runners, requires the underlying records. The [source scripts](source_workflow/README.md) reconstruct the running cohort from the original FitRec files and prepare the initial weather requests. The [record-level scripts](record_level/README.md) start from linked session and runner-day data and also cover the subgroup and sensitivity analyses. The steps that link weather to individual sessions and build the analysis panels are outside the source scripts. [Guide to the data](data/README.md).
