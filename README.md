# Behavioural and heart-rate responses to heat in recreational runners: evidence from wearable data

Christian García-Witulski

I used FitRec records to study how recreational runners changed their training on hotter days and whether the average heart rate for a whole run captured the rise during steady running. This repository contains the R scripts and aggregate data for the paper.

## Reproduce the tables and figures

With R installed, you can reproduce the tables and figures in one run. [run.R](run.R) calls the scripts in order and uses the aggregate data included here.

1. On GitHub, choose **Code → Download ZIP** and extract the ZIP, or clone the repository.
2. Open a terminal in the folder that contains `run.R` and run this command once. It installs `data.table` if needed. If the package is already installed, the same command works offline.

   ```sh
   Rscript --vanilla run.R --install
   ```

3. When the terminal prints `Completed.`, open `outputs/reproduced/`. The folders below follow the numbering in the paper.

| Where to look | Tables (HTML and CSV) | Figures (PDF or PNG) |
|---|---|---|
| Main paper | `tables/manuscript/` · Tables 1 and 2 | `figures/manuscript/` · Figures 1–6 |
| Supplement | `tables/supplementary/` · Tables S1–S29 | `figures/supplementary/` · Figures S1–S6 |

For example, open `tables/manuscript/table_2.html` to read Table 2 or `figures/manuscript/figure_6_participation.pdf` to see Figure 6. The [output guide](outputs/README.md) also lists the model estimates and checks.

### What runs, in order

| Step | Script called by `run.R` | What it produces |
|---|---|---|
| 1 | [temperature_models.R](scripts/temperature_models.R) | Contrasts at the 90th and 95th percentiles of daily maximum temperature relative to the median for the 11 outcomes in Table 2. It writes the estimates to `tables/manuscript/table2.csv` and checks all 81 reported values. |
| 2 | [paired_heart_rate.R](scripts/paired_heart_rate.R) | Three models comparing the heart rate response during steady running with the whole-run response in the same 18,626 runs. |
| 3 | [heart_rate_decomposition.R](scripts/heart_rate_decomposition.R) | Estimates of how much of the heart rate response is associated with changes in start time and speed. |
| 4 | [activity_models.R](scripts/activity_models.R) | Seven models of daily and weekly activity, including running around an isolated hot day. This script draws Figure 6. Panel b combines two recalculated estimates with nine reported in the paper. |
| 5 | [tables.R](scripts/tables.R) | HTML and CSV versions of Tables 1–2 and S1–S29, using the reported table values included in the repository. |
| 6 | [figures.R](scripts/figures.R) | Figures 2–5 and S2–S6, redrawn from aggregate data. It also copies the study diagram (Figure 1) and the map (Figure S1) into the figure folders. |
| 7 | [check_core_results.R](scripts/check_core_results.R) | Checks the paired heart rate estimates against Table S13 and the decomposition against Figure 3. The results are in `validation/`. |

## Optional analyses with individual records

Changing the sample or predictors, or resampling individual runners, requires the underlying records. The [source scripts](source_workflow/README.md) reconstruct the running cohort from the original FitRec files and prepare the initial weather requests. The [record-level scripts](record_level/README.md) start from linked session and runner-day data and also cover the subgroup and sensitivity analyses. The steps that link weather to individual sessions and build the analysis panels are outside the source scripts. [Guide to the data](data/README.md).
