# Aggregate source data for the tables and figures

These are the printed entries for the tables and the aggregate data used to draw the figures. The [main README](../../README.md) lists the models that can also be fitted from the data in the repository.

## Tables

`table_catalog.csv` gives each table's English title and note. `table_cells.csv` contains the displayed text for main Tables 1 and 2 and Supplementary Tables S1–S29. Columns identify the table, panel, row, row type, starting column and number of columns the cell spans. An entry such as ≥100 retains the reporting bound. The text preserves the decimal places, confidence limits and footnote symbols shown in the paper.

`scripts/tables.R` rebuilds each table as an HTML page and a CSV. The HTML contains the column groups, panels and note. In the CSV, the original headings appear as rows within each panel. To recalculate the estimates for Table 2, use the [temperature models](../../scripts/temperature_models.R). Their output is checked against the printed values.

## Figures

`scripts/figures.R` redraws nine numerical figures in R using the files in `figures/`. Estimates and confidence limits in these files are aggregate plotting values. They can have more decimal places than the corresponding tables. For the paired analyses, whole-session heart-rate models, and recording or sensitivity models with a protected runner count, the released t-based limits use the published lower bound on the number of runners. This applies to Figures 3–5, S3, S5 and S6. All affected limits agree with the original at two decimal places. The model estimates retain their original values. The Figure S4 histogram reports the share of sessions in each one-degree bin to the nearest 0.1 percentage point. The PDF layouts can differ from the figures in the paper. Figure 1 and Figure S1 are supplied as PNG images.

| Figure | Distributed source | What the script displays |
| --- | --- | --- |
| 1 | `fig1_design.png` | Copies the study design diagram. |
| 2 | `figure2_hourly.csv`, `figure2_daily_temperature.csv`, `figure2_daily_participation.csv` | Hourly temperature and start-time profiles, session temperatures, and weekday and weekend running. The disclosed bin counts retain their released rounding. |
| 3 | `figure3_stages.csv`, `figure3_stages_curve.csv`, `figure3_checks.csv` | The heart-rate decomposition, its change with temperature and the alternative measurements and specifications. The floating bars show the differences between adjacent stages. |
| 4 | `figure4_paired.csv` | Within-run heart-rate contrasts at the 90th and 95th temperature percentiles. |
| 5 | `figure5_primary.csv`, `figureS2_curves.csv`, `figure2_daily_temperature.csv` | Four outcome curves, the fixed-hour comparison for heat avoided, and the distribution across released two-degree temperature bins. |
| 6 | [Activity models](../processed/activity_model/README.md) | Daily probability of recording a run with heart rate data, cumulative change in days with such a run, and changes in the probability of recording other activities. |
| S1 | `supp_fig1_map.png` | Copies the geographic summary as a single image. The repository includes no spatial source file for this map. |
| S2 | `figureS2_curves.csv` | Heat avoided, start-time shifts and within-day temperature rank. |
| S3 | `figureS3_recording.csv` | Run recording and availability of heart-rate measurements. |
| S4 | `figureS4_history.csv`, `figureS4_distribution.csv`, `figureS4_temperatures.csv` | Effect modification by recent exposure history, the exposure distributions and their percentile markers. |
| S5 | `figureS5_weekends.csv` | Weekday and weekend contrasts in heat avoided and heart rate. |
| S6 | `figureS6_sensitivity.csv` | Four sets of estimates across the twelve reported specifications. |

The temperature series use degrees Celsius. Heart rate is in beats per minute and speed in metres per second. Each script labels probability contrasts in the scale shown on its axis. Figures 1 and S1 have no underlying participant-level points in this package.
