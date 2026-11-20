# Participation and other workouts

The analysis refits seven models using pooled normal equations and covariance matrices. The two daily activity models estimate the probability of recording any sport or an outdoor ride. The distributed-lag model estimates the probability of recording at least one run with heart rate data on a given day, using daily maximum temperatures from seven days before to three days after that date. For each complete ISO week, the four outcomes are the number of days with at least one recorded run with heart rate data, total running time, total running distance and the number of workouts across all sports. [Analysis script](../../../scripts/activity_models.R).

`pooled_matrices.csv` gives each distinct element of the within-runner design cross-product and the CR1-adjusted score covariance. `pooled_vectors.csv` gives the design–outcome cross-products. Daily models remove seasonal and weekday terms. Weekly models remove seasonal terms. The script solves the normal equations for the coefficients and applies the sandwich formula for their covariance. The calculation extends the same method used for Table 2 to the jointly estimated lag and temperature-band terms. [Formula and notation](../model_statistics/README.md).

In Figure 6a, each point shows the change in the probability of recording at least one run with heart rate data on that day, in percentage points. The line adds the daily changes from the hot day onwards. At day 7 it gives the change in the expected number of runner-days with at least one such run per 100 eight-day periods. A runner who records two runs on the same day contributes one day. The interval for the line accounts for the covariance between daily estimates.

In Figure 6b, filled points are the estimates fitted here for any sport and outdoor cycling. Open points show the other nine estimates, taken from Table S22 or the Results at their reported precision (`published_sports.csv`). The output file `figure6_activity_types.csv` identifies the source of each estimate.

The daily activity-type panel includes 729 runners. For the lag and weekly models, the manuscript reports at least 720 runners from that panel. The script calculates the confidence limits and p values over the resulting range of degrees of freedom. The weekly contrasts replace a day at 10–20 °C with a day in the stated band. The file `reported_activity.csv` gives the published point estimate and confidence limits for eight contrasts from the Results and Tables S19–S22. All 24 displayed fields agree across the possible confidence limits.

The preparation script takes four Parquet files as arguments. The input tables contain sessions, runner-days, workout metadata and daily weather for lagged exposures. It reconstructs workout dates, weekly totals and temperature lags. It checks that at least ten runners contribute to every released aggregate moment and to each nonempty overlap of eligibility and the two released activity categories. The generated files contain neither runner-level scores nor dates or locations. [Preparation code](../../../scripts/prepare/build_activity_statistics.R).

```sh
Rscript --vanilla scripts/prepare/build_activity_statistics.R sessions.parquet runner_days.parquet workout_metadata.parquet daily_weather.parquet data/processed/activity_model
```
