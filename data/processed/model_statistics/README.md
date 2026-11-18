# Aggregate inputs for Table 2

`moments.csv` contains eight pooled numerical values for each of the 11 outcome models. `temperature_percentiles.csv` provides the spline knots, the three comparison temperatures and the limits for the plotted range. [Estimation script](../../../scripts/temperature_models.R).

## Estimating a temperature contrast

For each outcome, the preparation first removes fixed effects by subtracting group means from the temperature spline, the nuisance covariates and the outcome. The groups are runners for the session models and runner–year-month combinations for the participation model. The session models also adjust for season and weekday, and the participation model adjusts for weekday. The preparation then residualises the two temperature-spline columns and the outcome against the remaining covariates.

Let $X$ be the two residualised spline columns, $y$ the residualised outcome and $n$ the number of observations in that model. The file gives the upper triangle of $A=X^\top X/n$ in `xx11`, `xx12` and `xx22`, and the two elements of $b=X^\top y/n$ in `xy1` and `xy2`. The public script solves $A\hat\beta=b$. For each requested temperature it subtracts the spline basis at the median from the basis at that temperature and multiplies the resulting contrast vector by $\hat\beta$. Residualising against the other predictors gives the same two temperature coefficients as fitting the full fixed-effects model.

## Clustered intervals and p values

For runner $g$, let $s_g=X_g^\top(y_g-X_g\hat\beta)$. The last three columns give the upper triangle of $M=f\sum_g s_gs_g^\top/n^2$, where $f=[G/(G-1)][(n-1)/(n-k)]$, $G$ is the number of runners and $k$ is the rank of the full within-group design. The public script calculates $V=A^{-1}MA^{-1}$ and the model-scale standard error $\sqrt{c^\top Vc}$ for each contrast $c$. This aggregate score covariance is specific to the fitted model. Deriving it for a new specification requires the underlying records.

The paper reports exact runner counts for the main cohort and participation panel. For the two restricted samples, it gives at least 810 runners with admissible heart-rate readings and at least 660 runners with steady-speed segments in the main sample. The respective parent cohorts contain 820 and 669 runners. The script calculates the limits of the t-based inference over those public bounds. The point estimate and clustered standard error are uniquely determined for all 22 contrasts. The columns ending in `_min` and `_max` give the range of each confidence limit or p value. Equal endpoints identify a unique value. Both ends of each range reproduce the corresponding printed field, including Holm adjustment across the four P95 primary comparisons.

## Preparation and checks

The builder accepts authorized Parquet files with one row per session, steady-speed session and runner-day, respectively. Run it from the repository root to generate the two aggregate CSVs. [Aggregate builder](../../../scripts/prepare/build_model_statistics.R).

```sh
Rscript --vanilla scripts/prepare/build_model_statistics.R sessions.parquet steady_segments.parquet runner_days.parquet data/processed/model_statistics
```

The builder verifies that each distinct pattern of inclusion across session models pools at least ten runners. For each numerical moment it also checks that at least ten runners contribute a nonzero quantity. It writes only the pooled model statistics and cohort-wide temperature percentiles. The estimation script checks all 81 displayed fields against Table 2.
