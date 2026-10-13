# Statistical methods of the missingness pipeline

This document describes what the code computes, in the order of the notebook. Every statement was written
from the code in `R/` (file names are given in brackets). Nothing here is a result: no real-data result
exists yet, so no estimate is quoted. Where a number is a convention from outside the code (for example
a rule of thumb for reading an effect size) it is marked **[convention]**.

---------------------------------------------------------------------------------------------------

## 0. Design in one paragraph

Fitbit activity, sleep and heart-rate data of four cohorts (Insomnia, Diabetes, Healthy, Prostate) are reduced to
one row per person per day. A **day is "observed" if a row exists**. Missingness is studied only **between**
a person's first and last record (window (a)), so what is measured are *interruptions that were followed by
resumption*; leading missingness and permanent drop-out are invisible by construction. The analysis has three
layers: (1) descriptive person-level missingness and its temporal structure; (2) a mixed-effects logistic
model of **gap onset** given the previous day's values; (3) two exploratory models that use **gap length**
(which is only known after a gap ends and is therefore never a predictor of onset). No disease-label prediction
and no mixture model is part of the analysis.

## 1. Definitions [`gaps.R`, `etl.R`, `persons.R`]

* **Window**: first to last observed date of a person *in that stream*; `window_days = last - first + 1`.
* **Gap**: a run of >= 1 consecutive calendar days with no row, lying strictly between two observed days.
  **Gap length** = number of missing days in the run. Gaps are found from the sorted observed dates; no
  calendar grid is built.
* **Selection already applied upstream** (not re-applied; only *checked* and reported in `selection_check`):
  >= 70% of window days observed, longest gap <= 14 days, >= 21 days. Consequences: gap lengths lie in 1-14 days,
  long-gap behaviour is under-represented, and all results are conditional on being a good wearer.
* **Cohort overlap**: a person in more than one cohort is, by default (`cfg$overlap$policy = "exclude"`), left
  out so that every person counts once (the cohort tests assume independent groups). `"priority"` keeps the person
  in the first cohort of Prostate > Diabetes > Insomnia > Healthy. Your Cell 8 reported 383 of 8,926 people in
  more than one cohort.
* **Covariates**: age in whole years **at the person's first record** (not at the run date), scaled as
  `age10 = (age - 60) / 10`; sex at birth (Female / Male / Other); race grouped from the race name; Charlson
  Comorbidity Index (CCI) groups 0, 1-2, 3-4, >=5; age groups <=44 (labelled "18-44"), 45-64, >=65.
  CCI is used only as a **stratifier** in descriptive plots, never as an adjuster (it contains age points and the
  index disease).

## 2. Person-level missingness: summaries, tests, effect sizes [`gaps.R`, `summarise.R`]

**Unit of analysis = the person** (one row per person and stream), so people with long follow-up do not dominate
and the tests are not pseudo-replicated over days.

Per person and stream: `n_observed`, `window_days`, `n_missing = window_days - n_observed`,
`pct_missing = 100 * n_missing / window_days`, `n_gaps`, `gap_max`, mean / median gap length, number of gaps
>= 2 and >= 7 days. An internal check stops the run if the missing days implied by the window differ from the
sum of the gap lengths (duplicated or unsorted dates).

**Descriptive table** (`cohort_missingness`): per stream x cohort, median (Q1-Q3) of follow-up days, record days,
`% missing`, longest gap and gaps per person; percentage of people with any gap and with a gap >= 7 days.

**Overall test** (`cohort_tests_overall`): for each of `pct_missing`, `gap_max`, `n_gaps` and each stream, a
**Kruskal-Wallis rank-sum test** across the four cohorts. Columns: `chisq` (the KW statistic, chi-square
distributed with `df = number of cohorts - 1` under H0), `df`, `p`. No effect size is stored for the omnibus test.

**Pairwise tests** (`cohort_tests_pairwise`): for every pair of cohorts, a **Wilcoxon rank-sum
(Mann-Whitney) test** (normal approximation, `exact = FALSE`). Columns: `n1, n2, median1, median2`, `p_raw`
and `p_holm` (**Holm adjustment over the pairs within one stream x variable**), and the effect size
**Cliff's delta**:

    delta = P(X > Y) - P(X < Y) = 2 U / (n1 n2) - 1 ,   range -1 ... +1

(U = Mann-Whitney statistic of group 1). Positive delta: group 1 tends to have larger values. delta is the
probability-of-superiority scale; it is not a difference in medians. **[convention]** Romano et al. label
|delta| < 0.147 negligible, < 0.33 small, < 0.474 medium, otherwise large. The Wilcoxon statistic itself is not stored;
its p value and delta are.

**Share of missing days** (`missing_fraction_models`): a **quasi-binomial GLM** (logit link) on
`cbind(n_missing, n_observed) ~ cohort` (reference level **Healthy**), unadjusted and adjusted for `age10` and
`sex` (sex "Other" is dropped in the adjusted model). Each person is weighted by the number of days in the window
(binomial denominator); the quasi-binomial dispersion (`dispersion` column) absorbs extra-binomial variation
between people. Effect size: **odds ratio** (`or`, `or_low`, `or_high`) for the odds that a day inside the window is
missing, versus Healthy. The CI in the table is `exp(estimate +/- 1.96 se)`; `p` is the model's t-test p value
(summary of the GLM with estimated dispersion). Estimates are *marginal over people* in the sense of a GLM
with a scale parameter; no person random effect is used here.

### 2b. Person-level % missing by race, sex, CCI group and age group [`summarise.R`]

For each stream and each stratum (sex, age group, race, CCI group; age is at the first record) the person-level `pct_missing` is
shown by cohort. Two families of rank tests are stored in `missing_pct_by_stratum_tests`, all on one value per person:
* `groups_within_cohort` (the one drawn on the figures): do the **levels of the stratum** differ within a cohort? **Kruskal-Wallis**
  for >= 3 levels (statistic H = `statistic`, `df`, effect size **epsilon-squared = H / (n - 1)**, range 0-1), **Wilcoxon rank-sum**
  for 2 levels (statistic W, effect size **Cliff's delta**, first level vs second). Levels with fewer than 5 people in a cohort are left out;
  fewer than 2 remaining levels: no test (`test = "none"`). Holm adjustment over the cohorts of one stream x stratum (`p_holm`).
* `cohorts_within_level`: do the **cohorts** differ within each level (the question of the earlier cohort figures)? Same tests, Holm over the levels.
These are unadjusted descriptive comparisons (age, sex and cohort are not adjusted for each other here; the adjusted comparison is the
quasi-binomial model above). Per-cell summaries (N, mean, SD, median, quartiles, P5, P95) are in `missing_pct_by_stratum`.

## 3. Temporal structure at the day level [`temporal.R`]

All uncertainty here is a **person-level (cluster) percentile bootstrap** (`cfg$boot$B = 200` replicates,
seeded), because days of one person are not independent. Counts come from the gap table and arithmetic on each
person's window.

1. **Onset vs continuation** (`transition`), pooled per stream x cohort:
   * `p_onset` = P(missing tomorrow | observed today) = sum(gap onsets) / sum(observed days that have a next day);
   * `p_continue` = P(missing tomorrow | missing today) = sum(missing days followed by a missing day) / sum(missing days);
   * `p_missing` = sum(missing days) / sum(window days);
   * `persistence_ratio` = `p_continue / p_missing`. If missing days were scattered independently the ratio is 1;
     > 1 means missing days come in runs. Each with a bootstrap 95% interval (`_low`, `_high`).
2. **Study week and calendar year** (`study_week_rates`, `calendar_year_rates`): pooled missing share
   (sum missing / sum window days in that week or year) per cohort with bootstrap 95% intervals; only groups with
   >= 30 contributing people are reported. No test is computed. Because the first and last day of every window are
   observed by construction, missingness is slightly under-estimated in the first and last week.
3. **Day of week** (`weekday_rates`, `weekday_tests`, `weekday_odds_ratios`): per weekday, the pooled missing share
   with bootstrap intervals. Formally, a **binomial GLMM** (`lme4::glmer`, logit, `nAGQ = 0`, bobyqa) on
   `cbind(missing days, observed days of that weekday) ~ ... + (1 | person_id)`, in three nested versions
   m0: cohort; m1: + weekday; m2: weekday x cohort (unadjusted and adjusted for age and sex).
   Tests are **likelihood-ratio tests** (statistic `chisq = max(0, 2 (logLik_big - logLik_small))`, `df` = difference in
   parameters, chi-square reference): "weekday effect" (m1 vs m0; 6 df) and "weekday x cohort" (m2 vs m1;
   6 x (number of cohorts - 1) df). Effect sizes in `weekday_odds_ratios`: odds ratio of a day being missing for each
   weekday **versus Sunday** (reference = first level), conditional on the person's random intercept, with Wald 95% CI.
4. **Attained age (continuous)** (`age_rates`, figure `missing_by_age`): exact birth dates are not kept, so the age at a day is
   `age at first record + 0.5 + elapsed days / 365.25` (uncertain by about +/- 0.5 year). Missing days and window days are counted
   by whole year of age (the same thresholds for both, so the days of a person add up exactly), and each age pools the days at ages
   a - 2 ... a + 2 (a 5-year moving window, `cfg$strata$age_half_width`). Pooled missing share per cohort with a person-level cluster
   bootstrap interval (shaded band); ages with fewer than 30 people are not drawn. No test.

**Statistics saved for these figures** (stage `missing_figures` or `temporal`): for study week (per 4 weeks), calendar year (per year) and age (per 10 years, one row per
whole year of age, not the moving window), `missing_trend_tests` holds, per stream and cohort, the log-odds slope (`estimate`), its robust
standard error, `z`, `p`, the odds ratio per unit with 95% CI and `p_holm` (Holm over the cohorts of one stream x variable), plus a row for all cohorts
with cohort as a covariate; `missing_trend_heterogeneity` holds the Wald chi-square test (`chisq`, `df = cohorts - 1`, `p`, `p_holm`) that the slope
differs between cohorts. They come from **pooled logistic regressions of the aggregated counts**, `cbind(missing days, observed days) ~ x`, with
**person-clustered sandwich covariance** (`sandwich::vcovCL`), so days of one person are not treated as independent; the trend is **linear on the logit
scale** and curvature is not tested (the figure shows the actual shape). Cohorts with fewer than 30 people or fewer than 2 distinct values are reported as
`not_estimable`. For the weekday: `missing_weekday_tests` has the Wald chi-square test of the 6 weekday terms per cohort (Holm over cohorts), for all
cohorts with cohort as a covariate, and the weekday x cohort interaction (6 x (cohorts - 1) df); `missing_weekday_odds_ratios` the odds ratio of each weekday
versus Sunday per cohort with 95% CI. (The mixed-model version in `weekday_tests` needs the temporal stage.) For the stratified % missing
figures, `missing_pct_by_stratum_models` adds quasi-binomial odds ratios (level vs reference, within each cohort, days-weighted) with the omnibus F test of the stratum.

All four "missing by ..." figures (study week, year, weekday, age) use one configuration: a line per cohort with a shaded 95% bootstrap band.

## 4. Descriptive statistics of the metric values [`descriptive.R`, `plots.R`]

Per person and metric: mean, SD and coefficient of variation (CV% = SD / mean x 100, defined only for mean > 0). The
summaries (`metric_summary*`) are then taken **over the person means**: N, mean, SD, median, Q1, Q3, P5, P95 by
cohort and by cohort x stratum (sex, age group, race, CCI group). These tables contain **no test**. The figures
print a Kruskal-Wallis p value between cohorts for each panel (and, in stratified figures, within each stratum
level); these are descriptive annotations, not corrected for the number of panels. (A Dunn-test helper with Bonferroni
brackets, `plot_single_metric_violin`, exists but is not called by the pipeline.) Columns that are > 90% NA or constant
were dropped by the ETL and are listed in `dropped_columns.tsv`.

## 5. Gap-onset model (the main inferential analysis) [`gaps.R`, `predictors.R`, `onset_model.R`]

### 5.1 Data structure (risk set) and outcome

One row per **observed day d** of a person, except the person's last observed day. Outcome
`y_ge{k} = 1` if a gap of **at least k days starts on day d+1** (`next observed date - d - 1 >= k`); otherwise 0.
Primary definition k = 1 ("the next day is missing"); alternatives k = 2, 3, 7. The predictors on row d are
the values of **day d itself**, i.e. the day before the possible gap; there is no information from d+1 or later
in the pairing (a unit test checks the pairing against an independently built calendar). The risk set is "common":
a row is kept only if **all** predictors of that stream are non-missing, so univariate, multivariate and joint
models use identical rows.

### 5.2 Predictors [`predictors.R`, `config.R`]

Predictors (name: source, transform): `act_steps`: steps, log(1 + x); `act_sedentary`: sedentary minutes, identity;
`act_mvpa`: MVPA minutes, log(1 + x); `sleep_inbed`: minutes in bed, log(1 + x); `sleep_eff`: sleep efficiency,
identity; `hr_cover`: total heart-rate minutes recorded (this measures **watch wear / coverage**, not behaviour),
identity; `hr_elev`: minutes in fat-burn + cardio + peak zones, log(1 + x). Sleep and heart-rate models use
the stream's own predictors plus the three activity predictors (the question "do gaps in stream S follow low
activity?").

For each predictor x (transformed) of person i on day d, per stream:

* **Within-person** `w = (x - mean_i) / sd_w`, `mean_i` = the person's mean over all their observed days,
  `sd_w = sqrt( sum (x - mean_i)^2 / (N_obs - N_persons) )` = one **pooled within-person SD** (degrees-of-freedom
  corrected). "How unusual was yesterday for this person."
* **Between-person** `b = (mean_i - mean of person means) / SD of person means`. "Is this person generally high or
  low." Keeping `w` and `b` as separate terms keeps within- and between-person effects apart.
* (Sensitivity only) **trailing** `tw = (x - mean of the previous 14 observed days) / sd_w`, requiring >= 7 of those days.

Known limitation: `mean_i` and `sd_w` use the whole follow-up, so the baseline for day d contains later days
(look-ahead in the baseline, not in the pairing). The trailing-baseline sensitivity removes it.

### 5.3 Model

Binomial mixed model, logit link, **random intercept per person**:

    logit P(y_ge{k}(i,d) = 1) = b0 + sum_p [ beta_p w_p(i,d) + gamma_p b_p(i) ] + (adjusters) + u_i ,   u_i ~ N(0, sigma_u^2)

fitted with `lme4::glmer(nAGQ = 0, optimizer bobyqa)` (default; `glmmTMB` is an option). `nAGQ = 0` is a faster
approximation (fixed effects are optimised inside the penalised least-squares step); it is about 10x faster than
glmmTMB and ~100x faster than `nAGQ = 1`, and a test requires its primary estimate to agree with glmmTMB within 0.05
on the log-odds scale. **Adjusted** models add: cohort, `age10`, sex, weekday of the predicted day, a natural
cubic spline of study day (3 df), calendar year of the predicted day (factor), and the gap-history terms
(log(1 + length of the most recent completed gap), any prior gap, log(1 + days since the last return, capped at 60)).
Adjusters that are constant in the data (for example sex in the all-male Prostate cohort) are dropped
automatically. The history terms use only gaps that ended on or before day d (they are earlier gaps, not the
gap being predicted).

Model grid per stream: every single predictor, every multi-predictor set and the joint set, each unadjusted and
adjusted (`role = exploratory`); **primary** = joint predictors, adjusted, k = 1; **model D**
(`secondary_threshold`) = the primary model with k = 2, 3, 7.

### 5.4 Effect size and test statistic

* Coefficients `beta` (log-odds) with standard error `se` from the observed Hessian; **Wald z = beta / se**, two-sided
  `p` from the standard normal (`z` itself is not stored: compute `estimate / se`); 95% CI = `beta +/- 1.96 se`.
* **Odds ratio**, as stored in `or`, `or_low`, `or_high`: for **within-person (`w_`) terms OR = exp(-beta), "per 1
  SD *below* the person's usual"** (OR > 1: a lower-than-usual day goes with higher odds that a gap starts the next day);
  for `b_` terms and adjusters OR = exp(beta) per +1 unit (SD for `b_`). These are odds ratios, **not** risk ratios, and are
  **conditional** on the person's random intercept.
* **Primary hypothesis**: `w_act_steps` in the joint adjusted k = 1 model, one test per stream, **Holm-adjusted across
  the streams** (`p_holm` in `primary_test`). This is the only pre-specified test. Everything else (model D, the other
  specifications, sensitivity, B and C) is secondary or exploratory and is **not** multiplicity-adjusted.
* **Heterogeneity between people**: `sigma2` (random-intercept variance, latent scale) and the **median odds ratio**
  `MOR = exp( sqrt(2 sigma2) x z_0.75 )` (typical OR between two random people with the same covariates).
* **Convergence and singularity** are reported, not hidden: `status` is `ok`, `not_converged`, `not_estimable` (fewer than 100 events
  or 30 people with events) or `failed`; a singular fit (person variance ~ 0) is flagged (`singular`, note) and counts as converged.
* **Collinearity** (`vif`): variance inflation factor of each `w_` predictor from the correlation matrix, `diag(R^-1)`.
* **Cohort heterogeneity** (`cohort_interaction`): LRT of the primary model against the same model plus
  `w_term x cohort` (`chisq`, `df`, `p`); `cohort_specific`: the primary term fitted **within each cohort** separately
  (cohorts without enough events are reported as not estimable).
* **Absolute risk** (`risk_unadjusted`): from the **unadjusted joint** k = 1 model, the population-average probability
  that a gap starts tomorrow when yesterday's `w` for the primary term is -1, 0, +1 SD (all other predictors at their average),
  `mean over u of plogis(b0 + b1 a + u)` with `u ~ N(0, sigma2)` integrated on a 201-point grid. The interval (`risk_low`, `risk_high`)
  comes from 2,000 draws of (`b0`, `b1`) from their joint normal distribution; `sigma2` is held fixed, so it reflects
  intercept and slope uncertainty only.

Note: in every adjusted model `cohort` is coded with **Healthy as the reference** (`cfg$reference_cohort`; changing it makes all stored model results stale).

## 6. Sensitivity analyses of the primary effect [`sensitivity.R`]

Same question (does a lower-than-usual prior day raise the odds of a gap?) under different assumptions; one row per analysis
in `sensitivity`, with `estimate, se, conf_low, conf_high, p, or, or_low, or_high`. The same fits also store the estimates of **every**
within-person predictor of the model (`w_` / `tw_` terms, not only the tested one) in `sensitivity_all_terms`; these are what the
per-stream forest plot `sensitivity_forest_<stream>` draws (the spline and random-slope analyses are tests and are not in that table):

| analysis | what changes |
|---|---|
| `main` | the primary model itself |
| `trailing_baseline` | `tw_` (previous 14 observed days) instead of `w_`: no look-ahead in the baseline; the first observed days lack a baseline and drop out |
| `exclude_partial_wear` | drops rows whose previous day looks like partial wear (steps < 500 or sedentary >= 1380 min; sleep: minutes in bed < 60; heart rate: total minutes < 1200) |
| `first_gap_only` | only days up to and including each person's first gap onset (at most one event per person; no gap history) |
| `main_sleep_only` | sleep only: a night counts as observed only if a main sleep exists (naps alone do not count) |
| `conditional_logit` | `survival::clogit`, person as stratum, Efron approximation, persons with >= 1 event only; time-invariant terms (cohort, age, sex, `b_`) cancel. No random-effect assumption |
| `cluster_robust_logit` | pooled logistic regression without random effect, person-clustered sandwich standard errors (`sandwich::vcovCL`) |
| `spline_vs_linear` | LRT replacing the linear `w` term by a natural spline with 3 df (2 extra df); small p: non-linear. Reported in `p` and the note (`chisq`, `df`), no OR |
| `random_slope` | LRT adding an uncorrelated person-specific slope for `w`; because the variance is tested on its boundary the p value is **halved** (50:50 mixture of chi-squares, 1 df). Small p: the effect differs between people |

Estimates from different estimators are comparable on the log-odds scale of the same term; the conditional logit
estimates a within-person effect and need not equal the random-intercept estimate exactly.

## 7. Gap-length models B and C (exploratory) [`gap_models.R`]

Both use **one row per gap** (all gaps; people with several gaps contribute several rows, hence a person random intercept).
They are **conditional on a gap having occurred**, so they describe gap-prone people, and they do not estimate the
effect of anything on gap onset. Because of the upstream filter (longest gap <= 14 days), the "15+" length group is
expected to be empty; sparse groups (< 20 gaps) are merged into the next shorter group.

**Model B - "is activity lower after longer gaps?"** Linear mixed model (`lme4::lmer`, REML):
`outcome ~ log2(gap length) [+ prior-day deviation] [+ cohort + age10 + sex] + (1 | person_id)`.
Outcome: the within-person z-score of log(1 + steps) on the **first full day back** (`next observed date + 1`, primary), or on the first
day back (sensitivity; the watch is often put on mid-day). Four specifications: unadjusted; adjusted; adjusted + the pre-gap
deviation (change from before to after); first-day-back outcome adjusted. **Effect size**: `estimate` = change in the
outcome per **doubling** of the gap length, in within-person SD units, `se`, `conf_low/high = est +/- 1.96 se`, and `p` from a
**normal approximation** to the t statistic (`2 pnorm(-|est/se|)`; lme4 gives no p value; with 10^4-10^5 gaps the difference to a t reference is negligible).
Also reported: `n_gaps`, `n_persons`, outcome mean/SD, `sigma_person`, `sigma_resid`, `singular`, `converged`, `status`.
The figure draws the fitted line (95% CI from the fixed-effect covariance, covariates at their sample mean) against the observed group means
with clustered 95% intervals.

**Model C - "does low activity the day before a gap predict a longer gap?"** Cumulative-logit (proportional-odds)
mixed model (`ordinal::clmm`, Laplace; sparse random-effect matrices) of the ordered length group (1, 2, 3, 4-7, 8-14, 15+)
on the prior-day within-person steps z-score (unadjusted, and adjusted for cohort, age, sex), with a person random
intercept: `logit P(Y <= j) = theta_j - beta w - ... - u_i`. It is run only if the **proportional-odds assumption** is
not rejected: a likelihood-ratio test of the model without random effect against the same model with a non-proportional (nominal) effect of `w`
(`po_p`, rejected if `po_p <= 0.05`; status `po_rejected`, no estimate; model D is then the analysis to report).
**Effect size**: `or = exp(-beta)`, the odds of a **longer gap per 1 SD below the person's usual**, with Wald CI and z-test p.
If the random-intercept variance is on the boundary or the Hessian is unusable, the model **without** the random effect is used and
**all uncertainty comes from a person-level cluster bootstrap** of the full parameter vector (`ok_bootstrap_se`; replicate b is
seeded with `seed + b`, so results do not depend on the number of workers); the SE is then the bootstrap SD and the interval
`est +/- 1.96 se`. Predicted probabilities of each length group at w = -1, 0, +1 are **population-average** (averaged over a seeded
sample of 500 gaps' covariate profiles and over the random intercept), with 95% intervals as percentiles of 300 parameter draws
(multivariate normal from the Hessian, or the bootstrap draws in the fallback).

## 8. Where each effect size and test statistic is stored

| analysis | effect size | test statistic / inference | table (`results/tables/`) |
|---|---|---|---|
| cohort comparison of missingness (omnibus) | none stored | Kruskal-Wallis `chisq`, `df`, `p` | `cohort_tests_overall` |
| cohort comparison (pairs) | Cliff's delta, medians | Wilcoxon `p_raw`, Holm `p_holm` | `cohort_tests_pairwise` |
| missing share by cohort | OR (+ 95% CI) vs Healthy | GLM t-test `p`, `dispersion` | `missing_fraction_models` |
| onset vs continuation | `p_onset`, `p_continue`, `persistence_ratio` (+ bootstrap CI) | none (intervals only) | `transition` |
| week / year / weekday shares | pooled missing % (+ bootstrap CI) | none | `study_week_rates`, `calendar_year_rates`, `weekday_rates` |
| weekday effect | OR vs Sunday (Wald CI) | LRT `chisq`, `df`, `p` (weekday; weekday x cohort) | `weekday_odds_ratios`, `weekday_tests` |
| metric values by cohort / stratum | median, mean, SD, quartiles | none in tables (KW p in figures) | `metric_summary*` |
| gap onset (all specs) | `or` per SD below usual (+ CI), `estimate` = log-odds | Wald z = `estimate / se`, `p`; primary: `p_holm` | `onset_coefficients` (every fixed effect), `primary_test`; `onset_deciles` = observed gap rate by decile of the primary term |
| onset heterogeneity | `sigma2`, `mor` | - | `onset_models` |
| cohort difference in the primary effect | per-cohort OR | LRT `chisq`, `df`, `p` | `cohort_specific`, `cohort_interaction` |
| absolute risk | % risk at -1, 0, +1 SD (+ interval) | - | `risk_unadjusted` |
| sensitivity | OR per estimator | Wald `p`; LRT `p` for spline / random slope | `sensitivity`, `sensitivity_all_terms` |
| model B | change per doubling of gap length (SD units) | normal-approximation `p` | `model_B_post_gap_activity` (+ `model_B_curve`, `model_B_observed`, `model_B_all_terms` = every fixed effect) |
| model C | OR of a longer gap per SD below usual; predicted probabilities | Wald `p`; proportional-odds LRT `po_p` | `model_C_gap_length`, `model_C_predicted_probs`, `model_C_all_terms` = every coefficient |

## 9. Multiplicity, assumptions and what the results can and cannot say

* **Multiplicity**: Holm is applied (a) to the pairwise cohort tests within each stream x variable and (b) to the primary onset
  test across streams. Nothing else is adjusted: the Kruskal-Wallis tests across 3 variables x 3 streams, all exploratory onset specifications,
  the weekday tests, the sensitivity analyses and models B and C should be read as descriptive or exploratory.
* **Independence**: person-level tests treat people as independent (hence the overlap policy). Within a person, the mixed models use a random
  intercept; serial dependence beyond the gap-history terms, calendar and study-time terms is not modelled, and the bootstrap/clustered
  analyses are the robustness checks for that.
* **Informative missingness is not MNAR proof**: lower activity before gaps shows that missingness is **not completely at random**
  (it depends on observed prior values); MNAR concerns the *unobserved* activity during the gap, which the data cannot show.
  Models B and C probe plausibility only. All associations are non-causal.
* **Selection**: the upstream filter limits gap length to <= 14 days and keeps good wearers; estimates are not generalisable to people with long
  gaps or drop-out. Prostate is all male and almost entirely older; adjusting for sex within it is impossible (the term is dropped).
* **Heart-rate "coverage"** (`hr_cover`) is a wear measure; its association with gaps is a wear-pattern signal, not physiology.
* **Reading an odds ratio**: with a rare outcome (a gap starts on a small share of days) the odds ratio approximates a rate ratio, but it is
  never reported as a risk ratio; absolute risks come only from `risk_unadjusted`.
* **Nothing is subsampled automatically**; if a model cannot fit in memory it is recorded as failed. `sample_persons` is an explicit option that
  widens intervals but does not change the model.
