# Statistical methods of the missingness pipeline

This document describes what the code computes. Every statement was written from the code in `R/` (file names are given in
brackets). Nothing here is a result: no real-data result exists yet, so no estimate is quoted (the only numbers are from
simulated data with a known truth, section 13). Where a number is a convention from outside the code (for example a rule of
thumb for reading an effect size) it is marked **[convention]**.

**Reading order of the analyses** (the notebook cells follow the same order except that the stratified analysis, which needs the
primary models, is the last model cell):

| step | analysis | sections | notebook cell |
|---|---|---|---|
| 1 | descriptive statistics of missingness and of the metric values | 2, 4 | 11, 11b |
| 2 | temporal structure of missingness | 3 | 11 |
| 3 | data quality: wear consistency (characterization, classification rule fixed before any model) | 10.1-10.2 | 11c |
| 4 | primary models: gap onset (primary, model A, model D), gap-level models B and C, exploratory specifications | 5, 7 | 12-15, 17 |
| 5 | stratified analyses by wear consistency, with the formal comparison of the strata | 10.3-10.5 | 18 |
| 6 | sensitivity analyses (incl. full-observation-period standardization) | 6 | 16 |
| 7 | model diagnostics (multicollinearity, joint term tests, convergence) and integrated interpretation | 11, 9 | printed by 12-18 |

Every figure and table belongs to one of these steps; `figures/figure_index.tsv` lists each figure with its step and model, and
README.md lists the tables.

---------------------------------------------------------------------------------------------------

## 0. Design in one paragraph

Fitbit activity, sleep and heart-rate data of four cohorts (Insomnia, Diabetes, Healthy, Prostate) are reduced to
one row per person per day. A **day is "observed" if a row exists**. Missingness is studied only **between**
a person's first and last record (window (a)), so what is measured are *interruptions that were followed by
resumption*; leading missingness and permanent drop-out are invisible by construction. The analysis has four
layers: (1) descriptive person-level missingness and its temporal structure; (2) a mixed-effects logistic
model of **gap onset** given the previous day's values, whose predictors are standardized only with data observed
**before** each day; (3) two exploratory models that use **gap length** (which is only known after a gap ends and is
therefore never a predictor of onset); (4) the same models within strata of **wear consistency**, sensitivity analyses
(including a full-observation-period standardization) and collinearity diagnostics. No disease-label prediction
and no mixture model is part of the analysis.

**Sample sizes at each stage** are written, not quoted here: people per cohort after the overlap rule (`cohort_overlap.tsv`,
Cell 5) -> people / person-days per stream (`L1/manifest.tsv`, Cell 8) -> the upstream selection check (`selection_check`) ->
people and gaps per cohort (`cohort_missingness`) -> rows, people, events and people with events of every model
(`onset_models`: `n_rows`, `n_persons`, `n_events`, `n_persons_event`; `model_B_post_gap_activity` / `model_C_gap_length`: `n_gaps`,
`n_persons`) -> people per wear stratum (`wear_consistency_classification`, `wear_strata_quality`). A day needs at least
`cfg$standardization$baseline_min_n` (7) earlier observed days to have a within- or between-person z-score, so the first days of each person are not
in the models; the number of rows lost this way is the difference between the risk set and `n_rows`.

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

For each predictor x (transformed) of person i on day d, per stream, the day's value is split into **the person's usual level `m`
and the deviation from it**, and each part is put on a standard scale:

    x(i,d) = m(i,d) + sd_w * w(i,d)          m(i,d) = mu_b + sd_b * b(i,d)

* **Within-person (the primary predictor)** `w = (x - m) / sd_w`: how many within-person SDs day d was **above (+) or below (-)** the
  person's own usual. **Effect-size contrast: +1 SD higher than the person's usual** (`or` = odds ratio for a +1 SD higher
  day; OR < 1 means that a more active day is followed by lower odds of a gap).
* **Between-person** `b = (m - mu_b) / sd_b`: how many SDs the person's usual level is above (+) the average person. Contrast: +1 SD higher
  usual level. Keeping `w` and `b` as separate terms keeps the within- and between-person effects apart.
* (Sensitivity) **trailing** `tw = (x - mean of the previous 14 observed days) / sd_w`, requiring >= 7 of those days.

**Primary standardization = `pre_gap`** (`cfg$standardization$method`; function `predictor_table()`):

* **Usual `m` = the mean of the person's observed days strictly BEFORE day d** (an expanding history), needing at least
  `baseline_min_n` = 7 of them; a day without a baseline has no z-score and is not in the model. Only data observed before day d - hence
  before any gap that starts on day d+1 - enter the person's usual, so the z-scores carry no information from the person's future
  (checked by a unit test that changes a person's later values and finds the earlier deviations unchanged).
* **Reference distributions.** `sd_w` = root mean square of `x - m` over all days that have a baseline (the pooled within-person SD about the
  prior usual; **reference distribution of w: within-person deviations from the prior usual, pooled over people and days**); `mu_b`, `sd_b` = mean and SD
  of `m` over the same days (**reference distribution of b: the person-level prior usual, weighted by person-days**). These are three constants per
  predictor and stream, stored in `scales`; they rescale all people identically and cannot move one person's values relative to their own history.
  They are computed from all rows of the stream, i.e. they also use later data of OTHER people and later days of the same person, but only as
  a scale; they carry no person-specific future information and cannot create an association between a person's earlier value and a later outcome.
* **What this changes relative to the earlier specification** (now the `full_period` sensitivity): the usual is a noisier estimate for early
  days (7+ days of history, not the whole period), which adds measurement error to `w` and `b` early in follow-up and removes the look-ahead.
  The estimand is "association of a day's value relative to what the person had shown before" rather than "relative to the person's whole-period
  mean". Because a person's whole-period mean contains the days of the gap period itself, the full-period version can mix the outcome into the
  predictor (a person who goes quiet before a gap has a lower mean, which shrinks the apparent deviation of the days before the gap); this
  information leakage is the reason the pre-gap version is primary.

**Sensitivity standardization = `full_period`** (`standardization = "full_period"`; sensitivity analysis `full_period_standardization`,
`B_back2_adj_fullperiod`, `C_adj_fullperiod`): `m` = the person's mean over ALL observed days; `sd_w` = pooled within-person SD
`sqrt(sum (x - m)^2 / (N_obs - N_persons))`; `mu_b`, `sd_b` = mean and SD of the person means over people. The potential bias of this version is
(i) look-ahead in the baseline (the usual contains days after d and after the gap) and (ii) a within-person constraint: deviations of one person sum
to zero over the period, so a run of low days is partly offset by high days elsewhere. Its estimate is reported next to the primary one; a large
difference would indicate that the primary result depends on the baseline.

Both standardizations use one SD of **day-to-day variation of the person** as the unit of `w`; they differ in the centring (history vs whole
period) and in the SD (about the prior usual vs about the whole-period mean). Odds ratios are therefore comparable in direction and approximately in
size, not identical quantities.

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
* **Odds ratio**, as stored in `or`, `or_low`, `or_high`: **OR = exp(beta) for every term, in the normal direction** - for
  within-person (`w_`) terms "per +1 SD *higher* than the person's usual" (OR < 1: a higher-than-usual day is followed by lower
  odds that a gap starts the next day; OR > 1: by higher odds); for `b_` terms per +1 SD of the person-level mean; for the adjusters
  per +1 unit (or against the reference level of a factor). `estimate` is always the log-odds `beta`. These are odds ratios, **not** risk ratios, and are
  **conditional** on the person's random intercept.
* **Primary hypothesis**: `w_act_steps` in the joint adjusted k = 1 model, one test per stream, **Holm-adjusted across
  the streams** (`p_holm` in `primary_test`). This is the only pre-specified test. Everything else (model D, the other
  specifications, sensitivity, B and C) is secondary or exploratory and is **not** multiplicity-adjusted.
* **Heterogeneity between people**: `sigma2` (random-intercept variance, latent scale) and the **median odds ratio**
  `MOR = exp( sqrt(2 sigma2) x z_0.75 )` (typical OR between two random people with the same covariates).
* **Convergence and singularity** are reported, not hidden: `status` is `ok`, `not_converged`, `not_estimable` (fewer than 100 events
  or 30 people with events) or `failed`; a singular fit (person variance ~ 0) is flagged (`singular`, note) and counts as converged.
* **Collinearity**: see section 11 (generalized VIF of every term of every fitted model, `model_collinearity`; the older `vif` table, the
  VIF of the `w_` predictors from their correlation matrix, is kept). **Joint tests of each term** (`onset_term_tests`, `model_term_tests`): for a
  term with several coefficients (cohort, weekday, calendar year, the study-time spline) the Wald chi-square test of all its coefficients at
  once; the spline coefficients are never read one by one.
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

Same question (is a higher-than-usual prior day followed by lower odds of a gap?) under different assumptions; one row per analysis
in `sensitivity`, with `estimate, se, conf_low, conf_high, p, or, or_low, or_high`. The same fits also store the estimates of **every** term of
the model (predictors and adjustment covariates; not the intercept, and not the spline basis, which is tested jointly) in `sensitivity_all_terms`;
`sensitivity_forest_<stream>` draws the within-person predictors and `sensitivity_all_terms_other_/year_<stream>` all terms, primary vs each
sensitivity analysis (the spline and random-slope analyses are tests and are not in that table). Multicollinearity of each refit: `sensitivity_collinearity`:

| analysis | what changes |
|---|---|
| `main` | the primary model itself (pre-gap standardization) |
| `full_period_standardization` | the same model with the person's usual (and the between-person reference) taken from the **whole observation period** (`prepare_stream_bundle(..., standardization = "full_period")`); compare with `main` to see whether the result depends on the baseline definition and how much look-ahead matters |
| `trailing_baseline` | `tw_` (previous 14 observed days) instead of `w_`: no look-ahead in the baseline; the first observed days lack a baseline and drop out |
| `exclude_partial_wear` | drops rows whose previous day looks like partial wear (steps < 500 or sedentary >= 1380 min; sleep: minutes in bed < 60; heart rate: total minutes < 1200) |
| `first_gap_only` | only days up to and including each person's first gap onset (at most one event per person; no gap history) |
| `main_sleep_only` | sleep only: a night counts as observed only if a main sleep exists (naps alone do not count) |
| `conditional_logit` | `survival::clogit` with `strata(person_id)`, Efron approximation, persons with >= 1 event only; time-invariant terms (cohort, age, sex, `b_`) cancel. No random-effect assumption. (Before this revision the formula was written `survival::strata(person_id)`, which survival does not recognise as the stratification special: the person ids then entered as dummy covariates. This was corrected; earlier conditional-logit numbers should not be compared.) |
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
Outcome: the z-score of log(1 + steps) on the **first full day back** (`next observed date + 1`, primary), or on the first
day back (sensitivity; the watch is often put on mid-day), **against the person's usual THROUGH the last observed day before the gap** (pre-gap
standardization: `(x_post - m_through_prev_day) / sd_w`; 0 = the person's usual before the gap, + = more steps) - so the outcome is the change relative to
data observed before the gap. The pre-gap day itself enters as `w_pre_day_steps` (against the usual before that day). Gaps with fewer than
7 observed days before them have no baseline and are not in the models. The full-period sensitivity (`B_back2_adj_fullperiod`, `C_adj_fullperiod`)
uses whole-period z-scores for the pre-gap day and the days back. Four specifications: unadjusted; adjusted; adjusted + the pre-gap
deviation (change from before to after); first-day-back outcome adjusted. **Effect size**: `estimate` = change in the
outcome per **doubling** of the gap length, in within-person SD units, `se`, `conf_low/high = est +/- 1.96 se`, and `p` from a
**normal approximation** to the t statistic (`2 pnorm(-|est/se|)`; lme4 gives no p value; with 10^4-10^5 gaps the difference to a t reference is negligible).
Also reported: `n_gaps`, `n_persons`, outcome mean/SD, `sigma_person`, `sigma_resid`, `singular`, `converged`, `status`.
Figures: `gap_dose_response` draws the fitted straight line (interval from the fixed-effect covariance, covariates at their sample mean) at gap lengths
1, 2, 3, 5 and 13 days; `gapB_forest_<stream>` every coefficient of each specification.

### 7b. Observed activity after a gap, next to the model [`gap_models.R`, `plots_effects.R`]

`post_gap_activity_observed_<stream>` (table `model_B_observed_trajectory`) shows the data behind model B without fitting anything:
* **Time origin**: day 0 = the first observed day after the gap (the return day); day k = the calendar day k days after it, **k = 0 ... 7** (an 8-day window).
  A day without an activity record contributes nothing to that k, so the number of gaps and people **decreases with k**; both are in the table and
  in each panel title (day 0).
* **Value**: z-score of log(1 + steps) against the person's usual through the last day before the gap, in within-person SD units.
* **Groups**: gap length 1, 2, 3, 4-7, 8-14, 15+ days and all gaps. **Summary**: the mean of the gaps' values, with the person-clustered ratio-estimator
  interval (a person with several gaps is one cluster; `se` = half-width / 1.96; the figure follows `cfg$plot$interval_style`).
* **Model-based estimates are drawn apart** (red open diamond, day 1 only: the day model B describes): the unadjusted model B value at the group's median
  gap length. It is not the same quantity as the observed group mean (a group spans several lengths, the model is a straight line in log2 length), and
  the observed means are unadjusted; the comparison shows whether the model's line passes through the data, not a test.

**Model C - "does low activity the day before a gap predict a longer gap?"** Cumulative-logit (proportional-odds)
mixed model (`ordinal::clmm`, Laplace; sparse random-effect matrices) of the ordered length group (1, 2, 3, 4-7, 8-14, 15+)
on the prior-day within-person steps z-score (unadjusted, and adjusted for cohort, age, sex), with a person random
intercept: `logit P(Y <= j) = theta_j - beta w - ... - u_i`. It is run only if the **proportional-odds assumption** is
not rejected: a likelihood-ratio test of the model without random effect against the same model with a non-proportional (nominal) effect of `w`
(`po_p`, rejected if `po_p <= 0.05`; status `po_rejected`, no estimate; model D is then the analysis to report).
**Effect size**: `or = exp(beta)`, the odds of a **longer gap per +1 SD higher than the person's usual** (OR > 1: longer gaps after more active days), with Wald CI and z-test p.
If the random-intercept variance is on the boundary or the Hessian is unusable, the model **without** the random effect is used and
**all uncertainty comes from a person-level cluster bootstrap** of the full parameter vector (`ok_bootstrap_se`; replicate b is
seeded with `seed + b`, so results do not depend on the number of workers); the SE is then the bootstrap SD and the interval
`est +/- 1.96 se`. Predicted probabilities of each length group at w = -1, 0, +1 are **population-average** (averaged over a seeded
sample of 500 gaps' covariate profiles and over the random intercept), with 95% intervals as percentiles of 300 parameter draws
(multivariate normal from the Hessian, or the bootstrap draws in the fallback).

### 7c. MVPA before the gap as a covariate of models B and C [`gap_models.R`, `plots_effects.R`]

Models B and C were also fitted with **moderate-to-vigorous physical activity (MVPA)** as an additional covariate (units `B_back2_adj_mvpa`, `B_back2_adj_pre_mvpa`, `C_adj_mvpa`).
All earlier specifications, variable definitions and standardization are unchanged and are still reported; the MVPA fits are **added next to them**, not instead of them, so each
estimate can be read with and without MVPA.

* **Definition (unchanged variable).** MVPA = fairly active + very active minutes of the day (`mvpa_minutes`), transformed `log(1 + x)` (`cfg$predictors$act_mvpa`).
  Covariate `w_pre_day_mvpa` = that value on the **last observed day before the gap**, minus the mean of the person's earlier observed days (needs >= 7), divided by the
  pooled within-person SD of MVPA - exactly the pre-gap standardization of `w_pre_day_steps` (steps), on the same day. Unit: **per +1 SD higher MVPA than the person's usual**.
* **Alignment with the observation period (checked by tests).** The covariate is measured on the same day as the exposure of model C and the pre-gap covariate of model B, which is
  before the gap starts and therefore before the gap length (model C outcome) and before the days after the gap (model B outcome): no information from the gap or later enters
  it (changing MVPA on or after the return day leaves it unchanged; only earlier days move it). It is **not** the MVPA of the days after the gap: that is a component of the
  outcome of model B (steps and MVPA of one day move together), and adjusting for it would remove the effect being studied.
* **Why adjusting is appropriate, and what changes.** In model B, vigorous activity before the gap is a plausible common correlate of gap length and of the activity level after
  the gap (fitness, health); in model C, MVPA is correlated with the exposure (steps) and may carry part of its association with gap length. Adding it answers "the association of gap
  length (B) / of steps (C) **at the same pre-gap MVPA**". The steps coefficient of model C and of `B_back2_adj_pre_mvpa` changes meaning accordingly: steps beyond the intensity of the
  day. Both versions are reported; neither replaces the other. Models are conditional on a gap having occurred, so MVPA is not a mediator of gap onset here.
* **Collinearity.** Steps and MVPA of one day are correlated; the generalized VIF of both terms is in `gap_mvpa_effects` (`vif_equiv`; >= 5 is "high", section 11) and in
  `model_collinearity`. A high value widens the intervals of both terms; nothing is removed.
* **Missing MVPA.** Gaps whose MVPA covariate is missing (no MVPA value, or fewer than 7 earlier observed days with MVPA) are dropped from the `*_mvpa` fits and counted in
  `n_gaps_no_mvpa`; the MVPA model therefore uses slightly fewer gaps than its counterpart (compare `n_gaps`).
* **Model C specifics.** MVPA enters the cumulative-logit model, the proportional-odds check (only the steps effect is tested for proportionality, MVPA is assumed proportional) and the
  predicted probabilities (averaged over the MVPA distribution of the analysed gaps) as an ordinary covariate.
* **Limitations.** MVPA minutes are zero-inflated, so the log(1 + x) z-score is skewed; only the day before the gap is used (not the person's habitual MVPA level, which the
  onset model carries as `b_act_mvpa`); the MVPA fits are exploratory and not multiplicity-adjusted.

**Effect estimates and intervals.** `gap_mvpa_effects` lists, per model, stream and pair of fits, the MVPA coefficient and the main terms with and without MVPA with 95%
intervals: model B in SD units of the activity on the first full day back (per +1 SD MVPA / steps, per doubling of the gap), model C as odds ratios of a longer gap per +1 SD;
`results_summary.md` prints the MVPA lines. Figures: `gapB_mvpa_effect_<stream>`, `gapC_mvpa_effect_<stream>` (the other covariates stay in `gapB_forest_` / `gapC_forest_`).

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
| gap onset (all specs) | `or` per +1 SD higher than usual (+ CI), `estimate` = log-odds | Wald z = `estimate / se`, `p`; primary: `p_holm` | `onset_coefficients` (every fixed effect), `primary_test`; `onset_deciles` = observed gap rate by decile of the primary term |
| onset heterogeneity | `sigma2`, `mor` | - | `onset_models` |
| cohort difference in the primary effect | per-cohort OR | LRT `chisq`, `df`, `p` | `cohort_specific`, `cohort_interaction` |
| absolute risk | % risk at -1, 0, +1 SD (+ interval) | - | `risk_unadjusted` |
| sensitivity | OR per estimator | Wald `p`; LRT `p` for spline / random slope | `sensitivity`, `sensitivity_all_terms` |
| model B | change per doubling of gap length (SD units) | normal-approximation `p` | `model_B_post_gap_activity` (+ `model_B_curve`, `model_B_observed`, `model_B_all_terms` = every fixed effect) |
| model C | OR of a longer gap per +1 SD higher than usual; predicted probabilities | Wald `p`; proportional-odds LRT `po_p` | `model_C_gap_length`, `model_C_predicted_probs`, `model_C_all_terms` = every coefficient |
| MVPA before the gap (models B, C) | B: SD change per +1 SD MVPA; C: OR of a longer gap per +1 SD MVPA (with and without MVPA in the model) | Wald `p`, `vif_equiv` | `gap_mvpa_effects` |
| observed activity after gaps | mean z-score by gap length and day since return, clustered interval, n gaps / people | none | `model_B_observed_trajectory` |
| full-period standardization (sensitivity) | OR per +1 SD, same terms as the primary model | Wald `p` | `sensitivity` / `sensitivity_all_terms` (`analysis = full_period_standardization`), `model_B_post_gap_activity` / `model_C_gap_length` (`*_fullperiod`) |
| wear consistency | coverage, strata, composition (SMD) | none | `wear_consistency_classification`, `wear_consistency_descriptive`, `wear_consistency_composition`, `wear_consistency_by_person` |
| stratified models | OR (onset, C) / SD change (B) per stratum, every term | strata compared: LRT of (steps x stratum) (onset); Wald test of the difference of independent estimates (B, C) | `wear_strata_onset_coefficients`, `wear_strata_gapB`, `wear_strata_gapC` (+ `*_all_terms`), `wear_strata_interaction`, `wear_strata_comparison`, `wear_strata_quality` |
| multicollinearity | GVIF, GVIF^(1/(2 df)), VIF-equivalent, flag | none (diagnostic) | `model_collinearity` (all models), `onset_collinearity`, `sensitivity_collinearity`, `gap_collinearity`, `wear_strata_onset_collinearity` |
| joint tests of terms | chi-square of all coefficients of a term | Wald chi-square, `df`, `p` | `model_term_tests`, `onset_term_tests`, `gap_term_tests`, `wear_strata_onset_term_tests` |

## 9. Multiplicity, assumptions and what the results can and cannot say

* **Multiplicity**: Holm is applied (a) to the pairwise cohort tests within each stream x variable and (b) to the primary onset
  test across streams. Nothing else is adjusted: the Kruskal-Wallis tests across 3 variables x 3 streams, all exploratory onset specifications,
  the weekday tests, the sensitivity analyses and models B and C should be read as descriptive or exploratory.
* **Independence**: person-level tests treat people as independent (hence the overlap policy). Within a person, the mixed models use a random
  intercept; serial dependence beyond the gap-history terms, calendar and study-time terms is not modelled, and the bootstrap/clustered
  analyses are the robustness checks for that.
* **Informative missingness is not MNAR proof**: a different activity level before gaps (higher or lower than usual) shows that missingness is **not completely at random**
  (it depends on observed prior values); MNAR concerns the *unobserved* activity during the gap, which the data cannot show.
  Models B and C probe plausibility only. All associations are non-causal.
* **Selection**: the upstream filter limits gap length to <= 14 days and keeps good wearers; estimates are not generalisable to people with long
  gaps or drop-out. Prostate is all male and almost entirely older; adjusting for sex within it is impossible (the term is dropped).
* **Heart-rate "coverage"** (`hr_cover`) is a wear measure; its association with gaps is a wear-pattern signal, not physiology.
* **Reading an odds ratio**: with a rare outcome (a gap starts on a small share of days) the odds ratio approximates a rate ratio, but it is
  never reported as a risk ratio; absolute risks come only from `risk_unadjusted`.
* **Nothing is subsampled automatically**; if a model cannot fit in memory it is recorded as failed. `sample_persons` is an explicit option that
  widens intervals but does not change the model.

## 10. Stratified analysis by wear consistency [`wear_consistency.R`]

**Purpose.** Do the findings (prior-day steps and the odds of a gap; activity after a gap; gap length) hold in people who wore the device consistently
and in people who wore it inconsistently? The analysis comes last because it needs the primary models and uses them unchanged.

### 10.1 Measure and classification (fixed before any stratified model)

* **Wear consistency of a person in a stream = observed-day coverage** = `observed days / days between the person's first and last record in that stream`
  (= `1 - pct_missing / 100` of the missingness summaries; `wear_coverage_table()`). It is reproducible from the day tables alone and uses the same
  observation window as the gap definition. People were selected upstream on >= 70% coverage, so coverage lies in [0.70, 1]; "inconsistent" means
  "more interruptions than the typical participant", not "barely wore it".
* **Rule** (`cfg$wear`, written to `wear_consistency_classification` **before** any stratified model is fitted, and part of every stratified checkpoint key):
  `method = "median"` (default): **Consistent = coverage >= the stream's median** over the analysed people, otherwise Inconsistent (ties at the median
  are consistent); `method = "fixed"`: Consistent = coverage >= `cfg$wear$threshold`. The median split gives strata of similar size and is a property
  of the sample, not tuned to any result. The stratum is defined in the **stream that is modelled**.
* **What the stratification can and cannot show.** Coverage is built from the same observed days that define the outcome: a person with many
  gaps has lower coverage and contributes more events. The strata therefore differ by construction in how often gaps occur, and a difference in the
  stratified estimates is **effect modification by wear pattern**, not an effect of wear consistency, and the strata are not exchangeable
  (see the composition table). Conditioning on a function of the outcome can also change which kind of person is in which stratum; nothing here
  corrects for that.

### 10.2 Descriptive statistics (stage `wear_consistency`, Cell 11c; no model)

`wear_consistency_descriptive`: per stream and group (All / Consistent / Inconsistent): people, coverage median (Q1-Q3), follow-up and observed days,
gaps per person and per 100 observed days, % with any gap, longest gap, longest run of observed days, age, % female, CCI, cohort shares.
`wear_consistency_composition`: **standardized mean differences** (Consistent minus Inconsistent) in age, sex, CCI, follow-up length and each
cohort share (`SMD = (m1 - m2) / sqrt((v1 + v2) / 2)`; for binary variables with the proportion variance `p(1-p)`); `|SMD| > 0.1` is read as a relevant
imbalance **[convention]**. `wear_consistency_by_person`: the person-level table behind the figures. Figures: `wear_consistency_distribution` (ECDF of
coverage by cohort with the cut-off), `wear_consistency_composition`.

### 10.3 Models within each stratum (stage `wear_strata`, Cell 18)

The existing models are refitted **unchanged** (same formulas, predictors, outcomes, standardization) on the people of one stratum:
* the **primary onset model** (joint, adjusted, k = 1; `cfg$wear$onset_roles` = `"primary"`, add `"secondary_threshold"` for model D);
* **model B** (adjusted; first full day back) and **model C** (adjusted, ordinal);
the stratum is defined per stream, and `cohort`, `age10`, `sex2` are dropped inside a stratum where they are constant. Strata with fewer than
`cfg$wear$min_persons_stratum` (100) people are described but not modelled. Not stratified: the exploratory grid, model D by default and the other sensitivity analyses.

### 10.4 Comparing the strata (the formal question)

* **Onset model**: a likelihood-ratio test of `w_act_steps x stratum` (the stratum main effect is in both models; `df = 1`), fitted on all people with the
  common adjusters (`wear_strata_interaction`). This is the test of "the same within-person association in both strata".
* **Models B and C**: the strata contain different people, so the two stratum-specific estimates are independent and
  `z = (b_inconsistent - b_consistent) / sqrt(se_c^2 + se_i^2)` tests the difference (`compare_strata()`). For the log-odds models (onset, C) the exponentiated
  difference is the **ratio of odds ratios, Inconsistent / Consistent** (`ratio_of_or` with 95% CI); for B it is a difference in SD units per doubling.
* `wear_strata_comparison` lists, per model and stream, both strata's effect (OR with 95% CI, or SD change), the difference or ratio with its interval, `p_wald`
  and `p_difference` (the LRT p for the onset model, otherwise the Wald p). **One significant and one non-significant stratum are never read as a difference.**
  The strata are not multiplicity-adjusted.

### 10.5 Sparse data, convergence and composition (`wear_strata_quality`)

Per fitted stratum: people, events (onset) or gaps (B, C), coefficients, **events per coefficient** (flag `sparse` if below `cfg$wear$epv_min` = 10 **[convention]**),
`status` (`ok`, `ok_bootstrap_se`, `not_converged`, `not_estimable`, `po_rejected`), `singular`, and `usable` (= estimable, converged and not sparse). The model
information (`wear_strata_onset_models`), collinearity and term tests of each stratum are saved as well. Differences in composition between the strata
(10.2) can explain a difference between the strata estimates that has nothing to do with wear; the adjusted models adjust only for cohort, age, sex and the
usual covariates.

## 10b. Fitbit device version: missing and composition by cohort [`device.R`]

**Data.** The CDR `device` table (exported by Cell 6c) has one row per person x device x record date with `device_version` (e.g. "Charge 5", "Versa 4", "MobileTrack") and `device_type` (TRACKER or SCALE). It is not one row per day, some rows have no date, and it also contains Fitbit scales; a scale is not a wearable.

**Definitions (fixed before looking at any result).**
* *Watch version known*: at least one row that is neither a scale nor a device that is not worn on the wrist, with a non-empty `device_version`. "Watch" means worn on the wrist: `cfg$device$forms` classifies a version as smartwatch (Versa, Sense, Ionic, Blaze, Surge, Pixel Watch), band (Charge, Alta, Inspire, Flex, Force, Luxe, Ace), clip-on (One, Zip, Ultra, Classic) or phone (MobileTrack); `cfg$device$wrist_forms` (default smartwatch + band) says which count. This classification is from product knowledge, not from the data (the export only distinguishes TRACKER and SCALE): check it, and change `wrist_forms` to "smartwatch" for smartwatches only. Everybody else is **missing**: no device record, scale only, phone or clip-on only, or rows without a version (counted separately in `device_missing_by_cohort`).
* *Main device*: the wearable version with the most records per person (ties: the most recent record, then alphabetical order). *Family*: the first matching pattern of `cfg$device$families` (Charge, Versa, Sense, Inspire, Alta, Luxe, MobileTrack), any other wearable is "Other tracker". The number of different wearable versions per person is kept (`device_switching_by_cohort`).
* *Device on a given day* (`device_at_dates`, for later use; not used by any model yet): the wearable record of the same person on or before that day (as-of join); days before the first record are NA unless `backfill = TRUE`.

**Statistics (stage `device`, Cell 11d; person level; seed `cfg$seed`).**
1. *Missing share per cohort*: proportion with the Wilson score 95% interval (`device_missing_by_cohort`).
2. *Do cohorts differ?* Pearson chi-square test of the cohort x (missing / known) table, with the bias-corrected Cramer's V (Bergsma) and its percentile bootstrap 95% interval (cohort sizes held fixed, each cohort resampled as a multinomial; `cfg$device$boot_B`). Each cohort against the reference cohort: Fisher exact odds ratio with its exact interval, Newcombe hybrid-score interval of the risk difference, Holm adjustment over the comparisons (`device_missing_vs_reference`).
3. *Adjusted for age and sex*: logistic regression of "missing" on cohort (reference Healthy), age (per 10 years) and sex, Wald intervals, joint likelihood-ratio test of cohort (`device_missing_model`). A cohort with no (or only) missing people cannot be estimated by maximum likelihood (separation): it is left out of the model and reported as "not estimable"; its exact Fisher result is in step 2.
4. *Composition*: share of each main-device family per cohort with Wilson intervals (`device_composition_by_cohort`). The cohort x category chi-square of independence is reported with the asymptotic p value AND a Monte-Carlo p value (`cfg$device$mc_B` replicates), because many cells are sparse (`min_expected`, `pct_cells_expected_lt5` are in `device_tests`); it is done with the missing category included and among people with a known device only; V with bootstrap interval. Adjusted standardized residuals (Haberman) per cell with Holm adjustment over all cells show where the cohorts differ (`device_composition_residuals`). Version-level counts: `device_version_by_cohort`.

**Figures.** `device_version_missing_and_composition`: one 100% bar per cohort, grey segment = missing, coloured segments = watch family, the tests in the subtitle (Wilson intervals are in the tables, not on the bars). `device_cohort_mix_by_version`: the same cohort x category table read the other way, one 100% bar per watch version (and Missing) coloured by cohort, with the cohort mix of all people as the top reference bar (`device_cohort_mix`); the test is the same independence test.

**What this does and does not show.** It describes who has which device in this sample; the p values test association between cohort and device information in these people, not causes. People with a device record but no version or with several devices are classified by the rule above, which is a choice. A significant test with a small V means a real but small difference; read V next to p. A missing device version is likely to be informative about data quality (the person's device history was not captured), so before adding the device as a covariate to a model, check how the missing group differs.

## 10c. Person-level % of days missing against every model variable [`outcome_variables.R`]

**Question.** How does the person-level outcome "% of days missing" (one value per person and stream) relate to each variable of the models, one variable at a time (stage `variables`, Cell 11e; descriptive, no model of the gap outcome)?

**Variables.** Categorical: cohort, sex, watch version (device family; "Missing" = no watch version, so nobody is dropped). Continuous: age at the first record and the between-person level of every model predictor of the stream (the person's mean over observed days; geometric mean for the log-transformed predictors: steps, MVPA, time in bed, elevated heart-rate minutes). The day-level variables of the models (weekday, study time, calendar year) are shown at the day level in the temporal figures; previous gap length, "any earlier gap" and days since return are built from the outcome itself and are not compared with it.

**Tests** (every test uses one value per person, so people with long follow-up do not dominate; levels with fewer than `cfg$outcome_vars$min_n` = 5 people are drawn but not tested).
* *Categorical*: Kruskal-Wallis test (three or more levels; effect size epsilon-squared = H/(n-1)) or Wilcoxon rank-sum test (two levels; Cliff's delta), with a percentile bootstrap 95% interval of the effect size (people resampled, `cfg$outcome_vars$boot_B`). Each level against the reference level (cohort: Healthy; sex: Female; watch: the most common watch family): Wilcoxon test with the Hodges-Lehmann difference in percentage points and its 95% interval, Cliff's delta, Holm adjustment within the variable (`outcome_variable_contrasts`).
* *Continuous*: Spearman rank correlation with the Bonett-Wright 95% interval; the figure also shows the mean % missing in quantile bins with its 95% interval so that non-monotone shapes are visible.
* *Multiplicity*: the overall p values of all variables of one stream are Holm-adjusted together (`p_holm` in `outcome_variable_tests`).
* *Adjusted association* (`outcome_variable_models`): quasi-binomial model cbind(missing days, observed days) ~ variable + the other person-level covariates (cohort, sex, age; the variable itself is left out of its own adjustment set), giving the odds ratio that a day is missing: each level against the reference, or per +1 SD of a continuous variable, with 95% intervals from the t distribution of the quasi-binomial fit and the F test of the variable. The unadjusted fit is reported next to it. The between-person levels and the age are correlated with cohort, which is why the adjusted estimate is given.

**Figures.** `outcome_by_variable_categorical` (columns: cohort, sex, watch version; rows: streams; violin + box, n per level, overall test and effect size in the corner) and `outcome_by_variable_continuous_<stream>` (one panel per variable: people as points coloured by cohort, bin means with 95% intervals, Spearman rho and the adjusted OR per SD in the corner; points outside the 0.5-99.5% range of x or above the 99th percentile of y are not drawn, the bins and tests use everybody).

**What this does and does not show.** Associations in this sample, one variable at a time; the adjusted odds ratios take cohort, sex and age into account but not the other variables. The between-person levels are computed from the observed days, which are the days that are not missing, so a relation with the missing share is partly mechanical (fewer observed days give a noisier mean). A small effect size with a very small p value reflects the large sample.

## 11. Multicollinearity diagnostics and joint term tests [`collinearity.R`]

**Method (all models with several fixed effects: onset grid, model D, sensitivity refits, models B and C, every stratified fit).** For each fit the
covariance matrix `V` of the fixed-effect estimates (intercept, and for the ordinal model the thresholds, removed) is converted to the correlation matrix
`R = cov2cor(V)`. This is the approach of `car::vif()` / `performance::check_collinearity()` and is valid for linear, logistic and ordinal mixed models,
because it measures the inflation of the uncertainty of an estimate through its correlation with the other estimates, **with the model's weights**, not only
the correlation of the raw columns. Terms that need several coefficients (cohort, weekday, calendar year, the study-time spline) are **one term**:

    GVIF(term) = det(R_term) * det(R_other terms) / det(R) ,  df = number of coefficients of the term (Fox & Monette 1992)
    GVIF^(1/(2 df))  is on the scale of a standard error ratio;   VIF-equivalent = GVIF^(1/df) = (GVIF^(1/(2 df)))^2  (= the ordinary VIF if df = 1)

**Thresholds [convention]**: VIF-equivalent >= 5 "high", >= 10 "severe" (`flag`); they are rules of thumb, not tests. **Interpretation**: a VIF-equivalent of v inflates
the standard error of that term by `sqrt(v)` and makes its estimate unstable against small changes of the data; it does **not** bias the coefficients or the predictions, and
the odds ratios of the **other** terms are unaffected except through their own VIF. **Nothing is removed or re-specified because of a high VIF**; flagged terms
are listed in `results_summary.md` and the figure `collinearity_gvif_<stream>` (one panel per model: primary onset, full-period standardization, models B and C,
the two wear strata). **Within- vs between-person terms**: `w` is a deviation from the person's usual and `b` is the usual itself, so they are nearly uncorrelated by
construction; high VIFs, if any, are expected among the **between-person** terms (people who walk more also have more MVPA) and among the gap-history terms
(`prev_gap_len_l`, `any_prior_gap`, `days_since_return_l`). **Interactions**: the interaction models (`w x cohort`, `w x stratum`) are used only for likelihood-ratio
tests; their product terms are collinear with their main effects by construction, so GVIF is reported for the main-effects model (the stratified fits are
main-effects fits). The conditional logit (only time-varying terms, no between-person terms) and the pooled cluster-robust fit (the same design as the primary fit)
are not diagnosed separately.

**Joint Wald tests** (`onset_term_tests`, `gap_term_tests`, `model_term_tests`): `chisq = b' V^-1 b` on `df` = number of coefficients of the term. They are the only
correct reading of the study-time spline (never interpret `ns(...)1`, `ns(...)2`, `ns(...)3` separately), of calendar year and of weekday; a within-person predictor
with one coefficient gives `chisq = (b / se)^2`, the squared Wald z.

## 12. Effect-size figures: what is shown and how to read it

* **Logistic models (onset, sensitivity, stratified onset, model D)** and the **ordinal model C**: odds ratios with confidence intervals on a log axis, reference 1.
  Model C is a cumulative-logit model; its odds ratio is "odds of a longer gap (higher length group) per +1 SD", the same at every cut-point (proportional odds);
  predicted probabilities of each length group (`gap_length_probs`) are the more intuitive summary. **Model B** is linear: estimates are changes in activity in SD
  units, reference 0. Model labels and units are on every axis / subtitle.
* **Interval**: dot +/- 1 SD of the estimate (its standard error; `cfg$plot$interval_style = "ci95"` gives the 95% interval). A bar that does not cross 1 (0)
  is **not** the test: p values and, for strata, the formal comparison are in the tables.
* **All terms**: `onset_forest_*`, `sensitivity_all_terms_*`, `gap_definition_all_terms_*`, `wear_strata_*_all_terms_*` show every coefficient of the model
  (predictors, cohort, age, sex, weekday, gap history; calendar years in the `_year_` figure; reference levels Healthy, Female, Sunday, the calendar year with the most gap starts; sparse years are pooled). The study-time
  spline and the intercept are not drawn (joint tests instead). A term whose interval cannot be estimated (standard error of the log odds ratio above 10) is left out of the figure with a note in its subtitle; the tables keep it. A term whose interval cannot be estimated (for example a calendar year with almost no gaps, standard error of the log odds ratio above 10) is left out of the figure with a note in its subtitle; the tables keep it. The first (dark) colour is always the primary analysis; the other colours are the sensitivity
  analysis, gap definition or stratum it is compared with. Only terms that exist in a model are drawn (the conditional logit has no cohort or between-person term).
* **Observed vs model-based** (`post_gap_activity_observed_*`): see 7b; the predicted curves and probabilities are population-average and model-based.

## 13. Validation: what was tested and what was not

Implemented and tested (`tests/testthat/`): the pre-gap baseline (`prior_mean`, `through_mean`: strictly earlier observations only, minimum count, NA handling); the
decomposition `x = m + sd_w w`, `m = mu_b + sd_b b` and the reference distributions; **no look-ahead** (changing a person's later values leaves the earlier deviations
unchanged for the pre-gap standardization and changes them for the full-period one; the same for the gap table: changes after the post-gap window leave a gap's z-scores
unchanged except for the common scale constant); that the full-period standardization reproduces the earlier specification; odds ratio = exp(coefficient) for every term;
GVIF against the textbook VIF / the determinant formula, flags and the joint Wald tests against direct matrix computation; classification by wear consistency (median, fixed,
ties, determinism), the standardized mean differences, the comparison of strata (independent-estimate test, ratio of odds ratios, the reference does not change the test), that
stratified data are a partition of the people, that the interaction test detects a planted difference between strata (p < 0.001) and not a common effect, resumable
stratified units whose checkpoints depend on the classification rule, and that every new figure builds with the intended layers, labels and sizes.

Regression check: the `full_period` standardization reproduces the previous version of `predictor_table()` exactly (all `w_`, `b_`, `tw_` columns and the scaling constants identical, maximum absolute difference 0 on simulated data).

Simulation check of the primary estimate (24 replicates of a simulated study of 300 people with a **known** effect, log-odds -0.5 per +1 within-person SD of log steps; script
`pipeline/validate_simulation.R`): see section 13.1. **Not validated**: any behaviour on the real data (distribution of coverage, size of the strata, real collinearity, real
effect of the baseline definition); the stratified models at full size (time and memory were measured for the primary models only, see README); the cut-off of the median split
(a design choice); the sensitivity of the conclusions to the 7-day minimum baseline.

### 13.1 Simulation results (not real data)

| standardization | mean log-OR (truth -0.5) | empirical SD of the estimate | mean model SE | 95% intervals covering the truth | events per fit |
|---|---|---|---|---|---|
| pre_gap (primary) | -0.500 | 0.068 | 0.070 | 24 / 24 | about 595 |
| full_period (sensitivity) | -0.516 | 0.064 | 0.066 | 24 / 24 | about 647 |

The model standard errors match the replicate-to-replicate SD (0.070 vs 0.068), so the uncertainty is not understated. The pre-gap estimate is unbiased in this
setting. The full-period estimate is 0.016 larger in magnitude, which with 24 replicates (standard error of the mean 0.013) cannot be separated from noise; it has about 9% more
events because it needs no minimum baseline. The simulated effect does not depend on the person's future, so this check shows that dropping the look-ahead costs no accuracy
here; it cannot show how much look-ahead matters in real data, which is what the `full_period_standardization` sensitivity analysis is for.
