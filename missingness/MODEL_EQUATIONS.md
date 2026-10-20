# Model equations, outcomes and predictors, by stream

Written from the code (`make_onset_grid`, `build_onset_formula`, `run_sensitivity_stage`, `fit_gap_dose_response`,
`fit_gap_length_ordinal`); the formulas below are the strings the code builds (checked by printing them).
A person is `i`, an observed day of that person in the stream is `d`.

## 0. Common definitions

**Observed day** (per stream): a row exists in that stream's day table (activity table; sleep table: a main sleep or a nap
counts, except in the `main_sleep_only` sensitivity where only `has_main_sleep` counts; heart-rate table: the four zone
rows collapsed to one row per day). **Gap**: >= 1 consecutive calendar days with no row between two observed days.

**Predictors** (day `d` itself, i.e. the day before the possible gap), each transformed and then split into
`w_` (within-person z-score) and `b_` (between-person z-score):

| name | source column | transform | stream |
|---|---|---|---|
| `act_steps` | `steps` | log(1 + x) | activity (and used by every stream) |
| `act_sedentary` | `sedentary_minutes` | identity | activity (and used by every stream) |
| `act_mvpa` | `mvpa_minutes` | log(1 + x) | activity (and used by every stream) |
| `sleep_inbed` | `minute_in_bed` | log(1 + x) | sleep |
| `sleep_eff` | `sleep_efficiency` | identity | sleep |
| `hr_cover` | `total_hr_min` (minutes of heart rate recorded: **wear / coverage**, not behaviour) | identity | hr |
| `hr_elev` | `elevated_min` (fat-burn + cardio + peak minutes) | log(1 + x) | hr |

**Standardization (primary = `pre_gap`, `cfg$standardization`).** For the person's day `d` of predictor `p` (transformed value `x`):

    m(d)   = mean of the person's observed x on days BEFORE d          (>= 7 such days, else the day has no z-score)
    w_p(d) = (x(d) - m(d)) / sd_w          sd_w = sqrt( mean over all days with m of (x - m)^2 )       (within-person, pooled)
    b_p(d) = (m(d) - mu_b) / sd_b          mu_b, sd_b = mean, SD of m over the same days                (between-person level)
    x(d)   = m(d) + sd_w * w_p(d)          m(d) = mu_b + sd_b * b_p(d)

`w_p` is in within-person SD units: **+1 = one SD above the person's own usual** (the effect-size contrast of every odds ratio is +1 SD HIGHER, never below).
`sd_w`, `mu_b`, `sd_b` are constants per predictor and stream (they only rescale). `tw_p = (x - mean of the previous 14 observed days) / sd_w` (needs >= 7 of them; sensitivity).
**Sensitivity `full_period`:** `m = mean over ALL the person's observed days` (look-ahead), `sd_w = sqrt( sum (x - m)^2 / (N_obs - N_persons) )`, `mu_b`, `sd_b` = mean and SD of the person means.
Odds ratios of every term are `OR = exp(beta)`; for `w_` terms they are **per +1 SD higher** than the person's usual (OR < 1: higher-than-usual days are followed by fewer gaps).

**Adjusters** (adjusted models; `cohort` has **Healthy as the reference level** (`cfg$reference_cohort`), so each cohort coefficient compares that cohort with Healthy; a term that is constant in the data, e.g. `sex2` in an all-male cohort, is dropped automatically):
`cohort`, `age10 = age at first record / 10` (per 10 years, not centred), `sex2` (Female/Male), `dow_t` (weekday of day d+1; reference Monday), `ns(study_day, 3)` (natural spline, 3 df, of days since the
first record), `device` (the person's main Fitbit device FAMILY, time-invariant: reference = the largest family that is not Missing, "Missing" = no usable watch device, small wrist families pooled as "Other wrist-worn", clip-on/phone its own level), `cci` (Charlson score per +1 point; people without a score get the median and `cci_missing = 1`), `prev_gap_len_l = log2(1 + length of the most recent completed gap, 0 if none)`,
`any_prior_gap`, `days_since_return_l = log2(1 + days since the last return from a gap, not capped, 0 if none)`. Calendar year is not a covariate. CCI is strongly tied to cohort: its GVIF is reported and the sensitivity refits `cci_removed`, `cci_as_groups` (reference = most common group, "Unknown" level) and `cci_within_cohort` (score minus the cohort mean) are run.

## 1. Gap-onset model (primary, model D and the exploratory grid)

**Rows**: every observed day `d` of a person in the stream except the person's last, with all predictors of the stream non-missing ("common risk set":
the stream's own predictors AND the three activity predictors on the same date; for adjusted models also `age10`, `sex2`, `cohort`).
**Outcome** `y_ge{k}(i,d) = 1` if a gap of at least `k` days starts on day d+1 (next observed date - d - 1 >= k), else 0.
Primary `k = 1`; model D `k = 2, 3, 7`.

**Estimation**: `lme4::glmer`, binomial, logit, `(1 | person_id)`, `nAGQ = 0`, bobyqa. Wald z tests; `OR` as above.

### Primary model (joint predictors, adjusted, k = 1) - the one pre-specified test: `w_act_steps`, Holm across streams

**Activity stream** (outcome: an *activity* gap starts tomorrow)

    y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa
          + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + device + cci + cci_missing
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Sleep stream** (outcome: a *sleep* gap starts tomorrow; predictors: yesterday's sleep and activity)

    y_ge1 ~ w_sleep_inbed + w_sleep_eff + w_act_steps + w_act_sedentary + w_act_mvpa
          + b_sleep_inbed + b_sleep_eff + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + device + cci + cci_missing
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Heart-rate stream** (outcome: a *heart-rate* gap starts tomorrow; predictors: yesterday's heart-rate coverage / elevated minutes and activity)

    y_ge1 ~ w_hr_cover + w_hr_elev + w_act_steps + w_act_sedentary + w_act_mvpa
          + b_hr_cover + b_hr_elev + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + device + cci + cci_missing
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

(`sex2` is dropped automatically wherever it is constant, e.g. in the cohort-specific fit of Prostate; `cohort` is dropped in cohort-specific fits.)

### Model D (secondary): the same three formulas with the outcome `y_ge2`, `y_ge3`, `y_ge7` (a gap of >= 2, 3, 7 days starts tomorrow).

### Model A: not a separate fit - the gap-history terms `prev_gap_len_l + any_prior_gap + days_since_return_l` inside the adjusted formulas above
(they use only gaps that ended on or before day d; gap length is never a predictor of the gap that is being predicted).

### Exploratory grid (all k = 1; each unadjusted = only the `w_` and `b_` terms + `(1 | person_id)`, and adjusted = those + the adjusters above)

| stream | sets of predictors |
|---|---|
| activity | `uni_act_steps`, `uni_act_sedentary`, `uni_act_mvpa` (one predictor each), `joint` (= the three activity predictors) |
| sleep | `own_uni_sleep_inbed`, `own_uni_sleep_eff`, `own_multi` (both sleep predictors), `act_uni_act_steps`, `act_uni_act_sedentary`, `act_uni_act_mvpa`, `act_multi` (the three activity predictors), `joint` (all five) |
| hr | the same pattern with `hr_cover`, `hr_elev` instead of the sleep predictors |

Each set enters as `w_p + b_p` for every predictor `p` of the set. Example (unadjusted joint, activity):
`y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa + b_act_steps + b_act_sedentary + b_act_mvpa + (1 | person_id)`.

## 2. Sensitivity analyses of the primary effect (term `w_act_steps`; all use the primary spec of the stream)

| analysis | outcome | model / change (shown for the activity stream; sleep and hr add their own predictors exactly as in section 1) |
|---|---|---|
| `main` | `y_ge1` | the primary model above (shares its fit; pre-gap standardization) |
| `full_period_standardization` | `y_ge1` | the primary formula with `w_`, `b_` computed with `standardization = "full_period"` (whole-period usual; same rows except that no minimum baseline is needed); term reported: `w_act_steps` |
| `trailing_baseline` | `y_ge1` | `w_` replaced by `tw_`: `y_ge1 ~ tw_act_steps + tw_act_sedentary + tw_act_mvpa + b_act_steps + b_act_sedentary + b_act_mvpa + <adjusters> + (1 | person_id)`; term reported: `tw_act_steps`; rows without a 14-day baseline drop out |
| `exclude_partial_wear` | `y_ge1` | primary formula, rows dropped whose day d is flagged partial wear (activity: `steps < 500 or sedentary_minutes >= 1380`; sleep: `minute_in_bed < 60`; hr: `total_hr_min < 1200`; sleep/hr rows are dropped if the stream's own day OR the activity day is flagged) |
| `first_gap_only` | `y_ge1` | primary formula on days up to and including each person's first gap onset (at most one event per person); the three gap-history terms are constant there and are dropped |
| `main_sleep_only` (sleep only) | `y_ge1` | primary formula, but a night is "observed" only if it has a main sleep (`has_main_sleep`); gaps and risk set are rebuilt on that definition |
| `conditional_logit` | `y_ge1` | `survival::clogit`, Efron: `y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa + dow_t + ns(study_day, 3) + prev_gap_len_l + any_prior_gap + days_since_return_l + strata(person_id)` (`strata` written without the `survival::` prefix so that survival recognises it) (sleep/hr: with their own `w_` terms too); only persons with >= 1 event; no `b_`, cohort, age, sex or year terms (they cancel or are time-invariant); `dow_t` kept only if it varies |
| `cluster_robust_logit` | `y_ge1` | the primary formula without the random intercept (pooled `glm`, binomial), person-clustered sandwich SEs (`vcovCL`) |
| `spline_vs_linear` | `y_ge1` | primary formula with `w_act_steps` replaced by `ns(w_act_steps, 3)`; likelihood-ratio test against the linear model (2 df) |
| `random_slope` | `y_ge1` | primary formula `+ (0 + w_act_steps | person_id)` (uncorrelated random slope); LRT against the primary model, p halved (boundary, 1 df) |

For the sleep and heart-rate streams the tested term is still `w_act_steps` (yesterday's activity), with the stream's own `w_` terms in the model as above.

## 3. Gap-length models B and C (exploratory; one row per gap; same equations in every stream)

The **stream** only decides which gaps are used (gaps of the activity / sleep / heart-rate series); the **activity** values come from the activity table in all three:
`pre = w_act_steps` on the last observed day before the gap (against the usual before that day, as in the onset model);
`back1`, `back2` = steps (log1p) on the first observed day after the gap (`next observed date`) and on the day after that (`next observed date + 1`, which must itself be an
observed activity day), **against the person's usual THROUGH the last observed day before the gap**: `back = (x_post - m_through_prev_day) / sd_w` (data observed before the gap only;
0 = usual, + = more steps). Gaps with fewer than 7 earlier observed days have no usual. With `standardization = "full_period"` (units `B_back2_adj_fullperiod`, `C_adj_fullperiod`) all of
these are whole-period z-scores. `post_k0 ... post_k7` = the same z-score on the return day and the 7 days after it (observed trajectory, `model_B_observed_trajectory`).
`log2_gap = log2(gap length)`.

**Model B** (`lme4::lmer`, REML, Gaussian; effect = change per doubling of the gap, in within-person SD units; normal-approximation p)

    B_back2_unadj    : back2 ~ log2_gap                                   + (1 | person_id)
    B_back2_adj      : back2 ~ log2_gap + cohort + age10 + sex2           + (1 | person_id)
    B_back2_adj_pre  : back2 ~ log2_gap + pre + cohort + age10 + sex2     + (1 | person_id)
    B_back1_adj      : back1 ~ log2_gap + cohort + age10 + sex2           + (1 | person_id)
    B_back2_adj_fullperiod : B_back2_adj with the whole-period standardization (sensitivity)
    B_back2_adj_mvpa     : back2 ~ log2_gap + pre_mvpa + cohort + age10 + sex2          + (1 | person_id)
    B_back2_adj_pre_mvpa : back2 ~ log2_gap + pre + pre_mvpa + cohort + age10 + sex2   + (1 | person_id)

**Model C** (`ordinal::clmm`, cumulative logit; outcome = ordered gap-length group `1 | 2 | 3 | 4-7 | 8-14 | 15+`, groups with < 20 gaps merged into the next shorter one)

    logit P(length group <= j) = theta_j - ( beta * pre [+ gamma1 * cohort + gamma2 * age10 + gamma3 * sex2] ) - u_i ,   u_i ~ N(0, sigma_u^2)
    C_unadj : only beta * pre        C_adj : with the covariates        C_adj_fullperiod : C_adj with the whole-period standardization (sensitivity)
    C_adj_mvpa : logit P(length group <= j) = theta_j - ( beta * pre + delta * pre_mvpa + gamma1 * cohort + gamma2 * age10 + gamma3 * sex2 ) - u_i

`pre_mvpa = w_pre_day_mvpa` = log(1 + MVPA minutes) on the last observed day before the gap (the same day as `pre`), minus the mean of the person's EARLIER observed days (>= 7), divided by the
pooled within-person SD of MVPA: **per +1 SD higher MVPA than the person's usual**. Gaps without it are dropped from the `*_mvpa` fits. `delta` and the MVPA coefficient of B are reported
with 95% intervals in `gap_mvpa_effects` (B: SD of activity; C: odds ratio `exp(delta)` of a longer gap); the specifications without MVPA are unchanged.

Run only if proportional odds is not rejected (LRT of the same model without random effect against one with a non-proportional effect of `pre`, p > 0.05). Effect: `OR = exp(beta)` for a
longer gap per +1 SD higher than usual activity. If the random-intercept variance is on the boundary: the same cumulative-logit model without `u_i` and person-clustered bootstrap SEs.
(`cohort`, `age10`, `sex2` are dropped where constant.)

## 4. The same models within strata of wear consistency (stage `wear_strata`)

**Stratum**: `coverage = observed days / (last - first + 1)` in the modelled stream; `Consistent` if `coverage >= median(coverage)` (default; `cfg$wear`), else `Inconsistent`.
**Models** (formulas as above, fitted separately on the people of one stratum; `cohort`, `age10`, `sex2` dropped where constant):

    onset  : y_ge1 ~ <primary joint adjusted formula>                       within stratum s in {Consistent, Inconsistent}
    model B: back2_dev ~ log2_gap + cohort + age10 + sex2 + (1 | person_id) within stratum s
    model C: ordinal cumulative logit as above, C_adj                        within stratum s
    test of the difference, onset (pooled people):  y_ge1 ~ <primary formula> + stratum + w_act_steps:stratum + (1 | person_id)  vs the same without the interaction (LRT, 1 df)
    test of the difference, models B and C: z = (b_Inconsistent - b_Consistent) / sqrt(se_C^2 + se_I^2)   (independent estimates; strata hold different people)

Effect sizes: onset and C: `OR = exp(beta)` per +1 SD higher prior-day steps and the **ratio of odds ratios** `exp(b_I - b_C)`; B: change in activity (SD) per doubling of the gap and its difference.

## 5. Diagnostics

    V = covariance of the fixed effects without intercept / thresholds,  R = cov2cor(V)
    GVIF(term) = det(R_term) det(R_other) / det(R),  df = coefficients of the term,  VIF-equivalent = GVIF^(1/df)      (>= 5 high, >= 10 severe)
    joint Wald test of a term:  chisq = b_term' V_term^-1 b_term  on df = coefficients of the term
