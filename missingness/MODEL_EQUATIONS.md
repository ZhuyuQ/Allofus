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

`w_p = (x - mean_i(x)) / sd_w`, `sd_w = sqrt( sum (x - mean_i)^2 / (N_obs - N_persons) )` (pooled within-person SD, per predictor and stream);
`b_p = (mean_i(x) - mean of person means) / SD of person means`; `tw_p = (x - mean of the previous 14 observed days) / sd_w` (needs >= 7 of them; sensitivity only).
Odds ratios of `w_` terms are reported **per 1 SD below** the person's usual: `OR = exp(-beta)`.

**Adjusters** (adjusted models; `cohort` has **Healthy as the reference level** (`cfg$reference_cohort`), so each cohort coefficient compares that cohort with Healthy; a term that is constant in the data, e.g. `sex2` in an all-male cohort, is dropped automatically):
`cohort`, `age10 = (age at first record - 60)/10`, `sex2` (Female/Male), `dow_t` (weekday of day d+1), `ns(study_day, 3)` (natural spline, 3 df, of days since the
first record), `factor(year_t)` (calendar year of day d+1), `prev_gap_len_l = log(1 + length of the most recent completed gap, 0 if none)`,
`any_prior_gap`, `days_since_return_l = log(1 + days since the last return from a gap, capped at 60, 0 if none)`.

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
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + factor(year_t)
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Sleep stream** (outcome: a *sleep* gap starts tomorrow; predictors: yesterday's sleep and activity)

    y_ge1 ~ w_sleep_inbed + w_sleep_eff + w_act_steps + w_act_sedentary + w_act_mvpa
          + b_sleep_inbed + b_sleep_eff + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + factor(year_t)
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Heart-rate stream** (outcome: a *heart-rate* gap starts tomorrow; predictors: yesterday's heart-rate coverage / elevated minutes and activity)

    y_ge1 ~ w_hr_cover + w_hr_elev + w_act_steps + w_act_sedentary + w_act_mvpa
          + b_hr_cover + b_hr_elev + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + factor(year_t)
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
| `main` | `y_ge1` | the primary model above (shares its fit) |
| `trailing_baseline` | `y_ge1` | `w_` replaced by `tw_`: `y_ge1 ~ tw_act_steps + tw_act_sedentary + tw_act_mvpa + b_act_steps + b_act_sedentary + b_act_mvpa + <adjusters> + (1 | person_id)`; term reported: `tw_act_steps`; rows without a 14-day baseline drop out |
| `exclude_partial_wear` | `y_ge1` | primary formula, rows dropped whose day d is flagged partial wear (activity: `steps < 500 or sedentary_minutes >= 1380`; sleep: `minute_in_bed < 60`; hr: `total_hr_min < 1200`; sleep/hr rows are dropped if the stream's own day OR the activity day is flagged) |
| `first_gap_only` | `y_ge1` | primary formula on days up to and including each person's first gap onset (at most one event per person); the three gap-history terms are constant there and are dropped |
| `main_sleep_only` (sleep only) | `y_ge1` | primary formula, but a night is "observed" only if it has a main sleep (`has_main_sleep`); gaps and risk set are rebuilt on that definition |
| `conditional_logit` | `y_ge1` | `survival::clogit`, Efron: `y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa + dow_t + ns(study_day, 3) + prev_gap_len_l + any_prior_gap + days_since_return_l + strata(person_id)` (sleep/hr: with their own `w_` terms too); only persons with >= 1 event; no `b_`, cohort, age, sex or year terms (they cancel or are time-invariant); `dow_t` kept only if it varies |
| `cluster_robust_logit` | `y_ge1` | the primary formula without the random intercept (pooled `glm`, binomial), person-clustered sandwich SEs (`vcovCL`) |
| `spline_vs_linear` | `y_ge1` | primary formula with `w_act_steps` replaced by `ns(w_act_steps, 3)`; likelihood-ratio test against the linear model (2 df) |
| `random_slope` | `y_ge1` | primary formula `+ (0 + w_act_steps | person_id)` (uncorrelated random slope); LRT against the primary model, p halved (boundary, 1 df) |

For the sleep and heart-rate streams the tested term is still `w_act_steps` (yesterday's activity), with the stream's own `w_` terms in the model as above.

## 3. Gap-length models B and C (exploratory; one row per gap; same equations in every stream)

The **stream** only decides which gaps are used (gaps of the activity / sleep / heart-rate series); the **activity** values come from the activity table in all three:
`pre = w_act_steps` on the last observed day before the gap, `back1 = w_act_steps` on the first observed day after the gap, `back2 = w_act_steps` on the day after that
(`next observed date + 1`, which must itself be an observed activity day). `log2_gap = log2(gap length)`.

**Model B** (`lme4::lmer`, REML, Gaussian; effect = change per doubling of the gap, in within-person SD units; normal-approximation p)

    B_back2_unadj    : back2 ~ log2_gap                                   + (1 | person_id)
    B_back2_adj      : back2 ~ log2_gap + cohort + age10 + sex2           + (1 | person_id)
    B_back2_adj_pre  : back2 ~ log2_gap + pre + cohort + age10 + sex2     + (1 | person_id)
    B_back1_adj      : back1 ~ log2_gap + cohort + age10 + sex2           + (1 | person_id)

**Model C** (`ordinal::clmm`, cumulative logit; outcome = ordered gap-length group `1 | 2 | 3 | 4-7 | 8-14 | 15+`, groups with < 20 gaps merged into the next shorter one)

    logit P(length group <= j) = theta_j - ( beta * pre [+ gamma1 * cohort + gamma2 * age10 + gamma3 * sex2] ) - u_i ,   u_i ~ N(0, sigma_u^2)
    C_unadj : only beta * pre        C_adj : with the covariates

Run only if proportional odds is not rejected (LRT of the same model without random effect against one with a non-proportional effect of `pre`, p > 0.05). Effect: `OR = exp(-beta)` for a
longer gap per 1 SD below usual activity. If the random-intercept variance is on the boundary: the same cumulative-logit model without `u_i` and person-clustered bootstrap SEs.
(`cohort`, `age10`, `sex2` are dropped where constant.)
