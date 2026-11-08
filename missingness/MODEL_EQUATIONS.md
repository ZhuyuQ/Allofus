# Model equations, outcomes and predictors, by stream

Written from the code (`make_onset_grid`, `build_onset_formula`, `run_sensitivity_stage`, `fit_gap_dose_response`,
`fit_gap_length_ordinal`); the formulas below are the strings the code builds (checked by printing them).
A person is `i`, an observed day of that person in the stream is `d`.

## 0. Common definitions

**Observed day** (per stream): a row exists in that stream's day table (activity table; sleep table: a main sleep or a nap
counts, except in the `main_sleep_only` sensitivity where only `has_main_sleep` counts; heart-rate table: the four zone
rows collapsed to one row per day). **Gap**: >= 1 consecutive calendar days with no row between two observed days.

**Wear** (`wear.R`; only the wear analyses use it): `wear_min(d)` = heart-rate minutes recorded on the calendar day `d` (`total_hr_min`, 0-1440), the same measure for the activity, sleep and
heart-rate streams; NA = no heart-rate record on `d` (unknown wear). Threshold analyses: `wear_min(d) >= 60 h` with `h` = 16, 20 and 22 (`cfg$wear$hours`).

**People and covariates come from the full-pool selection** (`../selection/`): age = whole years at the person's earliest record over the three streams; CCI = `cci_calc_func()` with that age.

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

    m(d)   = mean of the person's observed x on days BEFORE d          (>= 7 such days, else the day has no z-score; "observed" = observed in the MODELLED stream: for the
                                                                          activity predictors of the sleep / heart-rate onset models and for the day BEFORE a gap in models B / C there, days that also have
                                                                          a record in that stream; the days AFTER a gap in model B and their usual use all activity days)
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

**Rows**: every observed day `d` of a person in the stream except the person's last, **whose calendar day `d - 1` also has a record**, with all predictors of the stream non-missing **on day `d - 1`** ("common risk set":
the stream's own predictors AND the three activity predictors on that date; for adjusted models also `age10`, `sex2`, `cohort`). **Every `w_` and `b_` in the formulas below is the value of the last FULL day `d - 1`**
(`cfg$onset$predictor_day = "last_full_day"`): day `d`, the last day with a record before the gap, may be a partial day. The adjusters (weekday of d + 1, study time, gap history) are those of row `d`.
**Outcome** `y_ge{k}(i,d) = 1` if a gap of at least `k` days starts on day d+1 (next observed date - d - 1 >= k), else 0.
Primary `k = 1`; model D `k = 2, 3, 7`.

**Estimation**: `lme4::glmer`, binomial, logit, `(1 | person_id)`, `nAGQ = 0`, bobyqa. Wald z tests; `OR` as above.

### Primary model (joint predictors, adjusted, k = 1) - the one pre-specified test: `w_act_steps`, Holm across streams

**Activity stream** (outcome: an *activity* gap starts tomorrow)

    y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa
          + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + device + cci + cci_missing
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Sleep stream** (outcome: a *sleep* gap starts tomorrow; predictors: the sleep and activity of the last full day)

    y_ge1 ~ w_sleep_inbed + w_sleep_eff + w_act_steps + w_act_sedentary + w_act_mvpa
          + b_sleep_inbed + b_sleep_eff + b_act_steps + b_act_sedentary + b_act_mvpa
          + cohort + age10 + sex2 + dow_t + ns(study_day, 3) + device + cci + cci_missing
          + prev_gap_len_l + any_prior_gap + days_since_return_l + (1 | person_id)

**Heart-rate stream** (outcome: a *heart-rate* gap starts tomorrow; predictors: the heart-rate coverage / elevated minutes and the activity of the last full day)

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
| `day_before_gap` | `y_ge1` | the primary formula with every `w_` and `b_` read from day d (the day right before the gap, which may be a partial day) instead of the last full day d - 1; every row with complete predictors on d; term reported: `w_act_steps` (of d) |
| `last_full_day_same_rows` | `y_ge1` | the primary formula (predictors of d - 1) on the rows that have a record on d - 1 and complete predictors on both d - 1 and d (identical rows in the next two) |
| `day_before_gap_same_rows` | `y_ge1` | the same rows and formula, but every `w_` and `b_` read from day d; outcome unchanged: a gap starts on d + 1 |
| `both_days` | `y_ge1` | the same rows; the primary formula with day d, `+ lag1_w_act_steps` (steps on d - 1); reported: both estimates and the Wald test of `b_d - b_(d-1)`, variance `v_d + v_(d-1) - 2 cov` |
| `main_sleep_only` (sleep only) | `y_ge1` | primary formula, but a night is "observed" only if it has a main sleep (`has_main_sleep`); gaps and risk set are rebuilt on that definition |
| `wear_ge16h`, `wear_ge20h`, `wear_ge22h` | `y_ge1` | the `day_before_gap` model (every `w_` and `b_` read from day d, rows with complete predictors on d), restricted to the rows with `wear_min(d) >= 960` / `>= 1200` / `>= 1320` (rows with no heart-rate record on d excluded); the formula, the standardization (the usual = all earlier days, >= 7) and the outcome are the primary ones; term reported: `w_act_steps` (of d); the rows and events kept are in `wear_filter_flow` |
| `conditional_logit` | `y_ge1` | `survival::clogit`, Efron: `y_ge1 ~ w_act_steps + w_act_sedentary + w_act_mvpa + dow_t + ns(study_day, 3) + prev_gap_len_l + any_prior_gap + days_since_return_l + strata(person_id)` (`strata` written without the `survival::` prefix so that survival recognises it) (sleep/hr: with their own `w_` terms too); only persons with >= 1 event; no `b_`, cohort, age, sex or year terms (they cancel or are time-invariant); `dow_t` kept only if it varies |
| `cluster_robust_logit` | `y_ge1` | the primary formula without the random intercept (pooled `glm`, binomial), person-clustered sandwich SEs (`vcovCL`) |
| `spline_vs_linear` | `y_ge1` | primary formula with `w_act_steps` replaced by `ns(w_act_steps, 3)`; likelihood-ratio test against the linear model (2 df) |
| `random_slope` | `y_ge1` | primary formula `+ (0 + w_act_steps | person_id)` (uncorrelated random slope); LRT against the primary model, p halved (boundary, 1 df) |

For the sleep and heart-rate streams the tested term is still `w_act_steps` (the activity of the last full day), with the stream's own `w_` terms in the model as above.

## 3. Gap-length models B and C (exploratory; one row per gap; same equations in every stream)

The **stream** only decides which gaps are used (gaps of the activity / sleep / heart-rate series); the **activity** values come from the activity table in all three:
`pre = w_act_steps` on the last FULL day before the gap, d - 1 (against the usual before that day, as in the onset model; `cfg$gap_models$predictor_day`; the last recorded day d is the sensitivity analysis, see the table below);
`back1`, `back2` = steps (log1p) on the first observed day after the gap (`next observed date`) and on the day after that (`next observed date + 1`, which must itself be an
observed activity day), **against the person's usual THROUGH the last full day before the gap, d - 1**: `back = (x_post - m_through_d-1) / sd_w` (data observed before the gap only, the possibly partial day d not included;
0 = usual, + = more steps). Gaps with fewer than 7 earlier observed days have no usual. With `standardization = "full_period"` (units `B_back2_adj_fullperiod`, `C_adj_fullperiod`) all of
these are whole-period z-scores. **Predictor-day sensitivity of B and C** (the last recorded day d instead of the last full day d - 1; same formulas, the day-dependent columns `back*`, `w_pre_day_*` read from the `alt_` columns):
`B_back2_adj_dbg_all`, `B_back2_adj_pre_dbg_all`, `C_adj_dbg_all` (every gap with day d); `*_same_full` / `*_same_dbg` (exactly the gaps that have both days, d - 1 / d); `C_adj_both_days`
(`length_bin ~ w_pre_day_steps [d - 1] + w_dbg_steps [d] + adjusters + (1 | person_id)`; reported: both estimates and the Wald test of `b_(d-1) - b_d`, variance `v + v' - 2 cov`). `post_k0 ... post_k7` = the same z-score on the return day and the 7 days after it (observed trajectory, `model_B_observed_trajectory`).
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

`pre_mvpa = w_pre_day_mvpa` = log(1 + MVPA minutes) on the last full day before the gap, d - 1 (the same day as `pre`), minus the mean of the person's EARLIER observed days (>= 7), divided by the
pooled within-person SD of MVPA: **per +1 SD higher MVPA than the person's usual**. Gaps without it are dropped from the `*_mvpa` fits. `delta` and the MVPA coefficient of B are reported
with 95% intervals in `gap_mvpa_effects` (B: SD of activity; C: odds ratio `exp(delta)` of a longer gap); the specifications without MVPA are unchanged.

Run only if proportional odds is not rejected (LRT of the same model without random effect against one with a non-proportional effect of `pre`, p > 0.05). Effect: `OR = exp(beta)` for a
longer gap per +1 SD higher than usual activity. If the random-intercept variance is on the boundary: the same cumulative-logit model without `u_i` and person-clustered bootstrap SEs.
(`cohort`, `age10`, `sex2` are dropped where constant.)

## 4. The run-in analysis (stages `runin` = primary, `runin_sensitivity`): the same models on the days / gaps AFTER the first N = 28 days

**Run-in measure** (per person and stream; `runin.R`): `m = N - (observed days among study days 1..N)`, NA if the follow-up is shorter than N days.
`cov0 = log2(1 + m) - mean over the people analysed`  (per doubling of 1 + m);  `cov0_cat` in {0, 1, 2-3, >=4} (reference 0);  `cov0_lin = 100 m / N / 10` centred.
**Rows**: onset: rows whose predictor day (the last full day d - 1) has `study_day > N`; models B and C: gaps whose predictor day (the last full day before the gap, d - 1) has `study_day > N`. Identical rows in every variant.

    PRIMARY (stage runin):
    onset, reference           : y_ge1 ~ <primary joint adjusted formula>                                  + (1 | person_id)
    onset, interaction         : ... + cov0 + w_act_steps:cov0                                             + (1 | person_id)      H1: coefficient of cov0, H2: coefficient of w_act_steps:cov0 (Wald)
    SENSITIVITY (stage runin_sensitivity):
    onset, main                : ... + cov0                                                                + (1 | person_id)      H1: coefficient of cov0
    (H2 sensitivity: LRT of the interaction model vs the main-effect model, 1 df)
    onset, categorical         : ... + cov0_cat                                                            + (1 | person_id)      shape check (AIC; LRT vs reference, 3 df)
    onset, linear_interaction  : ... + cov0_lin + w_act_steps:cov0_lin                                     + (1 | person_id)      functional-form sensitivity
    model B  (reference | main | interaction): back2_dev ~ log2_gap [+ cov0 [+ log2_gap:cov0]] + cohort + age10 + sex2 + device + cci + cci_missing + (1 | person_id)
    model C  (reference | main | interaction): ordinal cumulative logit, length group ~ w_pre_day_steps [+ cov0 [+ w_pre_day_steps:cov0]] + covariates + (1 | person_id)

Effect of the key term at `m` missing days (interaction variants): `b_key + b_int (log2(1 + m) - centre)`, variance `v_key + L^2 v_int + 2 L cov(b_key, b_int)`, `L = log2(1 + m) - centre`.
Onset and C: `OR = exp(beta)`; B: change in activity (SD of the person's usual steps) per doubling of the gap. Hypotheses H1 (`cov0`), H2 (key term x `cov0`, Holm over the streams), H3 (primary
estimate with and without the run-in terms): STATISTICS.md section 10.

## 5. Diagnostics

    V = covariance of the fixed effects without intercept / thresholds,  R = cov2cor(V)
    GVIF(term) = det(R_term) det(R_other) / det(R),  df = coefficients of the term,  VIF-equivalent = GVIF^(1/df)      (>= 5 high, >= 10 severe)
    joint Wald test of a term:  chisq = b_term' V_term^-1 b_term  on df = coefficients of the term


## 6. Weekly rhythm of gap starts and gap ends (stage `weekday_rhythm`; `R/weekday_rhythm.R`)

A separate analysis with its own rows and models; nothing in sections 1-5 changes. The formulas below are the strings the code builds (`wr_formula()`, printed
for every specification; `device`, `cci` and `cci_missing` appear only when they vary in the data). Terms are those of STATISTICS.md 3b.

**Weekday and its coding.** `k = 0` (Monday) ... `6` (Sunday) of the day that carries the outcome: the gap start day `d + 1` (equal to the onset model's `dow_t` for
activity and heart rate) or the return day `t`. The **sleep stream is aligned to the night it started**: `k` and the holiday lookup use `stored date - 1`, done only in
`wr_align_date()`. `weekday` is a factor with **sum-to-zero contrasts** (`contr.sum`): its 6 coefficients `weekday1 ... weekday6` are the log-odds of Monday ... Saturday
against the average day, and Sunday is minus their sum. `c1 = cos(2 pi k / 7)`, `s1 = sin(2 pi k / 7)` (the first weekly harmonic) enter only the random effects.
`ns1, ns2, ns3` = the natural spline basis (3 df) of `study_day`, computed on the rows of each fit.

**Gap-start rows.** One row per observed day `d` except each person's last observed day. Outcome `start = 1` if a gap starts on `d + 1`, else 0. Neither a record on `d - 1` nor
a 7-day baseline is required. `dsr_f` = days since return (`d - r + 1`, the return day `r` = 1; before the first gap `d - first record + 1`), a factor with levels `1, ..., 21, 22+`
(reference level 1). `prev_gap_len_l = log2(1 + length of the last completed gap)`, `any_prior_gap` as in section 0.

**Gap-end rows.** For a gap of length `L` that starts on `s`, one row for each `t = s + 1, ..., s + L` with days so far in the gap `dsf = t - s`; outcome `end = 1` only on the
return day `t = s + L`. Rows with `dsf` above 13 (14 days so far always end, because the longest gap is 14 days) are dropped and counted. `dsf_f` = a factor with levels `1, ..., 13`.
Adjusters of gap-end rows: `cohort, age10, sex2, ns1-ns3, device, cci, cci_missing` (no gap-history terms; days since return is not included).

**Random effects.** `(1 + c1 + s1 | person_id)` (`correlated`): a person intercept and a personal first weekly harmonic. If the correlated fit of `S1` / `E1` does not converge or is
singular, **every fit of that stream and row type** uses `(1 | person_id) + (0 + c1 | person_id) + (0 + s1 | person_id)` (`independent_harmonics`); `S5` always has `(1 | person_id)`.
All fits: `lme4::glmer`, binomial, logit, `nAGQ = 0`, bobyqa (`cfg$model`).

**Formulas** (`<RE>` = the random effects above; `A_s` = `cohort + age10 + sex2 + ns1 + ns2 + ns3 + device + cci + cci_missing + prev_gap_len_l + any_prior_gap`;
`A_e` = `cohort + age10 + sex2 + ns1 + ns2 + ns3 + device + cci + cci_missing`; the terms of the fit that are constant in the data are dropped).

| fit | rows | formula | used for |
|---|---|---|---|
| `S0` | gap start | `start ~ A_s + <RE>` | null of H1 (starts) |
| `S1` | gap start | `start ~ A_s + weekday + <RE>` | H1 (starts), H3 and H4 null, Figures 1, 3 |
| `S2` | gap start | `start ~ A_s + dsr_f + <RE>` | null of H2; days-since-return curve |
| `S3` | gap start | `start ~ A_s + dsr_f + weekday + <RE>` | H2, H2r, Figure 2 |
| `S4` | gap start | `start ~ A_s + weekday + weekday:cohort + <RE>` | H4 (starts), Figure 4 |
| `S5` | gap start | `start ~ A_s + weekday + (1 | person_id)` | null of H3 |
| `E0` | gap end | `end ~ A_e + dsf_f + <RE>` | null of H1 (ends) |
| `E1` | gap end | `end ~ A_e + dsf_f + weekday + <RE>` | H1 (ends), H4 null, Figure 1 |
| `E4` | gap end | `end ~ A_e + dsf_f + weekday + weekday:cohort + <RE>` | H4 (ends), Figure 4 |
| `S1_holidays`, `S3_holidays`, `E1_holidays` | as `S1`, `S3`, `E1` | the same formulas on the rows whose outcome day is not within +/- 1 day of a US federal holiday or its observed date | sensitivity |
| `S1_before_first_gap`, `S3_before_first_gap` | gap start | the formulas of `S1`, `S3` on the rows up to and including each person's first gap start (`any_prior_gap` and `prev_gap_len_l` are then constant and dropped; `dsr_f` = days since the first record) | sensitivity |
| `S1_activity_same_rows` | gap start | the formula of `S1` on the rows of the primary onset model (a record on `d - 1` and complete predictors) | sensitivity |
| `S1_activity` | gap start | `start ~ A_s + weekday + w_p + b_p (p = the stream's primary predictors of d - 1) + <RE>` on the same rows | sensitivity |

**Likelihood-ratio tests** (`chisq = max(0, 2 (logLik_big - logLik_small))`, df = the difference in the numbers of parameters of the two fits, as counted by `logLik`):

| id | big vs small | df | Holm family |
|---|---|---|---|
| `H1s` | `S1` vs `S0` | 6 | H1 (3 streams x {start, end}) |
| `H1e` | `E1` vs `E0` | 6 | H1 |
| `H2` | `S3` vs `S2` | 6 | H2 (3 streams) |
| `H2r` | `S3` vs `S1` | 21 | none (reported) |
| `H3` | `S1` vs `S5` | 5 (correlated) or 2 (independent harmonics) | H3 (3 streams) |
| `H4s`, `H4e` | `S4` vs `S1`, `E4` vs `E1` | 6 x (cohorts - 1) | H4 (3 streams x {start, end}) |

**Population-average profile.** For a weekday `j`, with `x_i` the observed covariates of a seeded sample of `n_profile_rows` rows (weekday set to `j`; for `S4` / `E4` the rows of one cohort),
fixed effects `b` (the estimate, or one of `n_beta_draws` draws from `N(b, Vcov)`), random effects `u_m = (u_0, u_c, u_s)` (`n_re_draws` draws from `N(0, Sigma)`, one set shared by all weekdays and draws of `b`):

    p_j(b) = mean over rows i and draws m of  plogis( x_i(j)' b + u_0m + u_cm cos(2 pi j / 7) + u_sm sin(2 pi j / 7) )
    ratio_to_average_day_j = p_j / mean_j(p_j)          peak_to_trough = max_j p_j / min_j p_j
    percent_j = 100 p_j;  95% intervals = 2.5% and 97.5% quantiles over the draws of b (Sigma held fixed); the point estimate uses the estimated b.

These are ratios of predicted probabilities, not odds ratios or risk ratios of a person. The **days-since-return curve** is the same quantity with `dsr_f` set to each level and the weekday as observed
(the random-effect term uses each row's own `c1`, `s1`), `n_curve_draws` coefficient draws.

**H2 shrinkage.** `1 - log(peak_to_trough of S3) / log(peak_to_trough of S1)` (point estimate): the share of the weekday pattern explained by days since return.

**Personal rhythm (from `S1`).** With `e_k` the 7 weekday effects (`e_k = weekday_{k+1}` for `k = 0..5`, `e_6 = -sum`):

    a = (2/7) sum_k e_k cos(2 pi k / 7)          b = (2/7) sum_k e_k sin(2 pi k / 7)          shared peak day = atan2(b, a) 7 / (2 pi)  mod 7
    person rhythm = (a + u_c, b + u_s),  amplitude A = sqrt((a + u_c)^2 + (b + u_s)^2),  peak day = atan2(b + u_s, a + u_c) 7 / (2 pi) mod 7,  peak-to-trough = exp(2 A)  (odds scale)

Reported from `n_person_draws` (10,000) seeded draws of `(u_c, u_s)` from their estimated covariance (model based, not shrunk): the SDs of `u_0`, `u_c`, `u_s`, the median and the 10th and 90th percentiles of
`exp(2 A)`, and the share of people whose peak day lies within one day of the shared peak day (circular distance). The BLUPs of `(u_c, u_s)` are only binned for Figure 3 (shrunk toward the shared rhythm).
