# Stream states: which streams are present on each day

Descriptive analysis of the three Fitbit streams (activity, heart rate, sleep) day by day: how often each combination occurs, what the combinations mean, and how the removal and
filtering rules affect the cohorts. It reports **estimates and their uncertainty only, no significance tests**, and it never relates the steps of the previous day to the states.
Code: `R/stream_states.R`, `R/plots_stream_states.R`; run: `pipeline/run_stream.R` (notebook cells); tests: `tests/testthat/test-stream.R`. No existing file was changed.
Output folder: `<results>/stream/` (`tables/`, `figures/`, `methods.md`); the check without unknown-device days: `<results>/stream/no_unknown_device/`.

## Methods

**People and days.** The people of the existing selection (no re-selection). One row per person-day. A person's **window** runs from the first to the last record in *any* of the three streams.

**Removed and counted, not treated as missing.** A day is removed for the **first** of these reasons, in this order (so the counts add up):

1. *Stream not active*: the day lies before that stream's first record or after its last (any of the three streams).
2. *Device without heart rate*: the device **used on that day** has no heart-rate sensor. The device of a day is the tracker record of the person on or before the day (the latest record wins; two versions on the
   same record date: the alphabetically last). Records come from the **raw device export**, because clip-on trackers and phone tracking are not in the L1 watch timeline.
3. *Device unknown*: no tracker record on or before the day. The day is **kept if the person's usual device** (the version with the most records, ties: the most recent record) **measures heart rate**; otherwise it is
   removed. A person without any device record has no usual device, so these days are removed (setting `keep_unknown_if_no_usual`). **Check:** the whole analysis is repeated without *every* day of unknown device
   (`stream/no_unknown_device/`).

Which versions have no heart-rate sensor is **product knowledge, not data** (`st_config()$no_hr_versions`: One, Zip, Ultra, Classic, Flex, Flex 2, Force, Alta, Charge, Inspire, Ace, Ace 2, Ace 3, MobileTrack, exact names after
removing "Fitbit "); a version that is not listed is assumed to measure heart rate. The table `stream_device_hr_lookup` lists every version seen with the classification used: **check it**.

**Heart rate switched off.** There is no separate rule for heart rate being switched off. The selection already required the longest heart-rate gap to be at most 14 days, so long switched-off periods were removed before this analysis.

**Activity filter.** An activity row is set to missing when sedentary minutes cover the whole day and there is no heart-rate record that day. "The whole day" is 1,440 minutes, 1,380 on the second Sunday of March (clocks go forward) and
1,500 on the first Sunday of November. The time zone of a person is not known, so the **shortest possible length of the date** is used (1,380 on the March date, 1,440 otherwise) and the other active minutes (lightly, fairly, very active) must be
zero. If sedentary minutes cover the whole day but heart rate was recorded, the row is **kept and flagged** (worn but motionless).

**Presence after the filter.** Activity: a row remains. Heart rate: a record with **at least one minute**. Sleep: a **main** sleep record dated that day (naps do not count; Fitbit dates sleep by the wake-up day, so it describes the night before).

**States.** 1 all present | 2 all missing | 3 non-wear (activity present, heart rate and sleep missing) | 4 night missing (activity and heart rate present, sleep missing) | 5 other (every remaining combination). **Worn days** are states 1 and 4.
The expected nesting is sleep within heart rate within activity: the eight combinations are 4 nested ones (111, 110, 100, 000), "activity + sleep, heart rate missing" (101) and 3 other unexpected ones (011, 010, 001).
**Habitual night non-wear** = the share of a person's worn days that are state 4, for people with at least 30 worn days; **description only**: it is computed from the whole follow-up and would leak future information into a model.

**Uncertainty.** Shares of days: 95% intervals that account for repeated days of the same person (ratio estimator with people as clusters, `cluster_mean_ci`). Averages of person-level values: t intervals over people. Cells that represent
1-19 people are not shown (`<20`, values NA). Dates are used as stored (the local date against UTC is checked in Part C). All cohort differences are **unadjusted**: cohorts differ in age, device and length of follow-up, so differences are patterns, not effects.

## Parts, tables and figures (`stream_<name>.tsv`, `stream_<part>_<what>.png`)

| part | content | tables | figures |
|---|---|---|---|
| A data flow | by cohort: people and person-days; days removed by each reason (shares with clustered intervals); `n_days_at_step` = days with the characteristic, `n_window_days` = all days of the cohort | `flow`, `people` | `a_flow` (stacked shares), `a_reasons` (share per reason, with intervals) |
| B activity filter | activity rows set to missing (count, share of activity rows), people affected, worn-but-motionless days, spread of the person-level share | `activity_filter`, `activity_filter_people`, `activity_filter_spread` | `b_filter_share`, `b_filter_people` (ECDF over people) |
| C patterns and nesting | share of days in each of the eight combinations before and after the filter ("activity + sleep, heart rate missing" separate from the other unexpected ones); **before the filter**: heart rate present when activity is missing, sleep present when heart rate is missing; date check (streams shifted by -1, 0, +1 day); gate | `combos`, `nesting_checks`, `date_check`, `gate` | `c_combos`, `c_checks` |
| D what state 3 is | steps on state-3 days; state-3 against the person's own state-1 days **within the same type of day** (weekday / weekend), people with at least 3 days of each, steps and MVPA minutes, mean difference with a t interval over people | `steps_state3`, `state3_vs_state1` | `d_steps_state3`, `d_state3_vs_state1` |
| E composition | share of days in each state **pooled over days** (long follow-up weighs more) and as the **person average** (everyone weighs the same), overall, by age group, by device family; a random sample of people day by day; habitual night non-wear | `composition`, `night_habit` | `e_composition`, `e_composition_age`, `e_composition_device`, `e_day_by_day`, `e_night_habit` |
| F what "night missing" contains | nights where the day before and the day itself were worn days: share of state-4 nights (state-1 nights for comparison) with heart-rate wear of at least 20 and 22 hours on both days ("probably worn": likely, not proven) | `probably_worn` | `f_nights_hr` (two-day heat map), `f_probably_worn` |
| G partial wear | heart-rate hours on worn days (distribution, share below 10 and 20 hours); share of heart-rate records whose minutes exceed the length of the day | `hr_wear`, `hr_wear_below`, `minutes_check` | `g_hr_wear`, `g_hr_below`, `g_minutes_check` |

Person-level and day-level data that only feed a figure (person shares, steps of state-3 days, the sample of people, nights) are **never written to a table**. The figure `e_day_by_day` is an individual-level picture (one row per person, no identifiers): check the
disclosure rules before it leaves the workspace.

## The gate (Parts A-C first)

Cell 4 of `pipeline/run_stream.R` runs Parts A-C and the minutes check of Part G, then `st_require_gate()` **stops** unless both pass; Parts D-G are not run.
* **Date check** (Part C): the unexpected combinations other than "activity + sleep, heart rate missing" must be at most **1%** of analysed days, and the best date shift of sleep and of heart rate must be 0. A best shift of +1 means the table is one day late, -1 one day early (local date against UTC). `stream_date_check` shows the share for each shift.
* **Minutes check** (Part G): at most **0.1%** of heart-rate records may have minutes above the length of the day (1,440; 1,500 where clocks go back). The 0.1% is my choice for "almost zero" (`minutes_exceed_max`). If it fails, the heart-rate zones were combined incorrectly in the data preparation and heart-rate minutes cannot be used as a wear measure until this is fixed.

## If the minutes check fails

The gate stops Cell 4 and says which check failed. For the minutes check: Cell 4b (read-only) shows how far above the length of the day the records are (a ratio near 2 means rows counted twice), how concentrated they are in people and years, and, from the raw heart-rate export,
the share above the length of the day under three ways of combining the rows of a zone: sum of all rows (what the ETL does), drop exact duplicate rows then sum, largest row of each zone then sum. If one rule brings the share to about zero, Cell 4c uses those minutes **in this analysis only**
(`st_use_hr_minutes`: it redoes the activity filter, the presence flags and the removal reasons) and Cell 4 is run again. The L1 tables and the other analyses are not changed: their heart-rate minutes (for example the wear thresholds of the weekday analysis) would still contain the problem,
which is a reason to fix the ETL (`standardize_hr`) and rebuild L1 if the cause is the summing of repeated rows. Shares are clipped to 0-100%.

**Option 2 (Cell 4c): the raw minutes themselves are a little above the day.** If Cell 4b shows no duplicated zone rows (all three rules give the same share) and ratios of 1.0-1.1, the cause is not the ETL: a calendar day can be longer than 24 hours (travel across time zones,
daylight-saving dates outside the US rule used here). Then `minutes_tolerance = 0.1` accepts records up to 10% above the length of the day (about 26-27 hours) in the check, and `cap_hr_minutes = TRUE` caps the heart-rate wear of Parts F and G at the length of the day.
Records further above still count against the 0.1% limit. This is a judgement about the gate, made after looking at the data; state it in the Methods when it is used.

## Decisions to check

1. The list of device versions without heart rate (`stream_device_hr_lookup`); unlisted versions are assumed to measure heart rate.
2. A day of unknown device of a person without any device record is removed (literal reading of "keep if the usual device measures heart rate"); the check without unknown-device days removes all of them.
3. Concurrent devices: the alphabetically last version on the latest record date is "the device of the day".
4. Heart rate present = at least one minute; a heart-rate row with 0 or missing minutes counts as missing.
5. The whole-day rule uses the shortest possible day length of the date, because the time zone is unknown.
6. A day on which **any** stream is not active is removed for all streams (the state needs all three).
7. Part D: eligibility needs at least 3 days of each state **within the type of day**; the overall estimate is the average of the type-specific differences weighted by the person's state-3 days of each type.
8. Part E person average: every person with at least one analysed day; strata by `age_group` (age at the first record) and the main device family.
9. Part F: "probably worn" = the **smaller** of the two days' heart-rate hours is at least the cut-off (20 or 22 hours).

## Limits

Only descriptive; dates as stored; the time zone and the device in use on a day are inferred (as-of join of record dates); state 3 mixes a watch that was moved but not worn, steps counted by a phone, and a watch worn with heart rate switched off (Part D separates them by the person's own steps
and MVPA but cannot prove which); "probably worn" nights are likely, not proven; shares are unadjusted. Not run on real data yet. Tests use hand-counted data and simulated data (`tests/testthat/test-stream.R`).

# Stream model: does a high-step day d predict the status of day d + 1?

Code: `R/stream_model.R`, `R/plots_stream_model.R`; run: Cells 7-11 of `pipeline/run_stream.R`; tests: `tests/testthat/test-stream-model.R`. Output: `<results>/stream/model/` (`tables/stream_model_*.tsv`, `figures/`).

**Hypotheses.** H1 (primary, two-sided Wald test): more steps than usual on day d change the odds that d + 1 is **all missing** (expected OR < 1). H2: this association differs between all missing and activity only, and between all missing and sleep not logged.
H3 (negative control): steps are not associated with category 4, which cannot arise from non-wear.

**Time order.** Exposure = steps of calendar day d; outcome = the status of day d + 1 (the sleep ending on d + 1 usually starts on the evening of d: the 20-hour wear rule limits any overlap to at most 4 hours). Day d was chosen as the exposure day because the earlier per-stream analyses found the association concentrated on the last day.

**Risk set.** One row per person-day d with activity, heart rate and a main sleep record on d (fully recorded), at least 20 hours of heart-rate wear on d (sum of the zone minutes), and at least 7 earlier days that met both conditions. Days removed in the stream-state analysis are not days; rows whose d + 1 is an excluded day or that are the last day of the window are dropped and counted (`flow`).

**Outcome (d + 1).** 0 fully recorded (reference) | 1 all missing (main outcome) | 2 activity only | 3 sleep not logged (mixes nights without wear, nights worn but not detected, and dates with no sleep ending) | 4 any other pattern (activity + sleep without heart rate, sleep only, heart rate without activity; negative control; its composition is in `category4_composition`).

**Exposure and covariates (the onset model's settings).** x = log(1 + steps) on d; the person's usual m = the mean of x over their **earlier** qualifying days only (no future days); w = (x - m)/SD within, b = (m - mean)/SD between; the constants are taken over all risk-set rows.
Covariates: cohort, age, sex, CCI score (and a missing indicator), device family, weekday of d + 1, recording history (length of the last spell, any earlier spell, days since the last return), a natural spline of the day of follow-up. No calendar year.

**Model A.** One mixed logistic model per category (k against 0, random person intercept, the onset model's engine settings): subject-specific odds ratios per +1 SD of w. H1 = Wald test for category 1 (unadjusted). Equivalence for categories 2-4: two one-sided tests, 90% interval within 0.97-1.03; the Holm adjustment over the three is added (`p_tost_holm`), the unadjusted p is kept.
Negative-control verdict (`negative_control`): pass = 90% interval within the margin with at least 6,000 events; fail = the whole 90% interval outside the margin (main result possibly confounded); otherwise inconclusive.
**Model B.** One multinomial logit (reference 0, all rows) fitted by Newton-Raphson with person-clustered (sandwich) errors. H2 = the difference of the w coefficients of categories 1 and 2 and of 1 and 3 (ratio of odds ratios; Holm over the two, unadjusted p kept). Predicted probabilities (all five statuses) across w from -2 to +2 average the observed covariates; intervals from 2,000 draws of the coefficients with their cluster covariance.

**Absolute risk, two ways.** Model B (population-averaged) and model A (the random intercept integrated out on a quantile grid). Both give the probability across w and the risk difference between +1 and -1 SD with an interval (`risk_diff`).

**Event thresholds.** An effect claim needs at least 1,500 events of the category, an equivalence claim at least 6,000 (`event_counts`; look at it before the models). Cohort-specific fits are run only for cohorts with at least 1,500 all-missing events.

**Sensitivity.** More than 16 h and more than 22 h of heart-rate wear on d (the baseline days use the same rule), first 14 days excluded, first transition per person only, conditional logistic regression for category 1 (`checks_clogit`), lme4 against glmmTMB on a seeded subsample (`checks_engine`), and category 3 split into "probably worn" and "probably not worn" at 20/22 h on both days.

**Figures.** `stream_model_events`, `_observed_deciles`, `_forest` (odds ratios with the 0.97-1.03 band and the multinomial contrasts), `_absolute_risk` (models A and B), `_mix_by_steps` (stacked), `_sensitivity`.

**Limits.** Associations, not causes; subject-specific (A) and population-averaged (B) odds ratios differ by construction; category 3 is a mixture; a low-step evening can be the start of removing the watch (the wear rule and the 22 h variant address this partly).
