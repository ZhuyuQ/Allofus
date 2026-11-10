# Missingness pipeline (Fitbit activity, sleep, heart rate)

Day-level missingness for the four cohorts (Insomnia, Diabetes, Prostate, Healthy): summary
statistics and cohort comparisons, temporal structure, day-of-week clustering, a
mixed-effects model of **gap onset ~ activity of the last full day before the gap** (predictors standardized only with data observed before each
day), gap-level models of activity after a gap and of gap length, the same models on the days after a **run-in** (wear in study days 15-42; days 1-14 are skipped for the novelty effect), sensitivity analyses
and collinearity diagnostics. (Disease-label prediction and the mixture model are not part of this.)

**Input.** The people come from the full-pool selection in `../selection/` (see its `README.md`): CDR C2024Q3R9, >= 21 record days, longest gap <= 14 days and
compliance >= 70% (first to last record) in **each** of the three streams, CCI <= 10, with the age and the CCI taken at the person's **first record over the three
streams**. The list is `result/selection/20261031/final_pid_base.txt` (1,388 Insomnia, 1,306 Diabetes, 253 Prostate, 6,381 Healthy = 9,328 memberships) and the CCI
file `cci.txt` next to it; the 16 BigQuery exports of these people are in the bucket folder `20261031`. Nothing of the earlier exports or L1 tables is reused. On the
Workbench this folder is named `final_missing` (zip `final_missing.zip`); in the repository it is `missingness/`.

**Where to read what:** this file = layout, how to run, the analysis workflow and a map of every output (tables and figures) to its analysis;
`STATISTICS.md` = the statistical methods (definitions, standardization and its reference distributions, models, effect sizes, tests,
collinearity, run-in analysis, validation); `MODEL_EQUATIONS.md` = the exact formulas, predictors, outcomes and contrasts per stream.

## Layout

```
load_missingness.R        source everything in dependency order
R/config.R                make_config(): every setting that affects a result + cfg$resources (memory, workers)
R/utils.R                 checkpoints (run_unit), isolated workers, memory diagnostics, lazy objects, small helpers
R/io.R                    raw exports (L0), day tables (L1), the day store (reads one stream / some columns at a time)
R/bq_export.R             step 0: the SQL for the 16 BigQuery exports (activity, sleep, hr, meta); optional (Cell 6c): Fitbit device version export + probe + coverage check (not read by the ETL yet)
R/etl.R                   raw stream -> one row per person-day (activity, sleep, hr)
R/persons.R               cohort overlap, sex, age at first record, race, CCI groups
R/gaps.R                  gap runs, person summaries, the onset "risk set"
R/predictors.R            within-/between-person prior-day predictors (pre-gap or full-period standardization), per-stream bundle
R/collinearity.R          GVIF and joint Wald tests of every fitted model, the diagnostics tables
R/summarise.R             cohort tables and comparisons
R/temporal.R              onset vs continuation, study week, year, weekday
R/descriptive.R           summary statistics of metric values by cohort / strata
R/onset_model.R           primary model, model A (gap history), model D (k = 2, 3, 7); one checkpointed unit per model
R/sensitivity.R           sensitivity analyses of the primary model; one checkpointed unit per analysis
R/gap_models.R            models B and C (use gap length), their full-period variants and the observed post-gap trajectory; independent units per stream
R/outcome_variables.R     person-level % of days missing against each model variable: tests, adjusted odds ratios, figures
R/device.R                Fitbit watch version: person table + timeline from the device export, as-of assignment to days, missing / composition statistics
R/weekday_rhythm.R        weekly rhythm of gap starts and gap ends: rows, fits (S0-S5, E0, E1, E4 + sensitivity), tests H1-H4, profiles, personal rhythm, descriptions, results summary
R/plots_weekday_rhythm.R  its figures
pipeline/validate_weekday_rhythm.R   recovery of a planted weekday / cycle / personal rhythm on simulated data
R/wear.R                 wear time of the last day before a gap (heart-rate minutes recorded): description by stream and cohort, rows kept by the wear sensitivity analyses
R/runin.R                run-in analysis: wear in study days 15-42 (description), onset / B / C models after the run-in with and without it, tests, shape check
R/plots.R                 ggplot2 figures in the style of the earlier notebooks + figure captions
R/plots_effects.R         comparison forests (primary vs sensitivity / run-in variant), collinearity, observed post-gap, run-in figures
R/report.R                tables, captions, results summary, run metadata
R/pipeline.R              run_etl() and run_missingness_pipeline() (orchestration only)
R/simulate.R              seeded fake data with a planted effect (tests / demo)
pipeline/notebook_cells.R the notebook, Cells 1-19 (copy each "# ==== CELL n" block into a cell)
pipeline/validate_simulation.R  validation of the primary estimate on simulated data with a known effect (bias, SE calibration, coverage)
pipeline/run_exports.R    step 0, only if the CSV exports are missing (e.g. after a migration)
pipeline/run_etl.R        run once per export     |  pipeline/run_missingness.R   the analysis
pipeline/run_band.R       the pipeline on ANOTHER compliance group (50-60%): selection, exports, ETL, primary onset model, sensitivity analyses, and the comparison with the >= 70% run
pipeline/compare_bands.R  the comparison itself (tables and figures of two finished runs; used by run_band.R Cell 8)
pipeline/run_demo.R       end-to-end on simulated data
tests/                    Rscript tests/testthat.R   (about 5 minutes)
```

## Run

```r
# 0. only if the CSV exports do not exist: BigQuery -> 16 CSV exports in the workspace bucket
source("missingness/pipeline/run_exports.R")
# 1. once per export (workbench): raw CSVs -> compact day tables + persons table
source("missingness/pipeline/run_etl.R")
# 2. the analysis (reads only the compact tables)
source("missingness/pipeline/run_missingness.R")
```
In the notebook use `pipeline/notebook_cells.R` (Cells 1-19). Inputs besides BigQuery: `final_pid_base.txt` (person -> cohort) and `cci.txt`, both written by the
selection (`result/selection/20261031/`). Use the same export date (`20261031`) in `run_exports.R`, `run_etl.R` and the notebook. The notebook writes `L1` and the results to
`result/final_missing/{L1,results}`; `run_etl.R` and `run_missingness.R` still default to `result/final/missingness/...`: set the same folders in all of them.
Finished exports are skipped, so a failed run can simply be repeated.
No real data are needed to check the code: `Rscript pipeline/run_demo.R <dir>` runs everything on
simulated data whose true effect is known (OR 0.61 per +1 SD higher steps than usual).

## Data flow

```
L0 raw BigQuery CSV (16 exports)  ->  run_etl()  ->  L1 day_activity / day_sleep / day_hr + persons
L1  ->  run_missingness_pipeline()  ->  L2 (below)

result/final/missingness/
|-- L1/                       day_activity / day_sleep / day_hr, persons.rds, manifest.tsv, dropped_columns.tsv
`-- results/
    |-- tables/               every result table (*.tsv), incl. unit_failures, checkpoint_status, unit_run_status
    |-- figures/              *.png, figure_captions.md (what each figure shows and how to read it) and figure_index.tsv
    |                         (each figure with its analysis step, model and caption)
    |-- checkpoints/          one .rds per finished model (see below)
    |-- results_summary.md    the main numbers in plain language, written from the saved tables
    `-- run_config.json, session_info.txt, input_manifest.tsv
```
A day is **observed** when a row exists in its day table. Sleep: main sleep and naps are collapsed
to one row per day (a nap-only day counts as observed). Heart rate: the four zone rows (Out of
Range, Fat Burn, Cardio, Peak) become one wide row per day. Columns that are mostly NA (> 90%) or
constant are dropped, and the list is written to `dropped_columns.tsv`.

## Analysis workflow and output map

The analyses follow one logical order; every step reads the output of the earlier ones and the notebook cell that runs it is given.
Standardization: every within- and between-person z-score uses only data observed **before** the day (or gap) it is applied to (`cfg$standardization`,
see `STATISTICS.md` 5.2); the full-observation-period version is a sensitivity analysis. Effect sizes: odds ratios (OR = exp(coefficient)) **per +1 SD higher
than the person's usual** (model B: change in activity in SD units); reference levels Healthy (cohort), Female, Monday (weekday), the largest non-Missing Fitbit device family (`device`; small wrist families pooled as "Other wrist-worn"); age per 10 years (not centred), CCI score per +1 point, previous gap length and days since return per doubling of (1 + days) (not capped). Calendar year is **not** a covariate.

**Data sources and sample flow.** BigQuery exports (activity, sleep, heart rate, meta; Cells 6-7) -> day tables `L1/` (one row per observed person-day; Cell 8) ->
people per cohort after the overlap rule (Cell 5; `cohort_overlap.tsv`) -> the upstream selection (>= 70% observed, longest gap <= 14 days, >= 21 days; checked in
`selection_check`) -> models (the rows, people, events and gaps of every model are in `onset_models`, `model_B_post_gap_activity`, `model_C_gap_length`; days without
7 earlier observed days have no z-score and are not in a model) -> people, rows and events after the run-in (`runin_flow`, `runin_quality`).

| step | analysis (cell) | tables (`results/tables/`) | figures (`results/figures/`) |
|---|---|---|---|
| 1 | **Descriptive statistics** of missingness and of the metric values (Cell 11, 11b) | `cohort_missingness`, `cohort_tests_overall`, `cohort_tests_pairwise`, `missing_fraction_models`, `missing_pct_by_stratum*`, `metric_summary*`, `selection_check` | `heatmap_<stream>`, `missing_pct_by_cohort`, `longest_gap_by_cohort`, `gap_length_distribution`, `gap_ecdf`, `missing_pct_by_sex` / `_age_group` / `_race_plot` / `_cci_group`, `metrics_<stream>_by_cohort`, `metrics_<stream>_by_<stratum>` |
| 2 | **Temporal structure** (Cell 11, 11b) | `transition`, `study_week_rates`, `calendar_year_rates`, `weekday_rates`, `weekday_tests`, `weekday_odds_ratios`, `missing_trend_tests`, `missing_trend_heterogeneity`, `missing_weekday_*`, `age_rates` | `transition`, `missing_by_study_week`, `missing_by_year`, `missing_by_weekday`, `missing_by_age` |
| 1 | **Descriptive: person-level % of days missing against every model variable** (Cell 11e) | `outcome_variable_tests`, `outcome_variable_by_cohort`, `outcome_variable_contrasts`, `outcome_variable_models`, `outcome_variable_summary`, `outcome_variable_cells` | `outcome_by_variable_<variable>` (one per variable: cohort, sex, watch, age_first, act_steps, ..., streams side by side) |
| 3 | **Data quality: watch version** (Cells 6c, 8b, 11d; person level, no model) | `device_missing_by_cohort`, `device_tests`, `device_missing_vs_reference`, `device_missing_model`, `device_composition_by_cohort`, `device_composition_residuals`, `device_version_by_cohort`, `device_version_tests`, `device_switching_by_cohort` | `device_version_missing_and_composition`, `device_version_within_cohort` |
| 2b | **Weekly rhythm of gap starts and gap ends** (Cell 11g; stage `weekday_rhythm`; own results folder `<results>/weekday/`): shared and personal weekday patterns, calendar routine against charging cycle (days since return), cohort differences; tests H1-H4 | `weekday_rhythm_tests`, `_profile`, `_summary`, `_shrinkage`, `_dsr_curve`, `_personal`, `_coefficients`, `_fits`, `_flow`, `_holidays`, descriptions `_crude`, `_gap_length`, `_dsr_cells`, `_dsr_by_device`, `_first_record_weekday`; `results_summary.md` | `weekday_rhythm_profile`, `weekday_rhythm_calendar_vs_cycle_<stream>`, `weekday_rhythm_personal_<stream>`, `weekday_rhythm_cohort`, `_cohort_level`, `weekday_rhythm_sensitivity_start` / `_end`, supplements `weekday_rhythm_gap_length`, `_dsr_counts_<stream>`, `_dsr_by_device` |
| 3 | **Data quality: device version** (Cells 6c, 8b, 11d) | see the device row below | `device_version_missing_and_composition`, `device_version_within_cohort` |
| 3 | **Data quality: wear time of the last day before a gap** (Cell 11f; heart-rate minutes recorded on day d, all streams, no model) | `wear_last_day_distribution`, `wear_last_day` (median / quartiles, share under 16 h and 20 h with a person-clustered 95% interval, days with no heart-rate record), `wear_last_day_contrast` (last day minus ordinary days, cluster-robust test, Holm over cohorts) | `wear_last_day` (ECDF, streams in panels, cohorts in colour, solid = last day, dashed = ordinary days) |
| 4a | **Primary model: gap onset** (primary, model A inside it, model D, exploratory grid; Cells 12-14, 17) | `onset_coefficients` (every term), `onset_models`, `primary_test`, `onset_term_tests`, `vif`, `cohort_interaction`, `cohort_specific`, `risk_unadjusted`, `onset_deciles`, `onset_pre_gap_trajectory` | `onset_primary_forest`, `onset_forest_year_<stream>` + `onset_forest_other_<stream>` (all terms), `onset_deciles_<stream>`, `pre_gap_activity_observed_<stream>` (observed steps on the days before a gap, by gap length; only days followed by a gap), `risk_curve`, `cohort_forest`, `forest_gap_definition`, `gap_definition_all_terms_other_<stream>` + `_year_`, `exploratory_grid_effects_<stream>` |
| 4b | **Gap-level model B: activity after a gap** (Cell 15) | `model_B_post_gap_activity`, `model_B_all_terms`, `model_B_curve`, `model_B_observed`, `model_B_observed_trajectory` | `gap_dose_response`, `gapB_forest_<stream>`, `post_gap_activity_observed_<stream>` |
| 4c | **Gap-level model C: gap length** (Cell 15) | `model_C_gap_length`, `model_C_all_terms`, `model_C_predicted_probs` | `gap_length_probs`, `gapC_forest_<stream>` |
| 4d | **MVPA before the gap as a covariate of models B and C** (Cell 15; units `B_back2_adj_mvpa`, `B_back2_adj_pre_mvpa`, `C_adj_mvpa`) | `gap_mvpa_effects` (MVPA effect with 95% CI, the main terms with and without MVPA, n, VIF); the MVPA fits are also in `model_B_post_gap_activity`, `model_C_gap_length`, `model_B_all_terms`, `model_C_all_terms` | `gapB_mvpa_effect_<stream>`, `gapC_mvpa_effect_<stream>` |
| 5 | **Run-in analysis: wear in study days 15-42 and gaps later** (Cells 11c description, 18 primary models = without a run-in term / with the run-in and its interaction, 18b sensitivity = other forms; only days / gaps AFTER the run-in; H1 wear trait, H2 steps effect depends on run-in, H3 robustness) | `runin_flow`, `runin_distribution`, `runin_rate`, `runin_steps_by_group`, `runin_association`, `runin_onset_coefficients`, `runin_onset_models`, `runin_gapB` / `runin_gapC` (+ `_all_terms`), `runin_tests`, `runin_shape`, `runin_effect_by_missing`, `runin_primary_effect`, `runin_quality` | `runin_distribution`, `runin_rate_by_category` (H1), `runin_steps_by_group_<stream>` (H2), `runin_primary_effect` (H3), `runin_steps_by_missing_days` (H2 model), `runin_onset_all_terms_other_<stream>`, `runin_gapB_all_terms_<stream>`, `runin_gapC_all_terms_<stream>`, one GVIF figure per model `runin_gvif_<onset|model_B|model_C>_<stream>`; sensitivity (Cell 18b): `runin_sensitivity_primary_effect`, `runin_sensitivity_shape`, `runin_sensitivity_{onset,gapB,gapC}_all_terms_*` |
| 6 | **Sensitivity analyses** (incl. the wear analyses `wear_ge16h` / `wear_ge20h` / `wear_ge22h`: predictors of day d on the rows whose day d has >= 16 / >= 20 / >= 22 h of heart-rate minutes, table `wear_filter_flow`; and the predictor-day check of models B and C: `gap_day_check`, `gapB_day_check_<stream>`, `gapC_day_check_<stream>`) incl. full-period standardization (Cell 16; for models B and C the `*_fullperiod` units of Cell 15) | `sensitivity`, `sensitivity_all_terms`, `sensitivity_collinearity`, `day_check`, `day_descriptive`, `day_pattern` | `forest_sensitivity`, `sensitivity_forest_<stream>`, `sensitivity_all_terms_other_<stream>`, and the predictor-day check (the primary model reads the last full day d - 1; sensitivity: the day right before the gap d): `day_check_<stream>`, `day_pattern_<stream>` |
| 7 | **Model diagnostics and integrated summary** (written by every model cell) | `model_collinearity` (GVIF of every model of the main analysis), `model_term_tests`, `onset_collinearity`, `gap_collinearity`, `unit_failures`, `checkpoint_status`, `unit_run_status`; `results_summary.md` | `collinearity_gvif_<model>_<stream>` (one figure per model) |

`figures/figure_index.tsv` lists every figure present with its step, analysis, model and caption; `figures/figure_captions.md` has the captions;
`tables/table_index.tsv` lists every table with its step.
The all-terms figures come in two parts: `..._year_...` = the calendar-year terms (each year, or pooled years, against the year with the most gap starts), `..._other_...` = every other term
(within-person, between-person, cohort, age, sex, weekday, gap history); the study-time spline and the intercept are not drawn (joint tests: `onset_term_tests`,
`model_term_tests`). In all comparison figures the first (dark) colour is the primary analysis.

**Settings that change results** (`cfg`, part of every checkpoint key): `standardization` (`method = "pre_gap"`, `baseline_min_n = 7`), `reference_cohort`,
`runin` (`start` = 15, `end` = 42, `cat_breaks`, `deciles`, `epv_min`), the onset / gap-model settings as before. Results saved by an earlier
version of the code (other standardization, other odds-ratio direction, other stored content) are recognised as stale and refitted once.

## Another compliance group and its comparison (`pipeline/run_band.R`)

`run_band.R` (notebook cells) repeats the pipeline for people whose compliance lies in **[50, 60)% in all three streams**, every other filter as in the baseline
selection (>= 21 days and longest gap <= 14 days in each stream, CCI <= 10). The list is made by the `../selection` code from the wear metrics, persons and CCI of the
full-pool selection (`result/selection/20261031`); nothing in the model code or in the selection code is changed. The two lists share no person (the script checks it).

| | >= 70% run | 50-60% run |
|---|---|---|
| list | `result/selection/20261031/final_pid_base.txt` | `result/selection/<NEW_DATE>/final_pid_base50_60.txt` (+ `attrition_base50_60.tsv`, `cci.txt`) |
| exports | `<bucket>/bq_exports/<owner>/20261031/` | `.../<NEW_DATE>/` (a new list needs a new date; `20261104`) |
| day tables | `result/final_missing/L1` | `result/final_missing/L1/5060` |
| results | `result/final_missing/results` | `result/final_missing/results/5060` |
| comparison | `result/final_missing/results/compare_70_vs_5060/{tables,figures}` | |

Stages run on the group: `missingness` and `temporal` (the tables the comparison reads), the **primary** onset model (`onset`, role `primary`), and the sensitivity
analyses (including the wear analyses). The compliance check `selection_check` uses 50% as its lower limit for this run.

**Comparison (Cell 8; reads the saved tables of the two runs, refits nothing).** Figures: `compare_primary_effect` (the primary term per stream, both groups in
colour), `compare_sensitivity` (every sensitivity analysis of both groups), `compare_missing_pct`, `compare_longest_gap`, `compare_gaps_per_person`,
`compare_gap_length_distribution`, `compare_transition`. Statistics are in tables outside the plots: `compare_primary` (estimate, SE, OR and size of each
model), `compare_primary_ratio` (the ratio of the two odds ratios, second group relative to the first, with a 95% interval: the groups are different people, so the
two estimates are independent, log ratio = b2 - b1, SE = sqrt(se1^2 + se2^2); no p-value), `compare_within_sd`, `compare_sensitivity`, `compare_missingness`.
Each estimate is **per +1 SD of its own model** (every run standardizes with its own pooled within-person SD of log(1 + steps)); the two SDs are in
`compare_within_sd` and in the figure caption, so a ratio of 1 means equal odds per own SD, not necessarily per the same number of steps. The missingness figures are
descriptive: the % of days missing differs by design, because compliance defines the groups. The two groups also differ in who they are (age, cohort mix, CCI,
device), so a difference is an association between the groups, not an effect of compliance.

## Weekly rhythm of gap starts and gap ends (stage `weekday_rhythm`, Cell 11g)

Does the watch come off (gap start) or go back on (gap end, return) on particular weekdays; is that a shared calendar routine or a charging cycle (days since return); does it differ between people and cohorts?
Run it with its **own results folder** so that nothing of the other analyses is touched:
`run_missingness_pipeline(day_tables, persons, cfg, file.path(OUT_DIR, "weekday"), stages = "weekday_rhythm", model_streams = c("activity", "sleep", "hr"))`.
Each stream has about 16 fits (`S0-S5` gap starts, `E0, E1, E4` gap ends, and the sensitivity fits); every fit is a checkpointed unit, so a crashed run is simply repeated. Terms, rows, models, hypotheses H1-H4 (Holm within four families), effect sizes,
figures and limits: `STATISTICS.md` 3b; the formulas: `MODEL_EQUATIONS.md` 6. Sleep is aligned to the night it started. **No result exists yet on real data** (simulation checks: `STATISTICS.md` section 13, `pipeline/validate_weekday_rhythm.R`).

Settings (`cfg$weekday_rhythm`; they do not enter the keys of any other analysis; each weekday unit carries them in its own key): `dsr_max = 21` (days since return levels 1-21, 22+), `days_so_far_max = 13`, `holiday_window = 1`,
`n_profile_rows = 2000`, `n_re_draws = 100`, `n_beta_draws = 300`, `n_curve_draws = 100`, `n_person_draws = 10000`, `min_people = 20` (cells with 1-19 people are not shown), `min_cell_days = 200` (cells with fewer days at risk are drawn grey),
`sleep_alignment = "night_start"` (`"as_stored"` switches the shift off).

**Measured run time.** *Simulated data:* one gap-start fit with the correlated personal rhythm (`S1`) on 918,961 rows took 8.4 minutes and 1.3 GB (resident memory of the worker). *Real data, activity stream of the >= 70% group (first run on the Workbench, 24-26 GB free; the session never used more than 3.9 GB):*
all 16 fits and the descriptions finished with no failed unit in about **18 h 20 min** (the log times: `S1` 2 h 52 min, `S0` 30 min, `S2` 52 min, `S3` 4 h 4 min, `S4` 4 h 0 min, `S5` 4 min, the gap-end fits `E1`, `E0`, `E4` 1-2 min each because there are only about as many gap-end rows as missing days,
`S1_holidays` 41 min, `S3_holidays` 54 min, `S1_before_first_gap` 43 min, `S3_before_first_gap` 1 h 10 min, `S1_activity_same_rows` 1 h 31 min, `S1_activity` 54 min). The linear scaling of the simulation (about half an hour for `S1`) was too optimistic: expect several hours for each fit with days since return or
weekday x cohort. The sleep and heart-rate streams were not timed; plan on a day or more each. Every fit is checkpointed, so the stage can be stopped and resumed.

## Resources, memory, checkpoints and resuming

**Settings** (`cfg$resources` in `R/config.R`; they never change a result and are not part of any checkpoint key).
Defaults are for a 4-core / 32 GB runtime:

| setting | default | meaning |
|---|---|---|
| `memory_budget_gb` | 24 | most memory one model fit may use |
| `memory_reserve_gb` | 4 | always left free for the OS, Jupyter and your objects |
| `bytes_per_cell` | 56 | memory per (row x model column) of an lme4 fit, used by the guard |
| `isolate_fits` | TRUE | every heavy fit runs in its own worker process, whose memory is returned when it exits |
| `light_workers` | 3 | workers for light units: models B and C of a stream run side by side (each < 1 GB), and the bootstrap replicates of model C; capped so they fit in free memory |
| `ordinal_sparse` | TRUE | model C: sparse random-effect matrices in `ordinal::clmm` (identical model; logLik equal to 4 decimals, about 10x faster and smaller) |
| `memory_wait_sec` | 20 | the guard waits up to this long for a just-finished worker's memory to be returned by the OS before refusing a model |
| `memory_guard` | TRUE | a model that cannot fit in the free memory is **stopped with a diagnosis** (recorded as failed); nothing is subsampled |
| `stop_on_error` | FALSE | FALSE: a failed model is recorded and the run continues; TRUE: stop at the first failure |

Heavy fits run **one at a time** (`max` parallelism of a fit is 1): running several large fits at once would
multiply peak memory, and a fit is dominated by one single-threaded optimiser anyway.

**How memory is kept low**

* The day store reads **one stream, and only the columns a stage needs** from L1; nothing is read when all
  models of a stream are already finished. `load_inputs()` itself reads nothing big.
* Every stage processes one stream at a time; large tables (`inputs`, the person-day risk set `bundle`)
  are built lazily, only if some model of that stream has to run, and are dropped before the next stream.
* The risk set is built **once per stream**; columns are added by reference (no merge copies) and raw / trailing
  columns are only kept where used.
* Models B and C (Cell 15) never build the risk set: they use one row per gap and one activity predictor.
* Every fit, and every comparison of two fits (cohort interaction, spline, random slope), keeps one
  model in memory at a time and runs in a worker.
* The ETL writes each stream to disk and releases it; the notebook keeps no copy of the day tables.

**Checkpoints.** Every independently runnable model is a *unit* and writes
`result/final/missingness/results/checkpoints/<kind>__<stream>__<model>__<settings-tag>.rds`
(kinds: `onset`, `onset_extra`, `gapB`, `gapC`, `sensitivity`, `runin_onset`, `runin_extra`, `runin_gapB`, `runin_gapC`). A file is reused only if its stored key matches
the run asking: the result-changing settings, the **md5 of the L1 files** it was fitted on, the model
specification and the person sample. It also records warnings, run time, R and package versions. Files are written
under a temporary name and renamed, so a crash cannot leave a half-written file that looks valid. A unit that
fails (error, memory guard, or a fit that reports `failed`) is written to a separate `*.failed.rds` with its error
message, shown in `unit_failures.tsv` / `results_summary.md`, and **retried** on the next call; it never
touches other units. `not_estimable`, `not_converged` and `po_rejected` are results (reported as such) and are reused.
`checkpoint_report(dir, cfg)` lists valid / stale / failed / unreadable files; `prune_checkpoints()` removes the unusable ones.

**After a kernel crash:** run Cell 2 (setup), Cell 10 (opens L1), then the cell you were in. Finished models are
skipped (a stream whose models are all finished is not even read); only missing, failed or outdated ones run. ETL
(Cells 5-9) is not repeated while `L1/` exists. Results of all finished streams are collected into the tables on every call.

## Figures (from the saved result tables; no model is refitted)

`stages = "missing_figures"` redraws only the missingness-summary figures above and writes their statistics tables (about 1-2 minutes at
full size, no models, no checkpoints, nothing else is rewritten). Statistics saved next to the figures (`results/tables/`):

| table | content |
|---|---|
| `missing_pct_by_stratum` | N, mean, SD, median, quartiles, P5, P95 of the person-level % missing by cohort x level of sex / age group / race / CCI group |
| `missing_pct_by_stratum_tests` | Kruskal-Wallis (>= 3 groups) or Wilcoxon (2 groups) between the levels within each cohort, and between cohorts within each level; statistic, df, p, Holm p, epsilon-squared / Cliff's delta |
| `missing_pct_by_stratum_models` | quasi-binomial odds ratios of a day being missing, level vs reference, within each cohort (95% CI, dispersion, omnibus F test and its Holm p) |
| `missing_trend_tests` | trend of the missing share with study week (per 4 weeks), calendar year and age (per 10 years): odds ratio per unit with 95% CI, robust z test, per cohort (Holm) and for all cohorts |
| `missing_trend_heterogeneity` | Wald chi-square test that the trend differs between cohorts |
| `missing_weekday_tests`, `missing_weekday_odds_ratios` | Wald chi-square tests of the weekday effect (per cohort, all cohorts) and of weekday x cohort; odds ratios of each weekday vs Sunday per cohort |

The trend and weekday models are pooled logistic regressions of the person x week / year / age / weekday counts with person-clustered
(sandwich) standard errors; the trend is linear on the logit scale. Estimate figures from the models (forest plots, predicted probabilities) draw one dot +/- 1 SD of the estimate (its standard error; `cfg$plot$interval_style = "ci95"` draws the 95% interval, `"box"` one box per estimate: line = estimate, box = central 50% interval, whiskers = 95% interval; the uncertainty of the estimate, not the spread of the data). The cohort transition figure keeps its own setting, `cfg$plot$transition_style` (default `"box"`). Violin figures with
more than 12 violins in total (violins per panel x panels) put their panels in 2 rows (13-24) or 3 rows (more) with a matching image size.

Which figure belongs to which analysis: see the output map above (`figure_index.tsv` has it per figure). Short readings of the older figures:
`onset_primary_forest` = OR of a gap tomorrow per +1 SD higher than usual for each predictor of the primary model; `onset_deciles_<stream>` = observed
(model-free) % of days followed by a gap by decile of the within-person steps of the last full day (d - 1), one line per cohort; `risk_curve` = population-average % chance of a gap tomorrow at
-1 SD / 0 (usual) / +1 SD; `cohort_forest` = the primary effect per cohort with the LRT for a cohort difference; `forest_gap_definition` = the primary effect for gaps >= 1, 2, 3, 7 days;
`gap_dose_response` = model B estimates at gap lengths 1, 2, 3, 5, 13 days, unadjusted vs adjusted; `post_gap_activity_observed_<stream>` = observed steps after gaps by day since the return,
with the model B estimate drawn apart; `gap_length_probs` = model C predicted probability of each length group at -1 SD / 0 / +1 SD; `collinearity_gvif_<model>_<stream>` = VIF-equivalent of every term of ONE model (one figure per model: `onset_primary`, `onset_full_period`, `model_B`, `model_C`).

Odds ratios are never presented as risk ratios, the curves and probabilities are marginal (population-average,
stated in each caption), and all associations are non-causal.

## Measured on simulated data of the real size (9,312 people, 3.5M person-days per stream)

Sandbox with 4 cores and 16 GB. **Not run on the real data**; numbers will differ with the real gap structure.

| step | before | now |
|---|---|---|
| session memory with the inputs open (Cell 10) | 1.0 GB (all streams read; 1.7 GB more while the ETL object stays in the notebook) | 0.11 GB (a stream is read when a stage needs it) |
| session memory after building the activity bundle | 3.4 GB | 1.75 GB (sleep 2.8 GB, hr 2.7 GB; each released before the next stream) |
| models B and C, activity stream (86,017 gaps) | not finished after 47 min (CPU-bound `clmm`) | 262 s (sleep 171 s, hr 403 s), 6 units side by side, session 0.35 GB, each worker about 0.5 GB |
| rerun of Cell 15 when everything is finished | refits | 2 s, 18 units skipped |
| primary onset model, activity (3.56M rows) | 5 min | 3 min 46 s, OR 1.649 (CI 1.63-1.67) recovered |
| cohort-interaction test (two sequential fits) | two fits alive at once | peak 13 GB in total (11 GB worker + 1.8 GB session) |
| one cohort-specific fit (largest cohort) | - | peak 6.2 GB in total |

Not measured at full size: the primary / model-D fits of the sleep and heart-rate streams (the memory guard
estimates about 10 GB each, from rows x columns x 56 bytes; the measurement for activity was 9-11 GB, so the
estimate is about right for the heaviest fit), the 45 exploratory fits and the sensitivity suite. They run
one at a time in a worker and are each checkpointed, but their total time (hours) was not timed end to end.

## Known limitations

* **Wear of the last day.** The only wear measure is the heart-rate minutes recorded on the calendar day (a lower bound of the time worn; the activity table has no wear
  time of its own). A day with no heart-rate record has unknown wear and is excluded from the wear analyses, never counted as short. The wear analyses restrict the rows to
  full-wear last days, which removes the gaps that begin with a partial day: they estimate the association among gaps that start after a full-wear day, a different question
  from the primary one, not a correction of it (STATISTICS.md 6.2).

* The primary standardization uses only earlier observed days (>= 7), so the first week of every person is not in the models and `w`, `b` are noisier early in follow-up;
  the full-period standardization (look-ahead) is the sensitivity analysis, and the trailing 14-day baseline a second one. Keep them in view when interpreting the primary result.
* The run-in measure comes from days BEFORE the modelled days, but it is noisy (28 days) and overlaps with `any_prior_gap` by construction (a gap in the run-in); estimates describe
  association with early wear, not a causal effect of it. Whole-period coverage and strata of it are not used: they contain the outcome.
* Collinearity diagnostics are reported, not acted on: no predictor is removed because of a high VIF; the between-person terms of correlated behaviours are the usual suspects.
* The primary and model-D fits need about 10 GB each at full size. On a machine with less free memory the guard
  stops them with an explanation rather than subsampling; `cfg$model$sample_persons` is the documented, explicit
  alternative (same model, wider intervals, recorded in the checkpoint key).
* The weekly-rhythm fits have a personal weekly rhythm (random slopes), are slower and larger than the primary model, and the personal-rhythm quantities are model based (normal random effects) with BLUPs shrunk toward the shared rhythm;
  the pattern was first seen in the >= 70% compliance group, so no independent confirmation exists yet (planned: the 50-60% group).
* `ordinal::clmm` is single-threaded and is the slowest model (minutes per model with sparse matrices, `ordinal_sparse`);
  when its random-intercept variance is on the boundary the cluster-bootstrap fallback of the model without the random
  effect is used (status `ok_bootstrap_se`).
* A singular lme4 fit (person variance ~ 0) is reported as `singular = TRUE` with its note, and counts as converged.
* Models B and C are conditional on a gap having occurred and exploratory; their curves and probabilities are
  population-average (covariates at their analysed mix, random intercept integrated out), not predictions for one person.

## People in more than one cohort

`final_pid_base.txt` has one row per membership, so a person who is in two cohorts (e.g. insomnia and
diabetes) has two rows (9,328 rows, 8,939 people). The exports are per cohort and are not affected. In the current list 386 people are in more than one
cohort (Diabetes + Insomnia 296, Diabetes + Prostate 44, Insomnia + Prostate 43, all three 3); 8,553 people remain under the default policy. `run_etl()` resolves the
overlap once so that every person counts in exactly one cohort (the cohort tests assume
independent groups): `cfg$overlap$policy` is `"exclude"` (default; overlapping people are left
out), `"priority"` (kept in the first of their cohorts in `cfg$overlap$priority`) or `"error"`.
The number of people per overlap combination is written to `L1/cohort_overlap.tsv`.

## Definitions

* **Window**: first to last record, per person and stream. Only interior gaps are visible; leading
  and trailing missingness and permanent drop-out are not. All results mean "interruptions that
  were followed by resumption".
* **Gap**: >= 1 consecutive missing days. **Gap onset on day t**: day t-1 observed, day t missing.
* **Risk set** (one row per observed day except the last): outcome = a gap starts the next day.
* **Prior-day predictors**: `w_` = how many within-person SDs the day was above (+) the person's own usual, `b_` = how many SDs the person's usual is above the average person; the usual is the mean of
  the person's observed days **before** the day (pre-gap standardization). Odds ratios are *per +1 SD higher than the person's usual* (OR = exp(coefficient) for every term).
* **Run-in**: study days `cfg$runin$start`..`cfg$runin$end` (15-42) of a person's follow-up in a stream (days 1-14 are skipped: the novelty effect of a new device); **missing days in the run-in** = 28 - observed days among them; only days after day 42 are modelled; the models use `log2(1 + missing days)`, centred.

## Analyses and their role

| | Question | Role |
|---|---|---|
| Primary | Is a higher-than-usual prior day (steps) followed by lower odds that a gap starts? Joint adjusted model, gap >= 1 day | the only pre-specified test (`cfg$primary_term`), Holm across streams |
| A | Gap history (length of last gap, days since return) as adjustment | nuisance adjustment; first-gap-only check |
| B | Is activity on the first full day back lower after longer gaps? | MNAR *plausibility* |
| C | Does low prior-day activity predict a longer gap? (ordinal) | exploratory; only if proportional odds holds |
| D | Same as primary with gaps >= 2, 3, 7 days | main secondary analysis |
| Sensitivity | **full-observation-period standardization**, trailing baseline, partial wear, **last day with >= 16 / >= 20 / >= 22 h of heart-rate wear** (predictors of day d), first gap only, main-sleep-only, conditional logit, cluster-robust SEs, spline, random slope | same estimand, other assumptions |
| Run-in | primary onset model, model B and model C on the days / gaps after the run-in, without and with the run-in terms and their interaction | H1 wear trait -> later gaps; H2 effect modification of steps; H3 robustness |
| Diagnostics | generalized VIF of every term of every model; joint Wald tests of factor / spline terms | reported, no predictor is removed |

Gap length is **never** a predictor of whether a gap starts (it is only known after the gap).

## Read the results carefully

* Gaps preceded by a different activity level than usual show that missingness is **informative** (not completely at
  random). They do not prove MNAR: that depends on the unobserved activity *during* the gap.
* The data were filtered upstream (>= 70% observed, longest gap <= 14 days, >= 21 days, in every one of the three streams; CCI <= 10). Gap length
  is therefore limited to 1-14 days, long-gap effects are under-represented, and everything is
  conditional on being a good wearer. `selection_check.tsv` verifies the rule is met.
* Adjusting for CCI together with age and cohort over-adjusts (CCI contains age points and the
  index disease); CCI is used as a stratifier only. Prostate is all male and ~96% over 65.
* The default mixed-model engine is `lme4::glmer(nAGQ = 0)`: about 10x faster than glmmTMB and
  100x faster than `nAGQ = 1`, with fixed effects equal to 3 decimals in this setting. Test
  `verify_engine_agreement()` compares it with glmmTMB.

## Old function -> new

| old | new |
|---|---|
| `wear_compliance_func`, `miss_gap_func`, `filter_pid_func` | `get_gap_runs`, `summarise_person_gaps`, `cohort_missingness_table`, plots in `plots.R` |
| `read_bq_export_from_workspace_bucket` (16 copies) | `read_export_files` + `export_path` |
| `sleep_sum_main_nap_func` | `standardize_sleep` |
| `distinct(person_id, date)` on heart rate | `standardize_hr` (all four zones, deterministic) |
| `ac_*_plot_func` | `plot_metric_violin`, `plot_single_metric_violin`, `plot_metric_by_stratum` |
