# Missingness pipeline (Fitbit activity, sleep, heart rate)

Day-level missingness for the four cohorts (Insomnia, Diabetes, Prostate, Healthy): summary
statistics and cohort comparisons, temporal structure, day-of-week clustering, a
mixed-effects model of **gap onset ~ prior-day activity** (predictors standardized only with data observed before each
day), gap-level models of activity after a gap and of gap length, the same models by **wear consistency**, sensitivity analyses
and collinearity diagnostics. (Disease-label prediction and the mixture model are not part of this.)

**Where to read what:** this file = layout, how to run, the analysis workflow and a map of every output (tables and figures) to its analysis;
`STATISTICS.md` = the statistical methods (definitions, standardization and its reference distributions, models, effect sizes, tests,
collinearity, wear-consistency strata, validation); `MODEL_EQUATIONS.md` = the exact formulas, predictors, outcomes and contrasts per stream.

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
R/wear_consistency.R      wear-consistency classification, descriptive statistics, models within each stratum, comparison of strata
R/plots.R                 ggplot2 figures in the style of the earlier notebooks + figure captions
R/plots_effects.R         comparison forests (primary vs sensitivity / stratum), collinearity, observed post-gap, wear-consistency figures
R/report.R                tables, captions, results summary, run metadata
R/pipeline.R              run_etl() and run_missingness_pipeline() (orchestration only)
R/simulate.R              seeded fake data with a planted effect (tests / demo)
pipeline/notebook_cells.R the notebook, Cells 1-19 (copy each "# ==== CELL n" block into a cell)
pipeline/validate_simulation.R  validation of the primary estimate on simulated data with a known effect (bias, SE calibration, coverage)
pipeline/run_exports.R    step 0, only if the CSV exports are missing (e.g. after a migration)
pipeline/run_etl.R        run once per export     |  pipeline/run_missingness.R   the analysis
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
In the notebook use `pipeline/notebook_cells.R` (Cells 1-19). Inputs besides BigQuery: `final_pid.txt`
(person -> cohort) and `pre_final_cci.txt`. Use the same export date in `run_exports.R` and `run_etl.R`.
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
than the person's usual** (model B: change in activity in SD units); reference levels Healthy (cohort), Female, Sunday, the calendar year with the most gap starts (years with fewer than 10 gap starts are pooled with a neighbouring year, e.g. `2010-2012`).

**Data sources and sample flow.** BigQuery exports (activity, sleep, heart rate, meta; Cells 6-7) -> day tables `L1/` (one row per observed person-day; Cell 8) ->
people per cohort after the overlap rule (Cell 5; `cohort_overlap.tsv`) -> the upstream selection (>= 70% observed, longest gap <= 14 days, >= 21 days; checked in
`selection_check`) -> models (the rows, people, events and gaps of every model are in `onset_models`, `model_B_post_gap_activity`, `model_C_gap_length`; days without
7 earlier observed days have no z-score and are not in a model) -> wear strata (`wear_strata_quality`).

| step | analysis (cell) | tables (`results/tables/`) | figures (`results/figures/`) |
|---|---|---|---|
| 1 | **Descriptive statistics** of missingness and of the metric values (Cell 11, 11b) | `cohort_missingness`, `cohort_tests_overall`, `cohort_tests_pairwise`, `missing_fraction_models`, `missing_pct_by_stratum*`, `metric_summary*`, `selection_check` | `heatmap_<stream>`, `missing_pct_by_cohort`, `longest_gap_by_cohort`, `gap_length_distribution`, `gap_ecdf`, `missing_pct_by_sex` / `_age_group` / `_race_plot` / `_cci_group`, `metrics_<stream>_by_cohort`, `metrics_<stream>_by_<stratum>` |
| 2 | **Temporal structure** (Cell 11, 11b) | `transition`, `study_week_rates`, `calendar_year_rates`, `weekday_rates`, `weekday_tests`, `weekday_odds_ratios`, `missing_trend_tests`, `missing_trend_heterogeneity`, `missing_weekday_*`, `age_rates` | `transition`, `missing_by_study_week`, `missing_by_year`, `missing_by_weekday`, `missing_by_age` |
| 1 | **Descriptive: person-level % of days missing against every model variable** (Cell 11e) | `outcome_variable_tests`, `outcome_variable_contrasts`, `outcome_variable_models`, `outcome_variable_summary`, `outcome_variable_bins` | `outcome_by_variable_categorical`, `outcome_by_variable_continuous_<stream>` |
| 3 | **Data quality: watch version** (Cells 6c, 8b, 11d; person level, no model) | `device_cohort_mix`, `device_missing_by_cohort`, `device_tests`, `device_missing_vs_reference`, `device_missing_model`, `device_composition_by_cohort`, `device_composition_residuals`, `device_version_by_cohort`, `device_switching_by_cohort` | `device_version_missing_and_composition`, `device_cohort_mix_by_version` |
| 3 | **Data quality: wear consistency** (Cell 11c; classification rule written before any stratified model) | `wear_consistency_classification`, `wear_consistency_descriptive`, `wear_consistency_composition`, `wear_consistency_by_person` | `wear_consistency_distribution`, `wear_consistency_composition` |
| 4a | **Primary model: gap onset** (primary, model A inside it, model D, exploratory grid; Cells 12-14, 17) | `onset_coefficients` (every term), `onset_models`, `primary_test`, `onset_term_tests`, `vif`, `cohort_interaction`, `cohort_specific`, `risk_unadjusted`, `onset_deciles` | `onset_primary_forest`, `onset_forest_year_<stream>` + `onset_forest_other_<stream>` (all terms), `onset_deciles_<stream>`, `risk_curve`, `cohort_forest`, `forest_gap_definition`, `gap_definition_all_terms_other_<stream>` + `_year_`, `exploratory_grid_effects_<stream>` |
| 4b | **Gap-level model B: activity after a gap** (Cell 15) | `model_B_post_gap_activity`, `model_B_all_terms`, `model_B_curve`, `model_B_observed`, `model_B_observed_trajectory` | `gap_dose_response`, `gapB_forest_<stream>`, `post_gap_activity_observed_<stream>` |
| 4c | **Gap-level model C: gap length** (Cell 15) | `model_C_gap_length`, `model_C_all_terms`, `model_C_predicted_probs` | `gap_length_probs`, `gapC_forest_<stream>` |
| 4d | **MVPA before the gap as a covariate of models B and C** (Cell 15; units `B_back2_adj_mvpa`, `B_back2_adj_pre_mvpa`, `C_adj_mvpa`) | `gap_mvpa_effects` (MVPA effect with 95% CI, the main terms with and without MVPA, n, VIF); the MVPA fits are also in `model_B_post_gap_activity`, `model_C_gap_length`, `model_B_all_terms`, `model_C_all_terms` | `gapB_mvpa_effect_<stream>`, `gapC_mvpa_effect_<stream>` |
| 5 | **Stratified analyses by wear consistency** (Cell 18) | `wear_strata_onset_coefficients`, `wear_strata_onset_models`, `wear_strata_gapB` / `_gapC` (+ `_all_terms`), `wear_strata_interaction`, `wear_strata_comparison`, `wear_strata_quality` | `wear_strata_primary_effect`, `wear_strata_onset_all_terms_other_<stream>` + `_year_`, `wear_strata_gapB_all_terms_<stream>`, `wear_strata_gapC_all_terms_<stream>` |
| 6 | **Sensitivity analyses** incl. full-period standardization (Cell 16; for models B and C the `*_fullperiod` units of Cell 15) | `sensitivity`, `sensitivity_all_terms`, `sensitivity_collinearity` | `forest_sensitivity`, `sensitivity_forest_<stream>`, `sensitivity_all_terms_other_<stream>` + `_year_` |
| 7 | **Model diagnostics and integrated summary** (written by every model cell) | `model_collinearity` (GVIF of every model), `model_term_tests`, `onset_collinearity`, `gap_collinearity`, `wear_strata_onset_collinearity`, `unit_failures`, `checkpoint_status`, `unit_run_status`; `results_summary.md` | `collinearity_gvif_<stream>` |

`figures/figure_index.tsv` lists every figure present with its step, analysis, model and caption; `figures/figure_captions.md` has the captions;
`tables/table_index.tsv` lists every table with its step.
The all-terms figures come in two parts: `..._year_...` = the calendar-year terms (each year, or pooled years, against the year with the most gap starts), `..._other_...` = every other term
(within-person, between-person, cohort, age, sex, weekday, gap history); the study-time spline and the intercept are not drawn (joint tests: `onset_term_tests`,
`model_term_tests`). In all comparison figures the first (dark) colour is the primary analysis.

**Settings that change results** (`cfg`, part of every checkpoint key): `standardization` (`method = "pre_gap"`, `baseline_min_n = 7`), `reference_cohort`,
`wear` (`method = "median"` or `"fixed"`, `threshold`, `onset_roles`, `min_persons_stratum`, `epv_min`), the onset / gap-model settings as before. Results saved by an earlier
version of the code (other standardization, other odds-ratio direction, other stored content) are recognised as stale and refitted once.

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
(kinds: `onset`, `onset_extra`, `gapB`, `gapC`, `sensitivity`, `wear_onset`, `wear_extra`, `wear_gapB`, `wear_gapC`). A file is reused only if its stored key matches
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
(model-free) % of days followed by a gap by decile of yesterday's within-person steps, one line per cohort; `risk_curve` = population-average % chance of a gap tomorrow at
-1 SD / 0 (usual) / +1 SD; `cohort_forest` = the primary effect per cohort with the LRT for a cohort difference; `forest_gap_definition` = the primary effect for gaps >= 1, 2, 3, 7 days;
`gap_dose_response` = model B estimates at gap lengths 1, 2, 3, 5, 13 days, unadjusted vs adjusted; `post_gap_activity_observed_<stream>` = observed steps after gaps by day since the return,
with the model B estimate drawn apart; `gap_length_probs` = model C predicted probability of each length group at -1 SD / 0 / +1 SD; `collinearity_gvif_<stream>` = VIF-equivalent of every term.

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

* The primary standardization uses only earlier observed days (>= 7), so the first week of every person is not in the models and `w`, `b` are noisier early in follow-up;
  the full-period standardization (look-ahead) is the sensitivity analysis, and the trailing 14-day baseline a second one. Keep them in view when interpreting the primary result.
* Wear-consistency strata are defined from the same observed days as the outcome; the stratified estimates describe effect modification by wear pattern, not a causal effect of
  consistent wear, and the strata differ in composition (see the composition table).
* Collinearity diagnostics are reported, not acted on: no predictor is removed because of a high VIF; the between-person terms of correlated behaviours are the usual suspects.
* The primary and model-D fits need about 10 GB each at full size. On a machine with less free memory the guard
  stops them with an explanation rather than subsampling; `cfg$model$sample_persons` is the documented, explicit
  alternative (same model, wider intervals, recorded in the checkpoint key).
* `ordinal::clmm` is single-threaded and is the slowest model (minutes per model with sparse matrices, `ordinal_sparse`);
  when its random-intercept variance is on the boundary the cluster-bootstrap fallback of the model without the random
  effect is used (status `ok_bootstrap_se`).
* A singular lme4 fit (person variance ~ 0) is reported as `singular = TRUE` with its note, and counts as converged.
* Models B and C are conditional on a gap having occurred and exploratory; their curves and probabilities are
  population-average (covariates at their analysed mix, random intercept integrated out), not predictions for one person.

## People in more than one cohort

`final_pid.txt` has one row per membership, so a person who is in two cohorts (e.g. insomnia and
diabetes) has two rows. The exports are per cohort and are not affected. `run_etl()` resolves the
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
* **Wear consistency**: observed days / days between the first and last record, per stream; Consistent = at or above the stream's median (`cfg$wear`).

## Analyses and their role

| | Question | Role |
|---|---|---|
| Primary | Is a higher-than-usual prior day (steps) followed by lower odds that a gap starts? Joint adjusted model, gap >= 1 day | the only pre-specified test (`cfg$primary_term`), Holm across streams |
| A | Gap history (length of last gap, days since return) as adjustment | nuisance adjustment; first-gap-only check |
| B | Is activity on the first full day back lower after longer gaps? | MNAR *plausibility* |
| C | Does low prior-day activity predict a longer gap? (ordinal) | exploratory; only if proportional odds holds |
| D | Same as primary with gaps >= 2, 3, 7 days | main secondary analysis |
| Sensitivity | **full-observation-period standardization**, trailing baseline, partial wear, first gap only, main-sleep-only, conditional logit, cluster-robust SEs, spline, random slope | same estimand, other assumptions |
| Stratified | primary onset model, model B and model C within consistent / inconsistent wear; strata compared with an interaction test | effect modification by wear pattern; descriptive of strata first |
| Diagnostics | generalized VIF of every term of every model; joint Wald tests of factor / spline terms | reported, no predictor is removed |

Gap length is **never** a predictor of whether a gap starts (it is only known after the gap).

## Read the results carefully

* Gaps preceded by a different activity level than usual show that missingness is **informative** (not completely at
  random). They do not prove MNAR: that depends on the unobserved activity *during* the gap.
* The data were filtered upstream (>= 70% observed, longest gap <= 14 days, >= 21 days). Gap length
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
