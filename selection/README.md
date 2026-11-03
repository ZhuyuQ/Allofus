# Cohort selection (full pool -> final person list)

Builds the list of people (`final_pid_<tag>.txt`, columns `cohort`, `person_id`) that the missingness pipeline analyses, directly from the
All of Us Curated Data Repository, for any set of rules. Defaults = the rule of the earlier list; any rule can be changed on its own while the
others stay at their baseline (for example compliance 50-60% with everything else unchanged). Nothing in `missingness/` (`final_missing/` on the
Workbench) is modified: that folder is only read (its `age_on()` and `log_info()`, the export code of Cell 9, `race_group_from_name()` of Cell 10), and
`functions/cci_function.R` is sourced unchanged.

CDR **C2024Q3R9** (`wb-silky-artichoke-2408.C2024Q3R9`), bucket `gs://rw-migration-aou-rw-6f1dba72`. Both are resolved in Cell 1 with `wb resource resolve`.

## Files

| file | content |
|---|---|
| `selection_sql.R` | `SEL_COHORTS` (cohort definitions), `SEL_STREAMS`, `SEL_FLAGS`, the SQL string builders and `sel_fetch()`, the only function that talks to BigQuery |
| `selection_core.R` | the rules in plain R: wear metrics, age at the first record, CCI, `sel_select()`, comparison, writing, verification |
| `summary_tables.R` | Table 1, the Fitbit characteristics table and the Venn diagram of a selection |
| `selection_cells.R` / `selection_pipeline.ipynb` | the notebook cells in order (Cells 1-10; one `# ==== CELL n` block = one cell; the notebook holds the same cells) |
| `tests/test_selection.R`, `tests/test_summary.R` | tests (no BigQuery); run from `selection/` with `Rscript tests/test_selection.R` and `Rscript tests/test_summary.R` (the second writes `tests/venn_example.png`) |

Upload the repository folder `Allofus/` (`final_missing/`, `functions/`, `selection/`) to `ROOT` of Cell 2.

## The rule

A person is selected when **all** of the following hold.

1. **Cohort membership** (pool, below).
2. In **each of the three streams** (activity, sleep, heart-rate summary), with `n` = number of distinct dates with a record, `first` / `last` = first / last
   date, `F = last - first + 1` (follow-up span, inclusive):
   * record days `n >= min_days` (21);
   * longest interior gap `<= gap_max` (14 days): the longest run of calendar days without a record **between** two records (0 if none);
   * compliance `100 n / F` in `[comp_min, comp_max)` (70, Inf). It is compared exactly, `100 n >= comp_min F` and `100 n < comp_max F`, never through a
     rounded ratio. The upper bound is what a band such as 50-60% uses (`[50, 60)`, so 60% itself is in the next band).
3. **CCI** `<= cci_max` (10), computed with the age at the person's first record (below).

A person without data in one of the streams, or without a CCI score, is not selected. The attrition after each step is written for every run.

## Data used

* **Streams** (`SEL_STREAMS`): `activity_summary` (`date`), `sleep_daily_summary` (`sleep_date`), `heart_rate_summary` (`date`). A day is observed when a
  record exists; dates are counted once (`DISTINCT`). The same definition as `filter_pid_func` of the earlier notebooks (the test compares the two on
  random and boundary cases) and as the missingness pipeline.
* **Wear metrics** are computed in BigQuery per cohort and stream (`LAG` / `DATE_DIFF` over the sorted distinct dates), saved per cohort x stream in
  `OUT_DIR/units/wear_<cohort>_<stream>.tsv` (a finished file is read, not queried again) and verified in Cell 5 against the metrics recomputed in R
  (`wear_from_dates`) from the raw dates of a random sample of people per stream (the cell stops on any difference).

### Pool of each cohort (`SEL_COHORTS`, `sel_sql_pool`)

All cohorts need the Fitbit flags of `cb_search_person`: `has_fitbit`, `has_fitbit_activity_summary`, `has_fitbit_plus_device_information`,
`has_fitbit_heart_rate_summary`, `has_fitbit_heart_rate_level`, `has_fitbit_steps_intraday`, `has_fitbit_sleep_daily_summary`,
`has_fitbit_sleep_level` (all = 1). Concepts are expanded to their descendants as the dataset builder does (`cb_criteria` rank-1 parents, standard
concepts, events in `cb_search_all_events`).

| cohort | additional condition |
|---|---|
| Diabetes | an event in 201820 or a descendant |
| Insomnia | an event in 436962, 42538607, 37110488, 4102985 or 439013 (or descendants) |
| Prostate | also `has_fitbit_device_information` = 1; an event in 4163261 or 200970 (or descendants); sex at birth not 45878463 (Female) |
| Healthy | none of the 26 exclusion concepts (or descendants): 135930, 4101796, 4117930, 255919, 444112, 432453, 4024567, 4274025, 376106, 4011630, 4213101, 4038502, 4093991, 40480457, 4041284, 4042141, 43021226, 4041283, 4227253, 4023995, 4102111, 138239, 4293175, 4212577, 4302537, 4021667 |

A person can be in several disease cohorts; the list keeps **one row per membership** (the missingness pipeline resolves the overlap itself, policy `exclude`).

### Age

Whole years (`age_on()`) at the person's **earliest record over the three streams** (`first_date`), from the date of birth in `person`. The age at the date
a script is run is not used: in a longitudinal record the age must not depend on the calendar date of the analysis, and the same age is used by the
CCI here and by every model of the missingness pipeline (`persons.R` takes the earliest first date of the three day tables).

### CCI

`cci_calc_func()` of `functions/cci_function.R`, unchanged. The 785 SNOMED codes of its lookup are read from the function body (`sel_cci_codes()`), so the
two cannot differ; the conditions come from `condition_occurrence` joined to `concept` for those codes, per cohort pool. People without condition rows
get the age points only. The score depends on the cohort because the conditions are queried per pool; the file `cci.txt` has the format of
`pre_final_cci.txt` and is the CCI file of the missingness ETL.

## Run (Cells of `selection_cells.R`)

| cell | what | writes (in `OUT_DIR`) |
|---|---|---|
| 1-2 | CDR, bucket, paths, sources; run both again after a kernel restart (the CDR is an environment variable of the R session) | - |
| 3 | pool and date of birth of every cohort | - |
| 4 | wear metrics | `units/wear_*.tsv`, `wear_metrics.tsv` |
| 5 | verification of the SQL metrics against the raw dates | - |
| 6 | age at the first record, conditions, CCI | `persons.tsv`, `cci.txt` |
| 7 | selection with the default rules, and its comparison with the earlier list | `final_pid_base.txt`, `attrition_base.tsv` |
| 8 | another rule set, for example `sel_write(sel_select(wear, cci, pool, comp_min = 50, comp_max = 60), OUT_DIR, "comp50-60")` | `final_pid_comp50-60.txt`, `attrition_comp50-60.tsv` |
| 9 | the 16 BigQuery exports of the selected people (the export code of the missingness pipeline) | `<bucket>/bq_exports/<owner>/<EXPORT_DATE>/` |
| 10 | summary tables and the Venn diagram | `tables/`, `figures/` |

`OUT_DIR = <BASE>/result/selection/<EXPORT_DATE>` (`EXPORT_DATE = 20261031`, a fixed label of this selection, not the calendar date). Arguments of
`sel_select()`: `min_days = 21`, `gap_max = 14`, `comp_min = 70`, `comp_max = Inf`, `cci_max = 10`, `streams`; only the argument named changes. Cells 3-7 are needed
once; Cell 8 for another band reuses `wear` and `cci`. The export of Cell 9 skips folders that exist, so a failed run can be repeated; **a different
person list needs a new `EXPORT_DATE`** (otherwise the old exports are reused), and the missingness ETL then runs with the same date, `FINAL_PID` and `CCI_FILE`.

## Summary tables and Venn diagram (Cell 10, `summary_tables.R`)

Computed from the selection's own files (`final_pid_<tag>.txt`, `persons.tsv`, `cci.txt`, `wear_metrics.tsv`) plus one query of the person table (race,
ethnicity, sex at birth). Two variants are written: `all` (every cohort membership, N = rows of the list) and `excl_overlap` (people in more than one
cohort removed, the analysed sample).

* `tables/table1_demographics_<tag>_<variant>.tsv`: N; age group at the first record (18-44, 45-64, >=65); ethnicity (Hispanic or Latino, not, others); race
  (grouped by name with `race_group_from_name()`); sex at birth (Female, Male, others); CCI group (0, 1-2, 3-4, >=5). "n (%)", "-" for 0, "<1%" below one percent.
* `tables/table_fitbit_<tag>_<variant>.tsv`: median (IQR) of follow-up days, record days, missing rate (% = 100 - compliance) and longest gap, per stream.
* `tables/venn_counts_<tag>.tsv`, `figures/venn_<tag>.png` / `.pdf`: three equal circles (Diabetes, Insomnia, Prostate cancer); the regions are **not**
  proportional to the counts. Healthy is shown in a separate box because nobody in Healthy belongs to another cohort.

## Result of the baseline selection (`final_pid_base.txt`)

| cohort | selected |
|---|---|
| Insomnia | 1,388 |
| Diabetes | 1,306 |
| Prostate | 253 |
| Healthy | 6,381 |
| rows | 9,328 |

8,939 people, 386 of them in more than one cohort (Diabetes + Insomnia 296, Diabetes + Prostate 44, Insomnia + Prostate 43, all three 3); 8,553 remain when
the overlap is excluded. Against the earlier list (`final_pid.txt`, CDR R8): Healthy is identical; every earlier person passes the new wear rules; the
differences come only from the CCI recomputed with the age at the first record (9, 5 and 3 people added in the three disease cohorts and one Diabetes
person lost).

## Checks and limits

* The tests cover: the wear metrics against the earlier `filter_pid_func` (random cases and the boundaries 21 days, gap 14 / 15, 70%), exact compliance
  bands `[lo, hi)`, the all-streams rule, age at the first record, the CCI with the new age, the SQL strings (cohort definitions, flags, no partial-name
  bugs) and the summary tables. They run on simulated data; the SQL itself is verified only by Cell 5 on the real data.
* The pool depends on the CDR version; the selection is reproducible for C2024Q3R9 only.
* Compliance uses the span first to last record in **each** stream, so leading and trailing missingness is not counted, as in the earlier rule.
* The selection is conditional on being a good wearer in all three streams; this limits gap lengths to 1-14 days downstream.
