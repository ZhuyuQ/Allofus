# ==============================================================================
# config.R
# ------------------------------------------------------------------------------
# ONE place for every setting that affects a result. Nothing below is hard-coded
# elsewhere, so a run is reproducible from `cfg` + the input day tables.
#
# Usage:   cfg <- make_config()                      # defaults
#          cfg <- make_config(list(seed = 1L))       # override anything
# ==============================================================================

make_config <- function(overrides = list()) {
  cfg <- list(

    # ---- reproducibility -----------------------------------------------------
    seed = 20260511L,

    # ---- cohorts -------------------------------------------------------------
    # `cohorts` is the plotting / factor order used in the earlier notebooks.
    # `cohort_prefix` maps the cohort label to the prefix used in the bucket
    # paths (Prostate is stored as `pca_*`).
    cohorts       = c("Insomnia", "Diabetes", "Healthy", "Prostate"),
    reference_cohort = "Healthy",   # the reference level of `cohort` in every model (other cohorts are compared with it)
    cohort_prefix = c(Insomnia = "insomnia", Diabetes = "diabetes",
                      Healthy  = "healthy",  Prostate = "pca"),

    # ---- where the raw BigQuery exports live (used only by the ETL) ----------
    export = list(
      bucket  = Sys.getenv("WORKSPACE_BUCKET"),
      owner   = Sys.getenv("OWNER_EMAIL"),
      date    = "20260511",
      # BigQuery export only: the CDR dataset as "project.dataset" and the billing project.
      # They default to the notebook variables; set them here if those are blank or malformed.
      cdr     = Sys.getenv("WORKSPACE_CDR"),
      billing = Sys.getenv("GOOGLE_PROJECT")
    ),

    # ---- people who belong to more than one cohort ---------------------------
    # "exclude": leave them out (default; groups stay independent) | "priority": keep them in
    # the first cohort of `priority` they belong to | "error": stop and look at the overlap.
    overlap = list(policy = "exclude",
                   priority = c("Prostate", "Diabetes", "Insomnia", "Healthy")),

    # ---- raw stream registry -------------------------------------------------
    # sharded = TRUE  -> files are  <prefix>_<stem>_*.csv   (BigQuery wildcard)
    # sharded = FALSE -> one file   <prefix>_<stem>.csv
    # col_classes     -> forced column types (the old readers did this per stream)
    streams = list(
      activity = list(stem = "activity", sharded = TRUE,
                      date_col = "date",       col_classes = NULL),
      sleep    = list(stem = "sleep",    sharded = FALSE,
                      date_col = "sleep_date", col_classes = c(is_main_sleep = "character")),
      hr       = list(stem = "hr",       sharded = FALSE,
                      date_col = "date",       col_classes = c(zone_name = "character"))
    ),
    meta_stem = "meta",   # person table export (always sharded)

    # ---- storage of the compact day tables (L1) ------------------------------
    # "auto" = Parquet if the `arrow` package is installed, otherwise RDS.
    storage_format = "auto",

    # ---- selection already applied upstream (filter_pid_func) ----------------
    # We do NOT filter again. These values are only used to CHECK the assumption
    # and to label outputs (see check_selection_rule()).
    selection = list(min_pct_observed = 70, max_gap_days = 14L, min_days = 21L),

    # ---- ETL rules -----------------------------------------------------------
    etl = list(
      na_drop_frac = 0.90,        # drop a column if > 90% of its values are NA
      hr_zones     = c("out_of_range", "fat_burn", "cardio", "peak"),
      hr_elevated  = c("fat_burn", "cardio", "peak"),
      hr_max_min_per_day = 1440   # total HR minutes above this is reported
    ),

    # ---- predictors for the gap-onset models ---------------------------------
    # name -> which stream/column, how to transform it, and what it represents.
    # "coverage" predictors measure how much the watch recorded (wear), not behaviour.
    predictors = list(
      act_steps     = list(stream = "activity", col = "steps",             transform = "log1p",    kind = "behaviour"),
      act_sedentary = list(stream = "activity", col = "sedentary_minutes", transform = "identity", kind = "behaviour"),
      act_mvpa      = list(stream = "activity", col = "mvpa_minutes",      transform = "log1p",    kind = "behaviour"),
      sleep_inbed   = list(stream = "sleep",    col = "minute_in_bed",     transform = "log1p",    kind = "behaviour"),
      sleep_eff     = list(stream = "sleep",    col = "sleep_efficiency",  transform = "identity", kind = "behaviour"),
      hr_cover      = list(stream = "hr",       col = "total_hr_min",      transform = "identity", kind = "coverage"),
      hr_elev       = list(stream = "hr",       col = "elevated_min",      transform = "log1p",    kind = "behaviour")
    ),
    own_predictors = list(
      activity = c("act_steps", "act_sedentary", "act_mvpa"),
      sleep    = c("sleep_inbed", "sleep_eff"),
      hr       = c("hr_cover", "hr_elev")
    ),
    activity_predictors = c("act_steps", "act_sedentary", "act_mvpa"),

    # The ONE pre-specified coefficient that carries the primary hypothesis
    # (within-person prior-day steps, in the joint adjusted >=1-day model).
    primary_term = "w_act_steps",

    # ---- onset-model settings ------------------------------------------------
    onset = list(
      primary_k         = 1L,              # a gap = at least 1 missing day
      secondary_k       = c(2L, 3L, 7L),   # model D: stricter gap definitions
      trailing_n        = 14L,             # trailing baseline = last 14 observed days
      trailing_min_n    = 7L,              # ... needing at least 7 of them
      days_since_cap    = 60L,             # cap for "days since last gap ended"
      # a prior day is "partial wear" if any of these (R expressions on the day
      # table columns) is TRUE -> excluded in the partial-wear sensitivity
      partial_wear = list(activity = "steps < 500 | sedentary_minutes >= 1380",
                          sleep    = "minute_in_bed < 60",
                          hr       = "total_hr_min < 1200")
    ),

    # ---- mixed-model settings (these change RESULTS; resource settings are in `resources`) --------
    model = list(
      # "lme4" with nAGQ = 0 is ~10x faster than glmmTMB and ~100x faster than nAGQ = 1 here,
      # with the same fixed effects to 3 decimals (checked in tests / verify_engine_agreement()).
      engine                  = "lme4",      # "lme4" or "glmmTMB"
      nAGQ                    = 0L,          # only used when engine == "lme4"
      min_events              = 100L,        # below this a model is "not estimable"
      min_persons_with_events = 30L,
      sample_persons          = NULL,        # NULL = everybody. e.g. 2000 -> fit on a seeded random sample
                                             # of people (an explicit choice: it changes the estimate's precision)
      # speed knob: fit only the EXPLORATORY specs on a seeded sample of people (NULL = all people)
      exploratory_sample_persons = NULL
    ),

    # ---- computing resources (these never change a result, so they are not part of any checkpoint key) --
    # Defaults are for 4 CPU cores / 32 GB RAM: one heavy model at a time, ~8 GB kept free for the
    # operating system, Jupyter and the notebook's own objects.
    resources = list(
      memory_budget_gb  = 24,      # most memory (GB) the pipeline may use for one model fit
      memory_reserve_gb = 4,       # always leave at least this much free
      bytes_per_cell    = 56,      # lme4::glmer working memory per (row x model column), measured
      isolate_fits      = TRUE,    # run each heavy fit in its own worker process (memory freed on exit)
      light_workers     = 3L,      # workers for light units (models B and C, bootstrap replicates); capped by free memory
      ordinal_sparse    = TRUE,    # model C: sparse random-effect matrices (same model and estimates, ~10x faster, less memory)
      memory_wait_sec   = 20,      # the guard waits this long for a finished worker's memory to be returned
      memory_guard      = TRUE,    # stop ONE model with advice instead of letting the kernel die
      stop_on_error     = FALSE    # FALSE: a failed model is recorded and the run continues
    ),

    # ---- gap-level models (B and C) -----------------------------------------
    gap_models = list(
      length_bins   = c(1L, 2L, 3L, 4L, 8L, 15L),   # -> "1","2","3","4-7","8-14","15+"
      po_alpha      = 0.05,                         # proportional-odds test level
      min_gaps      = 100L,                         # fewer gaps than this -> "not_estimable"
      min_bin_n     = 20L                           # length groups with fewer gaps are merged downwards
    ),

    # ---- uncertainty ---------------------------------------------------------
    boot = list(B = 200L),   # person-level (cluster) bootstrap replicates

    # ---- checkpoints --------------------------------------------------------------------------------
    # One file per finished model, in <results>/checkpoints (see ?run_unit in utils.R).
    checkpoint = list(dirname = "checkpoints"),

    # ---- predictions drawn for the figures of models B and C ------------------------------------------
    predict = list(n_draws = 300L,          # parameter draws behind the intervals of predicted probabilities
                   n_covariate_rows = 500L, # seeded sample of gaps used to average over the covariate mix
                   gap_curve_max = 14L),    # model B curve is drawn for gap lengths 1 .. this

    # ---- strata used in descriptive plots ------------------------------------
    strata = list(
      age_breaks = c(-Inf, 44, 64, Inf), age_labels = c("18-44", "45-64", ">=65"),
      cci_breaks = c(-Inf, 0, 2, 4, Inf), cci_labels = c("0", "1-2", "3-4", ">=5"),
      age_half_width = 2L   # "missing by age" figure: each age pools the days at age -2 ... +2 years (a 5-year moving window)
    ),
    suppress_n_lt = 0L,      # small-cell suppression OFF (set to 20 to turn on)

    # ---- plotting palettes (copied from the earlier plotting code) -----------
    palette = list(
      cohort    = c(Diabetes = "#E74C3C", Prostate = "#3498DB",
                    Insomnia = "#F39C12", Healthy  = "#27AE60"),
      sex       = c(Female = "#E91E8C", Male = "#3498DB", Other = "#95A5A6"),
      age_group = c("18-44" = "#E91E8C", "45-64" = "#3498DB", ">=65" = "#95A5A6"),
      cci       = c("0" = "#27AE60", "1-2" = "#F39C12", "3-4" = "#E67E22", ">=5" = "#E74C3C"),
      race      = c(White = "#3498DB", `Black or African American` = "#E67E22",
                    Asian = "#27AE60", `More than one population` = "#8E44AD",
                    Other = "#95A5A6")
    ),
    plot = list(dpi = 150,
                interval_style = "se",       # model figures: "se" = dot +/- 1 SD (standard error); "ci95" = dot + 95% interval; "box" = line + 50% box + 95% whiskers
                transition_style = "box")    # the cohort transition-probability figure only (same three choices)
  )
  utils::modifyList(cfg, overrides)
}
