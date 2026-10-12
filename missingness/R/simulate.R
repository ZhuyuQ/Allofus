# ==============================================================================
# simulate.R
# ------------------------------------------------------------------------------
# A small, seeded simulator that produces RAW tables in the same shapes as the
# BigQuery exports (activity, sleep with main/nap rows, heart rate with four zone
# rows per day, meta, final_pid, cci). It is used by the tests and the demo to
# check that the whole pipeline recovers a KNOWN effect:
#
#   the log-odds that a gap starts tomorrow changes by `beta` for every 1 SD of
#   the person's within-person log-steps today   (beta < 0: low steps -> more gaps)
#
# The data are fake. Nothing here reads or writes real participant data.
# ==============================================================================

# Presence matrix (persons x days): TRUE = observed. Gaps start with probability
# plogis(base_logit + person_effect + beta * z[, t-1]) and last 1 + Geometric days.
simulate_presence <- function(z, base_logit, person_effect, beta, gap_p = 0.5, max_gap = 14L) {
  n <- nrow(z); T <- ncol(z)
  obs <- matrix(TRUE, n, T)
  gap_left <- integer(n)
  for (t in 2:T) {
    in_gap <- gap_left > 0L
    obs[in_gap, t] <- FALSE
    gap_left[in_gap] <- gap_left[in_gap] - 1L
    can <- which(!in_gap & obs[, t - 1L])
    if (length(can)) {
      onset <- stats::runif(length(can)) < stats::plogis(base_logit + person_effect[can] + beta * z[can, t - 1L])
      who <- can[onset]
      if (length(who)) {
        L <- pmin(1L + stats::rgeom(length(who), gap_p), max_gap)
        obs[who, t] <- FALSE
        gap_left[who] <- L - 1L
      }
    }
  }
  obs
}

# Drop the trailing missing days of each row so the window ends on an observed day
# (window (a): first to last record). Returns a list of observed day indices.
observed_days_by_person <- function(obs, window_len) {
  lapply(seq_len(nrow(obs)), function(i) {
    idx <- which(obs[i, seq_len(window_len[i])])
    idx
  })
}

simulate_study <- function(n_per_cohort = c(Insomnia = 120, Diabetes = 120, Healthy = 240, Prostate = 60),
                           n_days = c(60L, 140L), beta = -0.5, seed = 20260511L,
                           base_logit = c(activity = -3.9, sleep = -3.0, hr = -3.4)) {
  with_seed(seed, {
    cohort <- rep(names(n_per_cohort), n_per_cohort)
    n  <- length(cohort)
    pid <- 1000L + seq_len(n)
    T  <- max(n_days)
    win <- sample(n_days[1]:n_days[2], n, replace = TRUE)
    start <- as.Date("2022-01-01") + sample(0:365, n, replace = TRUE)

    # latent log steps: person level + AR(1) day-to-day noise; z = within-person z-score
    mu <- stats::rnorm(n, 8.8, 0.4)
    e  <- matrix(0, n, T); e[, 1] <- stats::rnorm(n)
    for (t in 2:T) e[, t] <- 0.5 * e[, t - 1] + sqrt(1 - 0.25) * stats::rnorm(n)
    logsteps <- mu + 0.35 * e
    z <- e                                   # unit-variance within-person deviation

    cohort_shift <- c(Insomnia = 0.15, Diabetes = 0.25, Healthy = -0.1, Prostate = 0)[cohort]
    pe <- stats::rnorm(n, 0, 0.5) + cohort_shift
    obs <- lapply(base_logit, function(b) simulate_presence(z, b, pe, beta))
    names(obs) <- names(base_logit)

    # helper: long table of observed (person, day index) for a stream, trimmed so that
    # the last day of each person's window is observed
    long_obs <- function(stream) {
      m <- obs[[stream]]
      rbindlist(lapply(seq_len(n), function(i) {
        idx <- which(m[i, seq_len(win[i])])
        if (!length(idx)) return(NULL)
        data.table(i = i, day = idx)
      }))
    }

    # ---- activity ------------------------------------------------------------------
    a <- long_obs("activity")
    a[, `:=`(person_id = pid[i], date = start[i] + day - 1L, steps = round(exp(logsteps[cbind(i, day)])))]
    a[, `:=`(very_active_minutes = stats::rpois(.N, pmax(steps / 1500, 0.1)),
             fairly_active_minutes = stats::rpois(.N, pmax(steps / 1000, 0.1)),
             lightly_active_minutes = stats::rpois(.N, 200),
             sedentary_minutes = pmin(1440L, round(stats::rnorm(.N, 900 - steps / 40, 60))),
             calories_out = round(1800 + steps / 6 + stats::rnorm(.N, 0, 80)),
             elevation = NA_real_, floors = NA_real_)]
    activity_raw <- a[, .(person_id, date, steps, calories_out, elevation, fairly_active_minutes,
                          floors, lightly_active_minutes, sedentary_minutes, very_active_minutes)]

    # ---- sleep (main + occasional nap rows; sleep_date naming as in the export) -------
    s <- long_obs("sleep")
    s[, `:=`(person_id = pid[i], sleep_date = start[i] + day - 1L,
             minute_in_bed = round(stats::rnorm(.N, 450, 60)))]
    s[, `:=`(minute_asleep = round(minute_in_bed * stats::runif(.N, 0.8, 0.95)),
             is_main_sleep = "true", minute_restless = NA_real_)]
    s[, `:=`(minute_awake = minute_in_bed - minute_asleep, minute_deep = round(minute_asleep * 0.2),
             minute_light = round(minute_asleep * 0.6), minute_rem = round(minute_asleep * 0.2),
             minute_wake = minute_in_bed - minute_asleep, minute_after_wakeup = 0)]
    naps <- s[stats::runif(.N) < 0.10]
    naps[, `:=`(is_main_sleep = "false", minute_in_bed = 40, minute_asleep = 35,
                minute_awake = 5, minute_deep = 5, minute_light = 25, minute_rem = 5, minute_wake = 5)]
    sleep_raw <- rbind(s, naps)[, .(person_id, sleep_date, is_main_sleep, minute_in_bed, minute_asleep,
                                    minute_after_wakeup, minute_awake, minute_restless, minute_deep,
                                    minute_light, minute_rem, minute_wake)]

    # ---- heart rate: four zone rows per day, Peak often absent -------------------------
    h <- long_obs("hr")
    h[, `:=`(person_id = pid[i], date = start[i] + day - 1L)]
    h[, `:=`(fat = round(stats::runif(.N, 20, 90)), car = round(stats::runif(.N, 0, 30)),
             peak = round(stats::runif(.N, 0, 8)))]
    h[, `:=`(oor = pmax(0, round(stats::rnorm(.N, 1405, 8)) - fat - car - peak))]
    zone_rows <- function(zone, mins, lo, hi, kcal, drop_zero = FALSE) {
      keep <- if (drop_zero) mins > 0 else rep(TRUE, nrow(h))
      h[keep, .(person_id, date, zone_name = zone, min_heart_rate = lo, max_heart_rate = hi,
                minute_in_zone = mins[keep], calorie_count = kcal * mins[keep])]
    }
    hr_raw <- rbind(
      zone_rows("Out of Range", h$oor,  30,  108, 1.5),
      zone_rows("Fat Burn",     h$fat,  108, 130, 5),
      zone_rows("Cardio",       h$car,  130, 160, 8),
      zone_rows("Peak",         h$peak, 160, 220, 10, drop_zero = TRUE))

    # ---- people ---------------------------------------------------------------------------------
    sex <- ifelse(cohort == "Prostate", "Male", sample(c("Female", "Male"), n, TRUE, prob = c(0.6, 0.4)))
    age <- ifelse(cohort == "Prostate", sample(66:85, n, TRUE), sample(22:85, n, TRUE))
    meta <- data.table(person_id = pid,
                       date_of_birth = format(start - round(age * 365.25), "%Y-%m-%d 00:00:00 UTC"),
                       race = sample(c("White", "Black or African American", "Asian",
                                       "More than one population", "None of these"), n, TRUE,
                                     prob = c(0.7, 0.1, 0.05, 0.05, 0.1)),
                       ethnicity = "Not Hispanic or Latino", sex_at_birth = sex,
                       self_reported_category = "White")
    final_pid <- data.table(person_id = pid, cohort = cohort)
    cci <- data.table(person_id = pid,
                      cci_score = stats::rpois(n, c(Insomnia = 3, Diabetes = 5, Healthy = 1, Prostate = 4)[cohort]))
    list(raw = list(activity = activity_raw, sleep = sleep_raw, hr = hr_raw),
         meta = meta, final_pid = final_pid, cci = cci, truth = list(beta = beta))
  })
}

# Convenience: simulate, run the ETL, build persons. Returns what the pipeline needs.
simulate_inputs <- function(cfg, ...) {
  sim <- simulate_study(...)
  day_tables <- lapply(setNames(names(sim$raw), names(sim$raw)),
                       function(s) build_day_table(sim$raw[[s]], s, cfg))
  persons <- build_person_table(sim$meta, sim$final_pid, sim$cci,
                                first_record_dates(day_tables), cfg)
  list(day_tables = day_tables, persons = persons, sim = sim)
}
