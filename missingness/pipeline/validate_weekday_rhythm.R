# ==============================================================================
# validate_weekday_rhythm.R - does the weekly-rhythm analysis recover a KNOWN truth? (simulated data only; nothing real is read)
# Run from missingness/:   Rscript pipeline/validate_weekday_rhythm.R [replicates] [people] [days] [outdir]
#
# Three scenarios (simulate_gap_process):
#   calendar  a shared weekday effect on gap starts (Friday low, Monday / Sunday high)       -> H1 detected; H2: the weekday pattern stays after days since return
#   cycle     the gap-start chance peaks 6-7 days since return, no weekday effect              -> H1 about null; H2 about null; the days-since-return curve recovers the peak
#   personal  a weekly rhythm with a random peak day per person, no shared effect             -> H1 about null; H3 detected; the person peak-to-trough percentiles cover the planted ones
# Per replicate the pipeline's own functions are used: wr_gap_start_rows, wr_fit_value (S0, S1, S2, S3, S5), wr_tests. "Truth" of the weekday ratio is
# exp(effect_k) / mean(exp(effect)) in the calendar scenario (exact in the low-rate regime: the gap-start chance is a few percent) and 1 for every weekday
# in the other two (no shared weekday effect).
# ==============================================================================
source("load_missingness.R"); load_missingness(".")
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
R <- as.integer(if (length(args) >= 1L) args[1L] else 6L); N <- as.integer(if (length(args) >= 2L) args[2L] else 200L)
T_days <- as.integer(if (length(args) >= 3L) args[3L] else 150L); outdir <- if (length(args) >= 4L) args[4L] else tempdir()

cfg <- make_config(list(boot = list(B = 30L)))
cfg$model$adjust_device <- FALSE
cfg$weekday_rhythm[c("n_profile_rows", "n_beta_draws", "n_re_draws", "n_curve_draws", "n_person_draws")] <- list(800L, 200L, 40L, 40L, 5000L)

WD <- c(0.4, 0.1, 0, 0.1, -0.9, -0.4, 0.3)
SCEN <- list(calendar = list(weekday_effect = WD),
             cycle = list(dsr_effect = function(d) ifelse(d %in% 6:7, 1.2, 0)),
             personal = list(personal_amp = 0.5))
TRUTH_RATIO <- list(calendar = exp(WD) / mean(exp(WD)), cycle = rep(1, 7), personal = rep(1, 7))

fit_set <- function(rows, specs, ids, seed) {                       # S1 decides the random-effects structure, the others follow it
  out <- list(); re <- NULL
  for (id in ids) {
    sp <- specs[fit_id == id]; dat <- wr_apply_filter(rows, sp, "activity", cfg)
    out[[id]] <- wr_fit_value(dat, sp, "activity", cfg, re, seed + match(id, specs$fit_id), decide = (id == "S1"))
    if (id == "S1") re <- out$S1$info$re_structure
  }
  out
}

validate_one <- function(scenario, rep) {
  sim <- do.call(simulate_gap_process, c(list(n = N, T = T_days, seed = 1000L * match(scenario, names(SCEN)) + rep), SCEN[[scenario]]))
  rows <- wr_gap_start_rows(sim$day_tables$activity, sim$persons, "activity", cfg)
  v <- fit_set(rows, wr_specs(), c("S1", "S0", "S2", "S3", "S5"), rep)
  fits <- rbindlist(lapply(v, `[[`, "info"))
  t <- wr_tests(fits[, stream := "activity"])
  pr <- v$S1$profile; cv <- v$S2$curve; tr <- TRUTH_RATIO[[scenario]]
  data.table(scenario = scenario, rep = rep, p_H1 = t[id == "H1s", p], p_H2 = t[id == "H2", p], p_H3 = t[id == "H3", p], p_H2r = t[id == "H2r", p],
             re_structure = fits[fit_id == "S1", re_structure], ptr = v$S1$summary$peak_to_trough,
             ratio_bias = mean(pr$ratio_to_average_day - tr), ratio_cover = mean(pr$ratio_low <= tr & pr$ratio_high >= tr),
             fri_ratio = pr[weekday == 4L, ratio_to_average_day], curve_peak_level = as.integer(cv[which.max(pct), sub("\\+", "", level)]),
             pers_median = v$S1$personal$person_ptr_median, pers_p10 = v$S1$personal$person_ptr_p10, pers_p90 = v$S1$personal$person_ptr_p90,
             shrinkage = 1 - log(v$S3$summary$peak_to_trough) / log(v$S1$summary$peak_to_trough))
}

grid <- CJ(scenario = names(SCEN), rep = seq_len(R), sorted = FALSE)
res <- rbindlist(parallel::mclapply(seq_len(nrow(grid)), function(i) tryCatch(validate_one(grid$scenario[i], grid$rep[i]), error = function(e) data.table(scenario = grid$scenario[i], rep = grid$rep[i], err = conditionMessage(e))),
                                    mc.cores = max(1L, min(4L, parallel::detectCores() - 1L))), fill = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!"err" %in% names(res)) res[, err := NA_character_]                 # no replicate failed
fwrite(res, file.path(outdir, "validation_weekday_rhythm.tsv"), sep = "\t")

truth_pers <- local({ set.seed(1); x <- exp(2 * 0.5 * stats::runif(1e5, 0.5, 1.5)); stats::quantile(x, c(0.1, 0.5, 0.9), names = FALSE) })
summ <- res[is.na(err), .(replicates = .N, detect_H1 = mean(p_H1 < 0.05), detect_H2 = mean(p_H2 < 0.05), detect_H3 = mean(p_H3 < 0.05), detect_H2r = mean(p_H2r < 0.05),
                          mean_ptr = mean(ptr), mean_ratio_bias = mean(ratio_bias), mean_ratio_cover = mean(ratio_cover), mean_fri_ratio = mean(fri_ratio),
                          mean_shrinkage = mean(shrinkage), curve_peak_at_6_or_7 = mean(curve_peak_level %in% 6:7),
                          pers_p10 = mean(pers_p10), pers_median = mean(pers_median), pers_p90 = mean(pers_p90)), by = scenario]
cat(sprintf("\nreplicates per scenario: %d requested; people %d, days %d. Truth: Friday ratio %.3f (calendar); personal peak-to-trough p10 / median / p90 = %.2f / %.2f / %.2f\n",
            R, N, T_days, TRUTH_RATIO$calendar[5L], truth_pers[1L], truth_pers[2L], truth_pers[3L]))
print(summ[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(summ, file.path(outdir, "validation_weekday_rhythm_summary.tsv"), sep = "\t")
