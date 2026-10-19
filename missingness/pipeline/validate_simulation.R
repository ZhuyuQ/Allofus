# validate_simulation.R - validation of the primary estimate on SIMULATED data with a KNOWN effect:
# log-odds -0.5 per +1 within-person SD of log steps. 24 replicates of a 300-person study, pre-gap vs full-period standardization:
# bias, empirical SD of the estimate, mean model SE, coverage of the 95% interval, and the VIF-equivalents of the within / between steps terms.
#   S_OUT=<dir> Rscript pipeline/validate_simulation.R      (about 25 minutes on 3 cores; writes <dir>/validate_sim.tsv)
root <- if (file.exists("load_missingness.R")) "." else "missingness"
source(file.path(root, "load_missingness.R")); load_missingness(root)
cfg <- make_config(list(boot = list(B = 30L), resources = list(isolate_fits = FALSE)))
one <- function(rep) {
  inp <- suppressMessages(simulate_inputs(cfg, n_per_cohort = c(Insomnia = 60, Diabetes = 60, Healthy = 120, Prostate = 30), beta = -0.5, seed = 1000 + rep))
  out <- lapply(c("pre_gap", "full_period"), function(m) {
    b <- suppressMessages(prepare_stream_bundle("activity", inp$day_tables, inp$persons, cfg, standardization = m))
    r <- fit_onset_spec(b, make_onset_grid(cfg)[role == "primary" & stream == "activity"], cfg)
    e <- r$coef[term == "w_act_steps"]
    data.table(rep = rep, method = m, est = e$estimate, se = e$se, n_events = e$n_events, sd_within = b$scales[predictor == "act_steps", sd_within],
               vif_between_steps = r$gvif[term == "b_act_steps", vif_equiv], vif_within_steps = r$gvif[term == "w_act_steps", vif_equiv])
  })
  rbindlist(out)
}
res <- rbindlist(parallel::mclapply(1:24, function(i) tryCatch(one(i), error = function(e) NULL), mc.cores = 3L))
fwrite(res, file.path(Sys.getenv("S_OUT", tempdir()), "validate_sim.tsv"), sep = "\t")
res[, covers := (est - 1.96 * se <= -0.5) & (est + 1.96 * se >= -0.5)]
print(res[, .(reps = .N, mean_est = round(mean(est), 3), bias = round(mean(est) + 0.5, 3), empirical_sd = round(sd(est), 3), mean_se = round(mean(se), 3),
              coverage_95 = round(mean(covers), 2), mean_events = round(mean(n_events)), sd_within = round(mean(sd_within), 3),
              vif_within = round(mean(vif_within_steps), 2), vif_between = round(mean(vif_between_steps), 1)), by = method])
