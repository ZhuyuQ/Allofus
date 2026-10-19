# ==============================================================================
# run_demo.R  - end-to-end smoke test on SIMULATED data (no real data needed).
#   Rscript pipeline/run_demo.R <output_dir>
# The simulated truth is: +1 SD higher steps than usual -> log-odds of a gap tomorrow
# -0.5 (odds ratio 0.61). The primary estimate should be close to that.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1L] else file.path(tempdir(), "missingness_demo")
root <- if (file.exists("load_missingness.R")) "." else "missingness"
source(file.path(root, "load_missingness.R")); load_missingness(root)

cfg <- make_config(list(boot = list(B = 50L)))
inp <- simulate_inputs(cfg)
res <- run_missingness_pipeline(inp$day_tables, inp$persons, cfg, outdir)
print(res$primary_test)
cat("results written to:", outdir, "\n")
