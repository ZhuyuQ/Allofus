# ==============================================================================
# load_missingness.R
# Source every file of the missingness pipeline in dependency order.
#   source("<path>/missingness/load_missingness.R"); load_missingness("<path>/missingness")
# ==============================================================================
load_missingness <- function(root = ".") {
  files <- c("utils", "config", "io", "etl", "persons", "covariates", "gaps", "predictors",
             "device", "outcome_variables", "collinearity", "onset_model", "sensitivity", "gap_models", "runin", "summarise", "temporal",
             "descriptive", "plots", "plots_effects", "report", "bq_export", "simulate", "pipeline")
  for (f in files) source(file.path(root, "R", paste0(f, ".R")), chdir = FALSE)
  invisible(files)
}
