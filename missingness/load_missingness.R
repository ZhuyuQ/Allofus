# ==============================================================================
# load_missingness.R
# Source every file of the missingness pipeline in dependency order.
#   source("<path>/missingness/load_missingness.R"); load_missingness("<path>/missingness")
# ==============================================================================
load_missingness <- function(root = ".") {
  files <- c("utils", "config", "io", "etl", "persons", "gaps", "predictors",
             "onset_model", "sensitivity", "gap_models", "summarise", "temporal",
             "descriptive", "plots", "report", "bq_export", "simulate", "pipeline")
  for (f in files) source(file.path(root, "R", paste0(f, ".R")), chdir = FALSE)
  invisible(files)
}
