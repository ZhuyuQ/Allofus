# ==============================================================================
# collinearity.R
# ------------------------------------------------------------------------------
# Multicollinearity diagnostics and term-level tests for every fitted model with several fixed effects.
#
# Method (the same for every model family: lmer, glmer, clmm, glmmTMB)
#   * V = covariance matrix of the estimated FIXED effects (intercept and, for the ordinal model, thresholds left out),
#     R = cov2cor(V). This is what car::vif() and performance::check_collinearity() use for GLMs / mixed models: it
#     measures how much the uncertainty of an estimate is inflated by the correlation with the other estimates, with the
#     weights of the model (not only the correlation of the raw columns).
#   * A predictor that needs several coefficients (a factor, a spline, calendar years) is ONE term: its generalized VIF is
#     GVIF = det(R_term) * det(R_other terms) / det(R)   (Fox & Monette, 1992), with df = number of coefficients of the term.
#   * To compare terms of different df, GVIF^(1/(2 df)) is on the scale of a standard error ratio; its square,
#       vif_equiv = GVIF^(1/df),
#     equals the ordinary VIF for a one-coefficient term and is read against the usual VIF thresholds.
# Thresholds (conventions, not tests; cfg-free constants below): vif_equiv >= 5 "high", >= 10 "severe".
# Interpretation: a high VIF inflates the standard error of that term by sqrt(VIF) (e.g. 5 -> 2.2 times) and makes its
# estimate unstable against small changes of the data; it does NOT bias the other coefficients or the model's predictions,
# and nothing is removed or re-specified because of it. Terms of an interaction would be collinear with their main effects by
# construction; the cohort-interaction models of this pipeline are only used for likelihood-ratio tests, so the diagnostic is reported for
# the main-effects model. The run-in interaction models are reported too: run-in is centred, so its product with the steps term is nearly
# uncorrelated with both.
# ==============================================================================

VIF_HIGH   <- 5
VIF_SEVERE <- 10

# Term labels of a model formula (the random-effect part "(1 | person_id)" removed).
formula_term_labels <- function(formula_text) {
  lab <- attr(stats::terms(stats::as.formula(formula_text)), "term.labels")
  lab[!grepl("|", lab, fixed = TRUE)]
}

# Which model term each coefficient belongs to: the longest term label the coefficient name starts with
# ("cohortDiabetes" -> "cohort", "splines::ns(study_day, 3)2" -> "splines::ns(study_day, 3)").
coef_term_groups <- function(coef_names, term_labels) {
  lab <- term_labels[order(-nchar(term_labels))]
  vapply(coef_names, function(cn) { hit <- lab[startsWith(cn, lab)]; if (length(hit)) hit[1L] else cn }, character(1))
}

# Generalized VIF per term from the covariance matrix of the fixed effects. `groups`: named character vector, coefficient -> term.
gvif_from_vcov <- function(V, groups) {
  V <- as.matrix(V)
  keep <- setdiff(rownames(V), "(Intercept)")
  keep <- keep[keep %in% names(groups)]
  if (length(keep) < 2L) return(data.table(term = character(), df = integer(), gvif = numeric(), gvif_adj = numeric(),
                                          vif_equiv = numeric(), flag = character()))
  V <- V[keep, keep, drop = FALSE]
  ok <- is.finite(diag(V)) & diag(V) > 0
  V <- V[ok, ok, drop = FALSE]; g <- unname(groups[rownames(V)])
  R <- stats::cov2cor(V)
  logdet <- function(M) if (nrow(M) == 0L) 0 else { d <- determinant(M, logarithm = TRUE); if (d$sign < 0) NA_real_ else as.numeric(d$modulus) }
  ld_all <- logdet(R)
  rows <- lapply(unique(g), function(tm) {
    ix <- which(g == tm)
    gv <- exp(logdet(R[ix, ix, drop = FALSE]) + logdet(R[-ix, -ix, drop = FALSE]) - ld_all)
    data.table(term = tm, df = length(ix), gvif = gv)
  })
  out <- rbindlist(rows)
  out[, `:=`(gvif_adj = gvif^(1 / (2 * df)), vif_equiv = gvif^(1 / df))]
  out[, flag := fifelse(is.na(vif_equiv), "not_computed", fifelse(vif_equiv >= VIF_SEVERE, "severe",
                                                                  fifelse(vif_equiv >= VIF_HIGH, "high", "ok")))]
  out[]
}

# Joint Wald test of every model term (all its coefficients at once): chisq = b' V^-1 b on df = number of coefficients.
# This is the only valid reading of a spline or of a factor with many levels; individual coefficients are not.
wald_term_tests <- function(beta, V, groups) {
  V <- as.matrix(V)
  keep <- intersect(setdiff(names(beta), "(Intercept)"), rownames(V))
  keep <- keep[keep %in% names(groups)]
  if (!length(keep)) return(data.table(term = character(), df = integer(), chisq = numeric(), p = numeric()))
  rows <- lapply(unique(groups[keep]), function(tm) {
    ix <- keep[groups[keep] == tm]
    b <- beta[ix]; Vi <- V[ix, ix, drop = FALSE]
    chisq <- tryCatch(as.numeric(t(b) %*% solve(Vi, b)), error = function(e) NA_real_)
    data.table(term = tm, df = length(ix), chisq = chisq, p = if (is.finite(chisq) && chisq >= 0) stats::pchisq(chisq, length(ix), lower.tail = FALSE) else NA_real_)
  })
  rbindlist(rows)
}

# Both diagnostics of one fitted model: list(gvif = , wald = ). `beta` named fixed effects, `V` their covariance matrix,
# `formula_text` the fixed-effect part of the model formula (term labels).
fixed_effect_diagnostics <- function(beta, V, formula_text) {
  groups <- coef_term_groups(setdiff(names(beta), "(Intercept)"), formula_term_labels(formula_text))
  list(gvif = gvif_from_vcov(V, groups), wald = wald_term_tests(beta, V, groups))
}

# Attach identifying columns to a diagnostics table (NULL / empty stay empty).
tag_diagnostics <- function(tab, ...) {
  if (is.null(tab) || nrow(tab) == 0L) return(NULL)
  ids <- list(...)
  cbind(as.data.table(ids), tab)
}

# ---- the tables written by the pipeline -----------------------------------------------------------------------------------

nz <- function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L

# Readable name of a run-in model (panel title of the collinearity figure).
runin_analysis_label <- function(variant, what) paste0("Run-in analysis (", RUNIN_LABELS[as.character(variant)], "): ", what)

# Readable name of an onset specification (also the panel title of the collinearity figure).
onset_analysis_label <- function(role, set_name, k, adjusted) {
  fifelse(role == "primary", "Primary onset model (gap >= 1 day)",
          fifelse(role == "secondary_threshold", sprintf("Model D: gap >= %d days", k),
                  sprintf("Exploratory onset: %s, %s", spec_set_label(set_name), fifelse(adjusted, "adjusted", "unadjusted"))))
}

# All GVIF results of all fitted models of the MAIN analysis in one table (long format), with a readable `analysis`, and `figure = TRUE` + `figure_id` for the
# models that get a figure of their own, collinearity_gvif_<figure_id>_<stream> (primary onset, full-period standardization, models B and C).
build_collinearity_table <- function(res) {
  keep <- c("stream", "model_family", "analysis", "stratum", "term", "df", "gvif", "gvif_adj", "vif_equiv", "flag")
  parts <- list()
  if (nz(res$onset_collinearity)) {
    o <- res$onset_collinearity[baseline == "whole"]
    parts[[1L]] <- o[, .(stream, model_family = "gap onset", analysis = onset_analysis_label(role, set_name, k, adjusted), stratum = "all", term, df, gvif, gvif_adj, vif_equiv, flag)]
  }
  if (nz(res$sensitivity_collinearity)) {
    sc <- res$sensitivity_collinearity
    parts[[2L]] <- sc[, .(stream, model_family = "gap onset (sensitivity)", analysis = paste0("Sensitivity: ", nice_label(analysis)), stratum = "all", term, df, gvif, gvif_adj, vif_equiv, flag)]
  }
  if (nz(res$gap_collinearity)) {
    g <- res$gap_collinearity
    parts[[3L]] <- g[, .(stream, model_family = paste("model", model), analysis = paste0("Model ", model, ": ", model_id), stratum = "all", term, df, gvif, gvif_adj, vif_equiv, flag)]
  }
  out <- rbindlist(parts, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  ids <- c("Primary onset model (gap >= 1 day)" = "onset_primary", "Sensitivity: Full Period Standardization" = "onset_full_period",
           "Model B: B_back2_adj_pre" = "model_B", "Model C: C_adj" = "model_C")
  out[, figure_id := unname(ids[analysis])]
  out[, figure := !is.na(figure_id)]
  setcolorder(out, c(keep, "figure", "figure_id"))
  out[]
}

# The same for the run-in analysis, kept apart from the main diagnostics (its own tables and figures). Every variant is listed; the models with the run-in
# and its interaction (the primary run-in analysis) are drawn, one figure per model.
build_runin_collinearity <- function(res) {
  keep <- c("stream", "model_family", "analysis", "variant", "term", "df", "gvif", "gvif_adj", "vif_equiv", "flag")
  parts <- list()
  if (nz(res$runin_onset_collinearity)) {
    w <- res$runin_onset_collinearity
    parts[[1L]] <- w[, .(stream, model_family = "gap onset", analysis = runin_analysis_label(runin_variant, "primary onset model"), variant = runin_variant, term, df, gvif, gvif_adj, vif_equiv, flag,
                         figure_id = fifelse(runin_variant == "interaction", "onset", NA_character_))]
  }
  if (nz(res$runin_gap_collinearity)) {
    w <- res$runin_gap_collinearity
    parts[[2L]] <- w[, .(stream, model_family = paste("model", model), analysis = runin_analysis_label(runin_variant, paste("model", model, model_id)), variant = runin_variant, term, df, gvif, gvif_adj, vif_equiv, flag,
                         figure_id = fifelse(runin_variant == "interaction", paste0("model_", model), NA_character_))]
  }
  out <- rbindlist(parts, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  out[, figure := !is.na(figure_id)]
  setcolorder(out, c(keep, "figure", "figure_id"))
  out[]
}

collinearity_for_figure <- function(tab) if (nz(tab)) tab[figure == TRUE] else tab

# The joint Wald tests of every model term (all coefficients of a factor / spline / year set together) in one table.
build_term_test_table <- function(res) {
  parts <- list()
  if (nz(res$onset_term_tests)) {
    o <- res$onset_term_tests[baseline == "whole"]
    parts[[1L]] <- o[, .(stream, model_family = "gap onset", analysis = onset_analysis_label(role, set_name, k, adjusted), stratum = "all", term, df, chisq, p)]
  }
  if (nz(res$gap_term_tests)) {
    g <- res$gap_term_tests
    parts[[2L]] <- g[, .(stream, model_family = paste("model", model), analysis = paste0("Model ", model, ": ", model_id), stratum = "all", term, df, chisq, p)]
  }
  rbindlist(parts, fill = TRUE)
}

# The joint Wald tests of the run-in models, apart from the main tests.
build_runin_term_tests <- function(res) {
  parts <- list()
  if (nz(res$runin_onset_term_tests)) {
    w <- res$runin_onset_term_tests
    parts[[1L]] <- w[, .(stream, model_family = "gap onset", analysis = runin_analysis_label(runin_variant, "primary onset model"), variant = runin_variant, term, df, chisq, p)]
  }
  if (nz(res$runin_gap_term_tests)) {
    w <- res$runin_gap_term_tests
    parts[[2L]] <- w[, .(stream, model_family = paste("model", model), analysis = runin_analysis_label(runin_variant, paste("model", model, model_id)), variant = runin_variant, term, df, chisq, p)]
  }
  rbindlist(parts, fill = TRUE)
}
