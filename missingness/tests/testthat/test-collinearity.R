test_that("gvif_from_vcov equals the textbook VIF for single-coefficient terms and the Fox-Monette GVIF for a factor", {
  set.seed(1); n <- 400
  d <- data.table(x1 = rnorm(n)); d[, x2 := 0.9 * x1 + sqrt(1 - 0.81) * rnorm(n)]; d[, x3 := rnorm(n)]
  d[, f := factor(sample(c("a", "b", "c"), n, TRUE))]
  d[, y := x1 + x2 + x3 + as.integer(f) + rnorm(n)]
  fit <- stats::lm(y ~ x1 + x2 + x3 + f, d)
  groups <- coef_term_groups(setdiff(names(coef(fit)), "(Intercept)"), formula_term_labels("y ~ x1 + x2 + x3 + f"))
  expect_equal(unname(groups[c("fb", "fc")]), c("f", "f"))
  g <- gvif_from_vcov(stats::vcov(fit), groups)
  vif_manual <- function(v) { r2 <- summary(stats::lm(stats::reformulate(setdiff(c("x1", "x2", "x3"), v), v), d))$r.squared; 1 / (1 - r2) }
  for (v in c("x1", "x2", "x3")) expect_equal(g[term == v, gvif], vif_manual(v), tolerance = 0.05)      # design-matrix VIF vs VIF from the estimates' covariance (same up to the factor)
  expect_gt(g[term == "x1", vif_equiv], 3); expect_lt(g[term == "x3", vif_equiv], 1.5)               # the planted correlation is found, an independent predictor is not flagged
  expect_equal(g[term == "f", df], 2L)
  # GVIF of the factor from the determinants of the correlation matrix of the coefficient estimates
  R <- stats::cov2cor(stats::vcov(fit)[-1, -1]); ix <- c("fb", "fc")
  expect_equal(g[term == "f", gvif], det(R[ix, ix]) * det(R[-match(ix, rownames(R)), -match(ix, rownames(R))]) / det(R), tolerance = 1e-8)
  expect_equal(g[term == "f", vif_equiv], g[term == "f", gvif_adj]^2)
  expect_equal(as.character(g[term == "x1", flag]), if (g[term == "x1", vif_equiv] >= 10) "severe" else if (g[term == "x1", vif_equiv] >= 5) "high" else "ok")
})

test_that("flags use the thresholds 5 and 10 and a near-duplicate predictor is flagged severe", {
  set.seed(2); n <- 300; a <- rnorm(n); b <- a + rnorm(n, sd = 0.05); cc <- rnorm(n)
  fit <- stats::lm(rnorm(n) ~ a + b + cc)
  g <- gvif_from_vcov(stats::vcov(fit), stats::setNames(c("a", "b", "cc"), c("a", "b", "cc")))
  expect_true(all(as.character(g[term %in% c("a", "b"), flag]) == "severe")); expect_equal(as.character(g[term == "cc", flag]), "ok")
  expect_equal(c(VIF_HIGH, VIF_SEVERE), c(5, 10))
  expect_equal(nrow(gvif_from_vcov(matrix(1, 1, 1, dimnames = list("a", "a")), c(a = "a"))), 0L)        # one term: nothing to diagnose
})

test_that("term-level Wald tests are the joint test of all coefficients of a term (spline / factor), not of single coefficients", {
  set.seed(3); n <- 500
  d <- data.table(x = rnorm(n), f = factor(sample(letters[1:4], n, TRUE)))
  d[, y := 0.8 * x + rnorm(n)]
  fit <- stats::lm(y ~ x + f, d)
  groups <- coef_term_groups(setdiff(names(coef(fit)), "(Intercept)"), formula_term_labels("y ~ x + f"))
  w <- wald_term_tests(coef(fit), stats::vcov(fit), groups)
  ix <- grep("^f", names(coef(fit)), value = TRUE)
  expect_equal(w[term == "f", chisq], as.numeric(t(coef(fit)[ix]) %*% solve(stats::vcov(fit)[ix, ix]) %*% coef(fit)[ix]))
  expect_equal(w[term == "f", df], 3L); expect_lt(w[term == "x", p], 1e-10); expect_gt(w[term == "f", p], 0.001)
  expect_equal(w[term == "x", chisq], (coef(fit)[["x"]] / sqrt(stats::vcov(fit)["x", "x"]))^2)
  expect_equal(formula_term_labels("y ~ a + splines::ns(t, 3) + (1 | person_id)"), c("a", "splines::ns(t, 3)"))
  expect_equal(unname(coef_term_groups(c("splines::ns(t, 3)2", "factor(year_t)2021", "cohortHealthy"), c("splines::ns(t, 3)", "factor(year_t)", "cohort"))),
               c("splines::ns(t, 3)", "factor(year_t)", "cohort"))
})

test_that("every fitted onset model reports GVIF per term (factors / splines as one term) and the joint Wald tests, from its own fit", {
  s <- small_study(-0.5); r <- fit_onset_spec(s$b, primary_spec(s$cfg), s$cfg)
  expect_true(all(c("w_act_steps", "b_act_steps", "cohort", "dow_t", "splines::ns(study_day, 3)", "cci", "prev_gap_len_l") %in% r$gvif$term))
  expect_equal(r$gvif[term == "cohort", df], 3L); expect_equal(r$gvif[term == "dow_t", df], 6L)
  expect_true(all(c("spec_id", "stream", "role", "standardization", "vif_equiv", "flag") %in% names(r$gvif)))
  expect_equal(r$wald[term == "w_act_steps", chisq], (r$coef[term == "w_act_steps", estimate] / r$coef[term == "w_act_steps", se])^2, tolerance = 1e-6)   # joint test of one term = Wald z^2
  expect_lt(r$gvif[term == "w_act_steps", vif_equiv], 5)                      # the within-person predictors of the simulation are not collinear
  tab <- build_collinearity_table(list(onset_collinearity = r$gvif))
  expect_true(all(c("analysis", "figure", "figure_id", "stratum") %in% names(tab))); expect_true(any(tab$figure)); expect_true(all(tab[figure == TRUE, figure_id] %in% c("onset_primary", "onset_full_period", "model_B", "model_C")))
  expect_equal(unique(tab$analysis), "Primary onset model (gap >= 1 day)")
})
