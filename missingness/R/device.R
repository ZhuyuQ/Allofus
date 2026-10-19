# ==============================================================================
# device.R
# ------------------------------------------------------------------------------
# The Fitbit device of each person (optional; exported by Cell 6c, built into L1 by run_etl_device()).
#
# WATCH data only: a person's watch version comes from rows that are neither a scale (device_type SCALE) nor a device that is
# not worn on the wrist (clip-on trackers, the phone-based "MobileTrack"); cfg$device$forms / wrist_forms define this.
# Source: the CDR `device` table, one row per person x device x record date, with device_version ("Charge 5", "Versa 4",
# "MobileTrack" ...) and device_type (TRACKER or SCALE). Scales are not watches, so
#   * a person's WATCH version is taken from wrist-worn rows only (no scale, no clip-on, no phone tracking);
#   * a person is "missing" the watch version when no such row with a version exists: no device record at all, scale
#     only, phone / clip-on only, or rows without a version.
# Person level (device_person): main version = the watch version with the most records (ties: the most recent record, then
# alphabetical), its family (cfg$device$families), the number of different wearable versions, first / last record date.
# Timeline (device_timeline): wearable rows with a record date, to give each DAY the device in use (device_at_dates()).
# The table has far fewer rows than days and some rows have no date, so a date match is not possible; the day takes the
# device of the most recent record ON OR BEFORE it ("as of" join).
#
# Statistics (device_summary), all person level, no model of the gap outcome:
#   * missing share per cohort with Wilson 95% intervals; Pearson chi-square across cohorts, bias-corrected Cramer's V with a
#     bootstrap interval, each cohort against the reference cohort (Fisher exact odds ratio, Newcombe interval of the risk
#     difference, Holm adjustment); a logistic model adjusted for age and sex;
#   * composition of the main device family per cohort (Wilson intervals), chi-square of independence with a Monte-Carlo p
#     value (many sparse cells), Cramer's V with a bootstrap interval, adjusted standardized residuals (Holm) per cell; the same
#     test among people WITH a known device.
# ==============================================================================

DEVICE_STATUS <- c(known = "Watch version known", scale_only = "Scale only", non_wrist = "Phone or clip-on only (not worn on the wrist)",
                   no_version = "Device row without version", no_record = "No device record")

# ---- small statistics helpers ----------------------------------------------------------------------------

# Wilson score interval of a proportion.
wilson_ci <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- ifelse(n > 0, x / n, NA_real_)
  den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.table(est = p, low = pmin(p, pmax(0, centre - half)), high = pmax(p, pmin(1, centre + half)))     # rounding must not put the estimate outside
}

# Newcombe's hybrid score interval of the difference of two proportions p1 - p2.
newcombe_diff_ci <- function(x1, n1, x2, n2, level = 0.95) {
  a <- wilson_ci(x1, n1, level); b <- wilson_ci(x2, n2, level)
  d <- a$est - b$est
  data.table(diff = d, low = d - sqrt((a$est - a$low)^2 + (b$high - b$est)^2),
             high = d + sqrt((a$high - a$est)^2 + (b$est - b$low)^2))
}

# Bias-corrected Cramer's V (Bergsma 2013) of a contingency table; rows / columns without counts are ignored.
cramers_v <- function(tab) {
  tab <- as.matrix(tab); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  n <- sum(tab); r <- nrow(tab); k <- ncol(tab)
  if (n < 2 || r < 2 || k < 2) return(NA_real_)
  e <- outer(rowSums(tab), colSums(tab)) / n
  phi2 <- sum((tab - e)^2 / e) / n
  phi2c <- max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
  kc <- k - (k - 1)^2 / (n - 1); rc <- r - (r - 1)^2 / (n - 1)
  sqrt(phi2c / min(kc - 1, rc - 1))
}

# Percentile bootstrap interval of Cramer's V: every row (cohort) is resampled as a multinomial of its own size (the cohort
# sizes are a feature of the design and stay fixed).
cramers_v_ci <- function(tab, B = 1000L, level = 0.95) {
  tab <- as.matrix(tab); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(c(NA_real_, NA_real_))
  v <- vapply(seq_len(B), function(i) {
    b <- t(vapply(seq_len(nrow(tab)), function(j) as.numeric(stats::rmultinom(1L, sum(tab[j, ]), tab[j, ] / sum(tab[j, ]))), numeric(ncol(tab))))
    cramers_v(b)
  }, numeric(1))
  unname(stats::quantile(v, c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE))
}

# Pearson chi-square of independence with the asymptotic p, a Monte-Carlo p, the smallest expected count and V with its interval.
chisq_table <- function(tab, B_mc = 10000L, B_boot = 1000L) {
  tab <- as.matrix(tab); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2L || ncol(tab) < 2L) {
    return(data.table(n = sum(tab), chisq = NA_real_, df = NA_integer_, p_asymptotic = NA_real_, p_monte_carlo = NA_real_,
                      min_expected = NA_real_, pct_cells_expected_lt5 = NA_real_, cramers_v = NA_real_, v_low = NA_real_, v_high = NA_real_))
  }
  a <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  mc <- suppressWarnings(stats::chisq.test(tab, simulate.p.value = TRUE, B = B_mc))
  ci <- cramers_v_ci(tab, B_boot)
  data.table(n = sum(tab), chisq = unname(a$statistic), df = unname(a$parameter), p_asymptotic = a$p.value, p_monte_carlo = mc$p.value,
             min_expected = min(a$expected), pct_cells_expected_lt5 = 100 * mean(a$expected < 5),
             cramers_v = cramers_v(tab), v_low = ci[1L], v_high = ci[2L])
}

# ---- the device of each person --------------------------------------------------------------------------

# Family of a device version: first matching pattern of cfg$device$families, else other_label; NA stays NA.
device_family <- function(version, cfg) {
  fam <- rep(NA_character_, length(version))
  ok <- !is.na(version)
  fam[ok] <- cfg$device$other_label
  for (nm in rev(names(cfg$device$families))) fam[ok & grepl(cfg$device$families[[nm]], version, ignore.case = TRUE)] <- nm
  fam
}

# How a device version is worn (cfg$device$forms, first match wins; anything else "unclassified"); NA stays NA.
device_form <- function(version, cfg) {
  frm <- rep(NA_character_, length(version))
  ok <- !is.na(version)
  frm[ok] <- "unclassified"
  for (nm in rev(names(cfg$device$forms))) frm[ok & grepl(cfg$device$forms[[nm]], version, ignore.case = TRUE)] <- nm
  frm
}

# Raw device rows -> clean rows (version / type trimmed, empty -> NA, date as Date, scale flag).
normalise_device_rows <- function(raw, cfg) {
  d <- data.table::as.data.table(raw)
  assert_columns(d, c("person_id", "device_version"), "raw device export")
  for (cl in c("device_date", "device_type")) if (!cl %in% names(d)) d[, (cl) := NA]
  out <- d[, .(person_id,
               device_date = suppressWarnings(as.Date(substr(as.character(device_date), 1L, 10L))),
               device_version = trimws(as.character(device_version)),
               device_type = toupper(trimws(as.character(device_type))))]
  out[!nzchar(device_version), device_version := NA_character_]
  out[, is_scale := device_type %in% toupper(cfg$device$scale_types)]
  unique(out)
}

# raw: the device rows of all cohorts; persons: the person table (every person in it gets a row).
# Returns list(person =, timeline =).
build_device_tables <- function(raw, persons, cfg) {
  persons <- data.table::as.data.table(persons)[, .(person_id, cohort)]
  d <- normalise_device_rows(raw, cfg)[person_id %in% persons$person_id]
  rows <- d[, .(n_rows = .N, n_scale_rows = sum(is_scale)), by = person_id]
  nw <- d[!is_scale & !is.na(device_version)]                                       # trackers with a version
  nw[, `:=`(device_form = device_form(device_version, cfg), device_family = device_family(device_version, cfg))]
  nw[, wrist := device_form %in% cfg$device$wrist_forms]
  nonwrist <- nw[wrist == FALSE, .(n_nonwrist_rows = .N), by = person_id]
  w <- nw[wrist == TRUE]                                                              # WATCH rows only: everything below uses these

  by_ver <- w[, .(n_rows = .N, last_date = if (all(is.na(device_date))) as.Date(NA) else max(device_date, na.rm = TRUE)),
              by = .(person_id, device_version, device_family, device_form)]
  data.table::setorderv(by_ver, c("person_id", "n_rows", "last_date", "device_version"), order = c(1L, -1L, -1L, 1L), na.last = TRUE)
  main <- by_ver[, .SD[1L], by = person_id][, .(person_id, device_version_main = device_version, device_family_main = device_family, device_form_main = device_form)]
  span <- w[, .(n_wearable_versions = uniqueN(device_version),
                first_device_date = if (all(is.na(device_date))) as.Date(NA) else min(device_date, na.rm = TRUE),
                last_device_date  = if (all(is.na(device_date))) as.Date(NA) else max(device_date, na.rm = TRUE)), by = person_id]

  p <- merge(persons, rows, by = "person_id", all.x = TRUE)
  p <- merge(p, main, by = "person_id", all.x = TRUE)
  p <- merge(p, span, by = "person_id", all.x = TRUE)
  p <- merge(p, nonwrist, by = "person_id", all.x = TRUE)
  p[is.na(n_nonwrist_rows), n_nonwrist_rows := 0L]
  p[is.na(n_rows), `:=`(n_rows = 0L, n_scale_rows = 0L)]
  p[is.na(n_wearable_versions), n_wearable_versions := 0L]
  p[, status := data.table::fcase(!is.na(device_version_main), DEVICE_STATUS[["known"]],
                                  n_rows == 0L, DEVICE_STATUS[["no_record"]],
                                  n_scale_rows == n_rows, DEVICE_STATUS[["scale_only"]],
                                  n_nonwrist_rows > 0L, DEVICE_STATUS[["non_wrist"]],
                                  default = DEVICE_STATUS[["no_version"]])]
  p[, device_missing := status != DEVICE_STATUS[["known"]]]
  p[, device_family := fifelse(device_missing, cfg$device$missing_label, device_family_main)]
  p[, device_family_main := NULL]
  data.table::setkey(p, person_id)

  tl <- unique(w[!is.na(device_date), .(person_id, device_date, device_version, device_family)])
  data.table::setorder(tl, person_id, device_date, device_version)
  list(person = p[], timeline = tl[])
}

# The device in use on each day: the wearable record of the same person on or before the day ("as of" join, the latest record
# wins; two versions on the same record date -> the alphabetically last).
# days: person_id, date.  backfill = TRUE gives days before the first record the first known device (default: NA).
device_at_dates <- function(days, timeline, backfill = FALSE) {
  d <- data.table::as.data.table(days)[, .(person_id, date)]
  tl <- data.table::as.data.table(timeline)[, .(person_id, device_date, device_version, device_family)]
  data.table::setorder(tl, person_id, device_date, device_version)
  tl <- tl[, .SD[.N], by = .(person_id, device_date)]
  d[, ord__ := .I]
  out <- tl[d, on = .(person_id, device_date = date), roll = TRUE]
  if (backfill) {
    first_known <- tl[, .SD[1L], by = person_id]
    miss <- is.na(out$device_version)
    if (any(miss)) {
      fk <- first_known[match(out$person_id[miss], person_id)]
      data.table::set(out, which(miss), "device_version", fk$device_version)
      data.table::set(out, which(miss), "device_family", fk$device_family)
    }
  }
  data.table::setnames(out, "device_date", "date")
  out[, date := d$date[ord__]]
  data.table::setorder(out, ord__)
  res <- out[, .(person_id, date, device_version, device_family)]
  data.table::setcolorder(res, c("person_id", "date", "device_version", "device_family"))
  res[]
}

# ---- ETL: raw device exports -> L1 ------------------------------------------------------------------------

# Read the device exports of every cohort (Cell 6c wrote them), build the person table and the timeline and save them next to
# the day tables: device_person.rds, device_timeline.rds, device_status.tsv. The day tables and persons.rds are NOT touched, so no
# model result is invalidated. `persons` defaults to <out_dir>/persons.rds (written by run_etl()).
run_etl_device <- function(cfg, out_dir, lister = gsutil_ls, reader = gsutil_fread, persons = NULL) {
  if (is.null(persons)) {
    pf <- file.path(out_dir, "persons.rds")
    if (!file.exists(pf)) stop("persons.rds not found in ", out_dir, ": run the ETL (Cell 8) first", call. = FALSE)
    persons <- readRDS(pf)
  }
  raw <- rbindlist(lapply(cfg$cohorts, function(co)
    tryCatch(read_export_files(export_path(cfg, co, "device"), c(device_version = "character", device_type = "character"),
                               lister = lister, reader = reader),
             error = function(e) stop("device export of ", co, " not found (", conditionMessage(e),
                                      "). Run Cell 6c (run_bq_exports(cfg, final_pid, streams = \"device\")) first.", call. = FALSE))),
    use.names = TRUE, fill = TRUE)
  tabs <- build_device_tables(raw, persons, cfg)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(tabs$person, file.path(out_dir, "device_person.rds"))
  saveRDS(tabs$timeline, file.path(out_dir, "device_timeline.rds"))
  status <- tabs$person[, .N, by = .(cohort, status)][order(cohort, status)]
  fwrite(status, file.path(out_dir, "device_status.tsv"), sep = "\t")
  log_info("device: %d people, %d with a wearable device version, %d timeline records", nrow(tabs$person),
           sum(!tabs$person$device_missing), nrow(tabs$timeline))
  invisible(c(tabs, list(summary = status)))
}

# ---- statistics --------------------------------------------------------------------------------------------

# dp: device_person; persons: the person table (age_first, sex). Returns the device_* tables (see the header).
device_summary <- function(dp, persons, cfg) {
  set.seed(cfg$seed)
  ref <- cfg$reference_cohort %||% cfg$cohorts[1L]
  miss_lab <- cfg$device$missing_label
  p <- merge(data.table::as.data.table(dp)[, .(person_id, cohort, status, device_missing, device_family, device_version_main, n_wearable_versions)],
             data.table::as.data.table(persons)[, .(person_id, age_first, sex)], by = "person_id", all.x = TRUE)
  p[, cohort := factor(as.character(cohort), levels = cfg$cohorts)]
  p <- p[!is.na(cohort)]
  lv <- levels(droplevels(p$cohort))

  # ---- (a) missing device version per cohort ----
  one <- function(x) x[, .(n_people = .N, n_missing = sum(device_missing),
                           n_no_record = sum(status == DEVICE_STATUS[["no_record"]]),
                           n_scale_only = sum(status == DEVICE_STATUS[["scale_only"]]),
                           n_non_wrist_only = sum(status == DEVICE_STATUS[["non_wrist"]]),
                           n_row_without_version = sum(status == DEVICE_STATUS[["no_version"]]))]
  mb <- rbind(p[, one(.SD), by = cohort][order(match(cohort, cfg$cohorts))], cbind(cohort = factor("All cohorts"), one(p)), fill = TRUE)
  ci <- wilson_ci(mb$n_missing, mb$n_people)
  mb[, `:=`(pct_missing = 100 * ci$est, pct_low = 100 * ci$low, pct_high = 100 * ci$high)]
  mb[, cohort := as.character(cohort)]

  # ---- (b) tests of the missing share ----
  t2 <- table(p$cohort, p$device_missing)[lv, , drop = FALSE]
  cs <- chisq_table(t2, cfg$device$mc_B, cfg$device$boot_B)
  tests <- cbind(test = "Missing device version differs between cohorts (Pearson chi-square, 4 x 2)", cs)
  pair <- rbindlist(lapply(setdiff(lv, ref), function(co) {
    a <- p[cohort == co]; b <- p[cohort == ref]
    x <- matrix(c(sum(a$device_missing), sum(!a$device_missing), sum(b$device_missing), sum(!b$device_missing)), 2L, byrow = TRUE)
    ft <- stats::fisher.test(x)
    dd <- newcombe_diff_ci(sum(a$device_missing), nrow(a), sum(b$device_missing), nrow(b))
    data.table(cohort = co, reference = ref, n_people = nrow(a), n_missing = sum(a$device_missing), reference_n_people = nrow(b),
               reference_n_missing = sum(b$device_missing), odds_ratio = unname(ft$estimate), or_low = ft$conf.int[1L], or_high = ft$conf.int[2L],
               risk_diff_pct_points = 100 * dd$diff, rd_low = 100 * dd$low, rd_high = 100 * dd$high, p_fisher = ft$p.value)
  }))
  if (nrow(pair)) pair[, p_holm := stats::p.adjust(p_fisher, "holm")]

  # ---- (c) adjusted model (age, sex); a cohort without events (or with only events) cannot be estimated ----
  dat <- p[!is.na(age_first) & !is.na(sex)]
  ev <- dat[, .(e = sum(device_missing), n = .N), by = cohort]
  usable <- ev[e > 0L & e < n, as.character(cohort)]
  model <- NULL
  if (ref %in% usable && length(usable) >= 2L) {
    d2 <- dat[as.character(cohort) %in% usable]
    d2[, `:=`(cohort = stats::relevel(factor(as.character(cohort), levels = intersect(cfg$cohorts, usable)), ref = ref),
              age10 = age_first / 10, sex = droplevels(factor(sex)))]
    f1 <- if (uniqueN(d2$sex) > 1L) device_missing ~ cohort + age10 + sex else device_missing ~ cohort + age10
    fit <- stats::glm(f1, data = d2, family = stats::binomial())
    fit0 <- stats::update(fit, . ~ . - cohort)
    sm <- stats::coef(summary(fit)); ci <- suppressMessages(stats::confint.default(fit))
    lrt <- stats::anova(fit0, fit, test = "Chisq")
    model <- data.table(term = rownames(sm)[-1L], odds_ratio = exp(sm[-1L, 1L]), or_low = exp(ci[-1L, 1L]), or_high = exp(ci[-1L, 2L]),
                        p = sm[-1L, 4L], n_people = nrow(d2), n_missing = sum(d2$device_missing), note = "")
    model <- rbind(model, data.table(term = "cohort (joint likelihood-ratio test)", odds_ratio = NA_real_, or_low = NA_real_, or_high = NA_real_,
                                     p = lrt$`Pr(>Chi)`[2L], n_people = nrow(d2), n_missing = sum(d2$device_missing),
                                     note = sprintf("chisq(%d) = %.2f", lrt$Df[2L], lrt$Deviance[2L])))
    model[, p_holm := fifelse(grepl("^cohort", term) & !grepl("joint", term), stats::p.adjust(p, "holm"), NA_real_)]
    p_holm_idx <- grepl("^cohort", model$term) & !grepl("joint", model$term)
    model[p_holm_idx, p_holm := stats::p.adjust(p, "holm")]
  }
  excluded <- setdiff(lv, usable)
  if (length(excluded)) {
    ex <- ev[cohort %in% excluded, .(term = paste0("cohort", cohort), odds_ratio = NA_real_, or_low = NA_real_, or_high = NA_real_, p = NA_real_,
                                      n_people = n, n_missing = e, note = sprintf("not estimable: %d of %d missing (see the Fisher exact test in device_missing_vs_reference)", e, n),
                                      p_holm = NA_real_)]
    model <- rbind(model, ex, fill = TRUE)
  }

  # ---- (d) composition of the main device family ----
  known <- p[device_missing == FALSE]
  fam_order <- known[, .N, by = device_family][order(-N), device_family]
  fam_order <- c(setdiff(fam_order, cfg$device$other_label), intersect(fam_order, cfg$device$other_label))
  cats <- c(miss_lab, fam_order)
  p[, category := factor(device_family, levels = cats)]
  ct <- table(p$cohort, p$category)[lv, , drop = FALSE]
  comp <- as.data.table(as.data.frame(ct, stringsAsFactors = FALSE))
  setnames(comp, c("cohort", "category", "n"))
  comp[, N := sum(n), by = cohort]
  cc <- wilson_ci(comp$n, comp$N)
  comp[, `:=`(pct = 100 * cc$est, pct_low = 100 * cc$low, pct_high = 100 * cc$high)]
  comp[, category := factor(category, levels = cats)]
  comp <- comp[order(match(cohort, cfg$cohorts), category)]

  tab_known <- ct[, setdiff(colnames(ct), miss_lab), drop = FALSE]
  tests <- rbind(tests,
                 cbind(test = "Composition differs between cohorts, missing included (chi-square, cohort x category)", chisq_table(ct, cfg$device$mc_B, cfg$device$boot_B)),
                 cbind(test = "Composition differs between cohorts, among people with a known device", chisq_table(tab_known, cfg$device$mc_B, cfg$device$boot_B)))
  tests[, `:=`(p_monte_carlo_holm = stats::p.adjust(p_monte_carlo, "holm"))]

  # adjusted standardized residuals (Haberman): (O - E) / sqrt(E (1 - row share)(1 - column share))
  nn <- sum(ct); rs <- rowSums(ct); csum <- colSums(ct)
  E <- outer(rs, csum) / nn
  Z <- (ct - E) / sqrt(E * outer(1 - rs / nn, 1 - csum / nn))
  resid <- as.data.table(as.data.frame(as.table(Z), stringsAsFactors = FALSE))
  setnames(resid, c("cohort", "category", "adj_residual"))
  resid[, `:=`(observed = as.vector(ct), expected = as.vector(E))]
  resid[, p := 2 * stats::pnorm(-abs(adj_residual))]
  resid[, p_holm := stats::p.adjust(p, "holm")]
  resid[, direction := fifelse(p_holm < 0.05, fifelse(adj_residual > 0, "more than expected", "fewer than expected"), "")]
  resid <- resid[, .(cohort, category, observed, expected = round(expected, 1), adj_residual = round(adj_residual, 2), p, p_holm, direction)]

  # ---- (e) full version detail, (f) people with more than one wearable version ----
  ver <- known[, .N, by = .(cohort, device_family, device_version = device_version_main)]
  ver[, N_known := sum(N), by = cohort]
  vc <- wilson_ci(ver$N, ver$N_known)
  ver[, `:=`(pct_of_known = 100 * vc$est, pct_low = 100 * vc$low, pct_high = 100 * vc$high)]
  setnames(ver, "N", "n_people")
  ver <- ver[order(match(cohort, cfg$cohorts), -n_people)]
  sw <- known[, .(n_known = .N, n_two_or_more = sum(n_wearable_versions >= 2L)), by = cohort][order(match(cohort, cfg$cohorts))]
  sc <- wilson_ci(sw$n_two_or_more, sw$n_known)
  sw[, `:=`(pct_two_or_more = 100 * sc$est, pct_low = 100 * sc$low, pct_high = 100 * sc$high)]

  # the same cohort x category table read the other way: the cohort mix of each category (column percentages), with the mix of all people on top
  mix <- as.data.table(as.data.frame(ct, stringsAsFactors = FALSE)); setnames(mix, c("cohort", "category", "n"))
  tot <- mix[, .(n = sum(n)), by = cohort][, category := "All people"]
  mix <- rbind(tot, mix)
  mix[, n_category := sum(n), by = category]
  mxc <- wilson_ci(mix$n, mix$n_category)
  mix[, `:=`(pct_of_category = 100 * mxc$est, pct_low = 100 * mxc$low, pct_high = 100 * mxc$high)]
  mix[, `:=`(cohort = factor(cohort, levels = cfg$cohorts), category = factor(category, levels = c("All people", cats)))]
  mix <- mix[order(category, cohort)]

  list(device_missing_by_cohort = mb[], device_missing_vs_reference = pair[], device_missing_model = model, device_cohort_mix = mix[],
       device_composition_by_cohort = comp[], device_composition_residuals = resid[], device_tests = tests[],
       device_version_by_cohort = ver[], device_switching_by_cohort = sw[])
}
