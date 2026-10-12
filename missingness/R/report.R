# ==============================================================================
# report.R
# ------------------------------------------------------------------------------
# Writing results and the information needed to reproduce a run.
#   write_tables()          every data.table in a named list -> <outdir>/tables/<name>.tsv
#   write_run_metadata()    config, package versions, input fingerprints (md5 of L1), git commit
#   write_figure_captions() what each figure shows and how to read it
#   write_results_summary() the main numbers in plain language, from the saved tables only
#   unit_status_table()     which models ran, were resumed from a checkpoint, or failed in this call
# ==============================================================================

write_tables <- function(tables, outdir) {
  dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) {
    x <- tables[[nm]]
    if (is.data.frame(x) && nrow(x) > 0L) {
      fwrite(x, file.path(outdir, "tables", paste0(nm, ".tsv")), sep = "\t")
    }
  }
  invisible(names(tables))
}

# Input manifest: which L1 files (md5) or in-memory tables the run used, so a result can be tied to
# the exact input it came from. Works from the day store; nothing big is read for file-backed stores.
input_manifest <- function(store) {
  rbindlist(lapply(c(store$streams, "persons"), function(nm) {
    fp <- store$fingerprint(nm)
    data.table(input = nm, md5 = fp$md5 %||% NA_character_, bytes = fp$bytes %||% NA_real_,
               n_rows = fp$n_rows %||% NA_integer_, n_persons = fp$n_persons %||% NA_integer_)
  }))
}

git_commit <- function() {
  tryCatch(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1L]),
           error = function(e) NA_character_, warning = function(w) NA_character_)
}

write_run_metadata <- function(cfg, outdir, store = NULL, extra = list()) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(c(list(config = cfg), extra), file.path(outdir, "run_config.json"),
                         auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    saveRDS(c(list(config = cfg), extra), file.path(outdir, "run_config.rds"))
  }
  writeLines(c(paste("run time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
               paste("git commit:", git_commit()), "", capture.output(sessionInfo())),
             file.path(outdir, "session_info.txt"))
  if (!is.null(store)) fwrite(input_manifest(store), file.path(outdir, "input_manifest.tsv"), sep = "\t")
  invisible(outdir)
}

# ---- execution summary: what is done, what failed, what was skipped ---------------------------------------

# One line per figure: file name + what it shows and how to read it.
write_figure_captions <- function(names, outdir) {
  cap <- vapply(names, caption_for, character(1))
  dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("# Figure captions", "", paste0("* **", names, ".png** - ", ifelse(is.na(cap), "(no caption)", cap))),
             file.path(outdir, "figures", "figure_captions.md"))
  invisible(cap)
}

# Which units ran, were resumed from a checkpoint, or failed in this call (a data.table).
unit_status_table <- function(status_list) {
  if (length(status_list) == 0L) return(data.table(stream = character(), unit = character(), status = character()))
  rbindlist(lapply(names(status_list), function(s) data.table(stream = s, unit = names(status_list[[s]]), status = unname(status_list[[s]]))))
}

# Short plain-language summary of the main findings, written from the saved result tables only.
write_results_summary <- function(res, cfg, outdir) {
  L <- c("# Results summary", "", "Written from the saved result tables (no model is refitted).", "")
  fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
  if (!is.null(res$primary_test) && nrow(res$primary_test) > 0L) {
    L <- c(L, "## Gap onset (primary model, gap >= 1 day)",
           "Odds that a gap starts tomorrow per 1 SD below the person's usual steps yesterday (OR > 1 = lower activity before gaps).",
           res$primary_test[, sprintf("* %s: OR %s, p=%s, Holm p=%s [%s, %s events]", stream, fmt_or(or, or_low, or_high),
                                      signif(p, 3), signif(p_holm, 3), status, format(n_events, big.mark = ","))], "")
  }
  if (!is.null(res$sensitivity) && nrow(res$sensitivity) > 0L) {
    r <- res$sensitivity[!is.na(or), .(lo = min(or), hi = max(or), n = .N), by = stream]
    L <- c(L, "## Sensitivity analyses", "Range of the primary-term odds ratio over the estimable sensitivity analyses:",
           r[, sprintf("* %s: %.2f to %.2f over %d analyses", stream, lo, hi, n)], "")
  }
  B <- res$model_B_post_gap_activity
  if (!is.null(B) && nrow(B) > 0L) {
    L <- c(L, "## Model B (exploratory): activity on the first full day back, per doubling of gap length (SD units)",
           B[model_id == "B_back2_adj", sprintf("* %s: %.3f (%.3f to %.3f), p=%s [%s]", stream, estimate, conf_low, conf_high,
                                                signif(p, 3), status)], "")
  }
  C <- res$model_C_gap_length
  if (!is.null(C) && nrow(C) > 0L) {
    L <- c(L, "## Model C (exploratory): odds of a longer gap per 1 SD below usual steps",
           C[model_id == "C_adj", sprintf("* %s: %s [%s]", stream,
                                          ifelse(is.na(or), paste0("not estimated (", note, "); proportional-odds p=", signif(po_p, 2)), sprintf("OR %s", fmt_or(or, or_low, or_high))), status)], "")
  }
  fl <- rbindlist(list(res$onset_failures, res$sensitivity_failures, res$gap_failures), fill = TRUE)
  if (nrow(fl) > 0L)
    L <- c(L, "## FAILED units (re-run the same cell to retry them)", fl[, sprintf("* %s / %s / %s: %s", kind, stream, id, error)], "")
  writeLines(L, file.path(outdir, "results_summary.md"))
  invisible(L)
}
