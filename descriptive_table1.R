#################################################################################
# descritibe_table1.R
# summarize tabuleer data and generate descritive statistics in a table
#
#################################################################################

library(readxl)

work.dir <- '/Users/zqiu/Documents/graduate_school/Bui Lab/AllOfUs/scripts/'
path.in <- '/Users/zqiu/Documents/graduate_school/Bui Lab/AllOfUs/'
path.out <- '/Users/zqiu/Documents/graduate_school/Bui Lab/AllOfUs/summary/'

setwd(work.dir)

generate_table1 <- function(
    data,
    group_var,
    vars = list(continuous  = c("age", "BMI"),
                        categorical = c("sex", "race", "ethnicity", "education")),
    var_labels   = NULL,
    digits_cont  = 1,
    digits_pct   = 1,
    p_digits     = 3,
    overall      = TRUE,
    missing_val  = c("show", "omit")
) {
    missing_val <- match.arg(missing_val)

    # ── 0. Input validation ────────────────────────────────────────────────────
    stopifnot(
        is.data.frame(data),
        is.character(group_var), length(group_var) == 1,
        group_var %in% names(data)
    )

    grp        <- factor(data[[group_var]])
    grp_levels <- levels(grp)
    if (length(grp_levels) != 2)
        stop("`group_var` must have exactly 2 levels. Found: ",
            paste(grp_levels, collapse = ", "))

    g1 <- grp_levels[1]
    g2 <- grp_levels[2]

    # ── Helper: display label ──────────────────────────────────────────────────
    label_for <- function(v) {
        if (!is.null(var_labels) && v %in% names(var_labels)) var_labels[[v]] else v
    }

    # ── Helper: format p-value ─────────────────────────────────────────────────
    fmt_p <- function(p) {
        if (is.na(p)) return("NA")
        if (p < 0.001) return("<0.001")
        formatC(p, digits = p_digits, format = "f")
    }
 
    # ── Helper: normality check (tiered by n) ──────────────────────────────────
    #
    #  n < 10          : S-W has no power → always non-parametric
    #  10  ≤ n < 30    : Visual zone — S-W used but visual inspection recommended
    #  30  ≤ n < 200   : S-W reliable; use p > 0.05
    #  200 ≤ n < 5000  : S-W over-sensitive; switch to skewness/kurtosis heuristic
    #  n ≥ 5000        : CLT guarantees normality of means → always parametric
    #
    is_normal <- function(x) {
        x <- x[!is.na(x)]
        n <- length(x)

        # Zone 1: too few observations — S-W has essentially no power
        if (n < 10)  return(FALSE)

        # Zone 2 & 3: S-W reliable range (10–199)
        if (n < 200) {
        return(tryCatch(shapiro.test(x)$p.value > 0.05, error = function(e) FALSE))
        }

        # Zone 4: S-W overpowered (200–4999) — trivial deviations cause rejection
        # Use distributional shape heuristic instead:
        #   |skewness| < 1  AND  |excess kurtosis| < 1
        if (n < 5000) {
        mu   <- mean(x); s <- sd(x)
        skew <- mean(((x - mu) / s)^3)
        kurt <- mean(((x - mu) / s)^4) - 3
        return(abs(skew) < 1 && abs(kurt) < 1)
        }

        # Zone 5: n ≥ 5000 — CLT guarantees approximate normality of means;
        # non-parametric tests also lose their robustness advantage at large n.
        return(TRUE)
    }
 
  # ── Helper: continuous row ─────────────────────────────────────────────────

    make_cont_row <- function(v) {
        x_all <- data[[v]]
        x1    <- x_all[grp == g1]
        x2    <- x_all[grp == g2]
    
        n_miss_all <- sum(is.na(x_all))
        n_miss_1   <- sum(is.na(x1))
        n_miss_2   <- sum(is.na(x2))
    
        norm1 <- is_normal(x1[!is.na(x1)])
        norm2 <- is_normal(x2[!is.na(x2)])
        use_parametric <- norm1 && norm2
    
        if (use_parametric) {
        # Mean ± SD
        fmt_cont <- function(x) {
            x <- x[!is.na(x)]
            if (length(x) == 0) return("—")
            sprintf("%s ± %s",
                    formatC(mean(x), digits = digits_cont, format = "f"),
                    formatC(sd(x),   digits = digits_cont, format = "f"))
        }
        summary_label <- "Mean ± SD"
        test_name     <- "t-test (Welch)"
        p_val <- tryCatch(
            t.test(x1, x2, var.equal = FALSE)$p.value,
            error = function(e) NA_real_
        )
        } else {
            # Median [Q1–Q3]
            fmt_cont <- function(x) {
                x <- x[!is.na(x)]
                if (length(x) == 0) return("—")
                q  <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
                sprintf("%s [%s–%s]",
                    formatC(q[2], digits = digits_cont, format = "f"),
                    formatC(q[1], digits = digits_cont, format = "f"),
                    formatC(q[3], digits = digits_cont, format = "f")
                )
            }
            summary_label <- "Median [Q1–Q3]"
            test_name     <- "Wilcoxon rank-sum"
            p_val <- tryCatch(
            wilcox.test(x1, x2, exact = FALSE)$p.value,
            error = function(e) NA_real_
            )
        }
    
        miss_str <- if (missing_val == "show" && n_miss_all > 0)
        sprintf(" [%d missing]", n_miss_all) else ""
    
        row <- data.frame(
            Variable    = paste0(label_for(v), miss_str),
            Summary     = summary_label,
            Overall     = fmt_cont(x_all),
            G1          = fmt_cont(x1),
            G2          = fmt_cont(x2),
            P_value     = fmt_p(p_val),
            Test        = test_name,
            stringsAsFactors = FALSE
        )
        return(row)
    }
 
    # ── Helper: categorical row ────────────────────────────────────────────────
    make_cat_rows <- function(v) {
        x_all <- data[[v]]
    
        # Coerce to factor — preserve existing levels if already a factor
        if (!is.factor(x_all)) x_all <- factor(x_all)
        lvls <- levels(x_all)
    
        n_miss_all <- sum(is.na(x_all))
        x1  <- x_all[grp == g1]
        x2  <- x_all[grp == g2]
    
        n_all <- sum(!is.na(x_all))
        n1    <- sum(!is.na(x1))
        n2    <- sum(!is.na(x2))
    
        # Contingency table (exclude NA)
        tbl <- table(x_all[!is.na(x_all)], grp[!is.na(x_all)])
    
        # Choose test: Fisher if any expected cell < 5, else chi-square
        exp_cells <- tryCatch(
            suppressWarnings(chisq.test(tbl)$expected),
            error = function(e) matrix(0)
        )
        use_fisher <- any(exp_cells < 5, na.rm = TRUE) || ncol(tbl) < 2
    
        if (use_fisher) {
            p_val     <- tryCatch(fisher.test(tbl, simulate.p.value = TRUE,
                                                B = 2000)$p.value,
                                    error = function(e) NA_real_)
            test_name <- "Fisher's exact"
        } else {
            p_val     <- tryCatch(chisq.test(tbl)$p.value,
                                    error = function(e) NA_real_)
            test_name <- "Chi-square"
        }
    
        miss_str <- if (missing_val == "show" && n_miss_all > 0)
        sprintf(" [%d missing]", n_miss_all) else ""
    
        # Header row
        header <- data.frame(
            Variable = paste0(label_for(v), miss_str, ", n (%)"),
            Summary  = "",
            Overall  = "",
            G1       = "",
            G2       = "",
            P_value  = fmt_p(p_val),
            Test     = test_name,
            stringsAsFactors = FALSE
        )
    
        fmt_pct <- function(k, n) {
            if (n == 0) return("—")
            sprintf("%d (%s%%)", k,
                formatC(100 * k / n, digits = digits_pct, format = "f"))
        }
    
        # One row per level
        level_rows <- lapply(lvls, function(lv) {
        k_all <- sum(x_all == lv, na.rm = TRUE)
        k1    <- sum(x1    == lv, na.rm = TRUE)
        k2    <- sum(x2    == lv, na.rm = TRUE)
        data.frame(
            Variable = paste0("  ", lv),
            Summary  = "",
            Overall  = fmt_pct(k_all, n_all),
            G1       = fmt_pct(k1, n1),
            G2       = fmt_pct(k2, n2),
            P_value  = "",
            Test     = "",
            stringsAsFactors = FALSE
        )
        }
        )
    
        do.call(rbind, c(list(header), level_rows))
  }
 
    # ── 1. N row ────────────────────────────────────────────────────────────────
    n_all <- nrow(data)
    n1    <- sum(grp == g1, na.rm = TRUE)
    n2    <- sum(grp == g2, na.rm = TRUE)

    n_row <- data.frame(
    Variable = "N",
    Summary  = "",
    Overall  = as.character(n_all),
    G1       = as.character(n1),
    G2       = as.character(n2),
    P_value  = "",
    Test     = "",
    stringsAsFactors = FALSE
    )

    # ── 2. Build continuous rows ───────────────────────────────────────────────
    cont_rows <- lapply(vars$continuous, make_cont_row)
    cont_rows <- do.call(rbind, cont_rows)

    # ── 3. Build categorical rows ──────────────────────────────────────────────
    cat_rows <- lapply(vars$categorical, make_cat_rows)
    cat_rows <- do.call(rbind, cat_rows)

    # ── 4. Assemble final table ────────────────────────────────────────────────
    tbl <- rbind(n_row, cont_rows, cat_rows)
    rownames(tbl) <- NULL

    # Rename group columns to actual group names
    names(tbl)[names(tbl) == "G1"] <- g1
    names(tbl)[names(tbl) == "G2"] <- g2

    if (!overall) tbl$Overall <- NULL

    # Re-order: Variable | Summary | [Overall] | G1 | G2 | P_value | Test
    return(tbl)
}

dob_to_age <- function(dob, ref_date = Sys.Date()) {
    dob <- as.Date(dob)           # handles "YYYY-MM-DD" strings
    ref_date <- as.Date(ref_date)

    age_raw <- as.integer(format(ref_date, "%Y")) -
                as.integer(format(dob,      "%Y"))

    # Subtract 1 if the birthday hasn't occurred yet this calendar year
    had_birthday <- (
        as.integer(format(ref_date, "%m")) >  as.integer(format(dob, "%m"))
    ) | (
        as.integer(format(ref_date, "%m")) == as.integer(format(dob, "%m")) &
        as.integer(format(ref_date, "%d")) >= as.integer(format(dob, "%d"))
    )

    return(age_raw - as.integer(!had_birthday))
}


#### read data #####################
df <- read_excel(
    paste0(path.in, 'fitbit_meta.xlsx')
);


# get age
df$age <- dob_to_age(df$date_of_birth)

# remove duplicated
df <- df[!duplicated(df$person_id),]

tbl <- generate_table1(
    data        = df,
    group_var   = "disease",
    vars        = list(
        continuous  = c("age"),
        categorical = c("sex_at_birth", "race", "ethnicity")
    ),
    var_labels  = labels,
    digits_cont = 1,
    digits_pct  = 1,
    p_digits    = 3,
    overall     = TRUE,
    missing_val = "show"
)


# ── Pretty print ────────────────────────────────────────────────────────────
cat("\n", paste(rep("=", 100), collapse = ""), "\n")
cat("TABLE 1. Baseline Characteristics by Disease Group\n")
cat(paste(rep("=", 100), collapse = ""), "\n\n")
 
# Determine column widths
col_widths <- sapply(names(tbl1), function(col) max(nchar(c(col, tbl1[[col]])), na.rm = TRUE) + 2)
fmt_row <- function(row_vals) {
  paste(mapply(function(v, w) formatC(v, width = -w, flag = "-"), row_vals, col_widths), collapse = "")
}
 
cat(fmt_row(names(tbl1)), "\n")
cat(paste(rep("-", sum(col_widths)), collapse = ""), "\n")
for (i in seq_len(nrow(tbl1))) {
  cat(fmt_row(as.character(tbl1[i, ])), "\n")
}
cat(paste(rep("=", 100), collapse = ""), "\n\n")
 
# Optional: kable output (uncomment if knitr is available)
# if (requireNamespace("knitr", quietly = TRUE)) {
#   print(knitr::kable(tbl1, format = "markdown",
#                      caption = "Table 1. Baseline Characteristics by Disease Group"))
# }
 
cat("Table 1 generated successfully.\n")
cat("Object 'tbl1' is a plain data.frame — pass to knitr::kable(), gt::gt(),\n")
cat("or flextable::flextable() for further formatting.\n\n")

