library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(patchwork)   
library(scales)
library(rstatix)   # install.packages("rstatix")
library(ggpubr)    # install.packages("ggpubr")


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# ac_active_plot_func
# get comparing for activity for each cohort
# Violin plots; CV distribution; Mean vs SD scatter of each variable
# ───────────────────────────────────────────────────────────────────────────────────────────────────
ac_active_plot_func <- function(df, path.out) {
  cohort_colors <- c(
    "Diabetes"       = "#E74C3C",
    "Prostate" = "#3498DB",
    "Insomnia"       = "#F39C12",
    "Healthy"        = "#27AE60"
  )
  
  # ── Plot 1: Activity minutes by cohort ───────────────────────────────────────
  p_v1_data <- df |>
    select(person_id, cohort,
           sedentary_minutes__mean, lightly_active_minutes__mean,
           fairly_active_minutes__mean, very_active_minutes__mean) |>
    pivot_longer(ends_with("__mean"), names_to = "metric", values_to = "mean_val") |>
    mutate(
      metric = str_remove(metric, "__mean"),
      metric = str_replace_all(metric, "_", " ") |> str_to_title(),
      cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy", "Prostate")) # "Healthy", 
    )
  summary(p_v1_data)
  # Kruskal-Wallis p-value per metric
  kw_v1 <- p_v1_data |>
    group_by(metric) |>
    kruskal_test(mean_val ~ cohort) |>
    mutate(label = paste0("KW p", ifelse(p < 0.001, "<0.001", paste0("=", round(p, 3)))))
  
  
  p_v1 <- ggplot(p_v1_data, aes(x = cohort, y = mean_val, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.6,
                 fill = "white", color = "gray30") +
    geom_text(data = kw_v1, aes(x = 2.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.5, size = 3, fontface = "italic") +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = cohort_colors) +
    labs(title = "Patient-Level Mean Activity Minutes by Cohort",
         x = NULL, y = "Mean minutes/day", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  ggsave(p_v1,      filename = paste0(path.out, "activity_zoneminutes_violin.png"), width = 9, height = 8, dpi = 150)
  
  # ── Plot 2: Steps by cohort ───────────────────────────────────────────────────
  p_v2_data <- df |>
    filter(!is.na(cohort)) |>
    select(person_id, cohort, steps__mean) |>
    rename(mean_val = steps__mean) |>
    mutate(cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy","Prostate")))
  
  step_upper <- p_v2_data |>
    group_by(cohort) |>
    summarise(upper = quantile(mean_val, 0.75, na.rm = TRUE) +
                3 * IQR(mean_val, na.rm = TRUE)) |>
    pull(upper) |> max()
  
  outlier_counts <- p_v2_data |>
    group_by(cohort) |>
    summarise(
      n_out  = sum(mean_val > step_upper, na.rm = TRUE),
      label  = ifelse(n_out > 0, paste0("n=", n_out, " above"), ""),
      y_pos  = step_upper * 0.98
    )
  
  # ── KW test ───────────────────────────────────────────────────────────────────
  kw_v2 <- p_v2_data |> kruskal_test(mean_val ~ cohort)
  kw_label <- paste0("Kruskal-Wallis p",
                     ifelse(kw_v2$p < 0.001, "<0.001",
                            paste0("=", round(kw_v2$p, 3))))
  dunn_v2 <- p_v2_data |>
    dunn_test(mean_val ~ cohort, p.adjust.method = "bonferroni") |>
    add_xy_position(x = "cohort")
  # ── Plot ──────────────────────────────────────────────────────────────────────
  p_v2 <- ggplot(p_v2_data, aes(x = cohort, y = mean_val, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8,
                data = p_v2_data |> filter(mean_val <= step_upper)) +
    geom_boxplot(width = 0.1, outlier.shape = NA,
                 fill = "white", color = "gray30") +
    geom_hline(yintercept = 10000, linetype = "dashed",
               color = "gray40", linewidth = 0.6) +
    geom_text(data = outlier_counts,
              aes(x = cohort, y = y_pos, label = label),
              inherit.aes = FALSE, size = 3, color = "#C0392B",
              vjust = 1, fontface = "italic") +
    annotate("text", x = 4.4, y = 10000, vjust = -0.5,
             hjust = 1, size = 3, color = "gray40",
             label = "10,000 step target") +
    annotate("text", x = 3.9, y = step_upper * 0.05,
             hjust = 0, vjust = 0, size = 2.2, fontface = "italic",
             color = "gray40", label = kw_label) +
    stat_pvalue_manual(dunn_v2, label = "p.adj.signif", tip.length = 0.01,
                       hide.ns = TRUE, size = 3, y.position = c(18000, 21000, 23000, 20000, 19000, 22000)) +
    coord_cartesian(ylim = c(0, step_upper)) +
    scale_fill_manual(values = cohort_colors) +
    labs(title    = "Patient-Level Mean Daily Steps by Cohort",
         x = NULL, y = "Mean steps/day", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  
  ggsave(p_v2,      filename = paste0(path.out, "activity_steps_violin.png"),       width = 7, height = 5, dpi = 150)
  
  # ── CV plot by cohort ─────────────────────────────────────────────────────────
  key_metrics <- c("steps", "calories_out", "sedentary_minutes",
                   "lightly_active_minutes", "fairly_active_minutes",
                   "very_active_minutes")
  p_cv <- df |>
    select(person_id, cohort, all_of(paste0(key_metrics, "__cv"))) |>
    pivot_longer(ends_with("__cv"), names_to = "metric", values_to = "cv") |>
    mutate(
      metric = str_remove(metric, "__cv"),
      metric = str_replace_all(metric, "_", " ") |> str_to_title(),
      cohort = factor(cohort, levels = c("Healthy", "Insomnia", "Diabetes", "Prostate"))
    ) |>
    filter(!is.na(cv), cv < 200) |>
    ggplot(aes(x = cohort, y = cv, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.4,
                 fill = "white", color = "gray30") +
    geom_hline(yintercept = 30, linetype = "dashed",
               color = "#C0392B", linewidth = 0.5) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = cohort_colors) +
    labs(title    = "Within-Person Variability (CV%) by Cohort",
         subtitle = "Red dashed = 30% CV reference",
         x = NULL, y = "CV (%)", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  ggsave(p_cv, filename = paste0(path.out, "activity_cv_violin.png"),width = 9, height = 8, dpi = 150)
  
  
  # ── Mean vs SD scatter by cohort ──────────────────────────────────────────────
  p_mean_sd <- df |>
    select(person_id, cohort, n_obs,
           ends_with("__mean"), ends_with("__sd")) |>
    pivot_longer(-c(person_id, cohort, n_obs),
                 names_to  = c("metric", ".value"),
                 names_sep = "__") |>
    mutate(metric = str_replace_all(metric, "_", " ") |> str_to_title(),
           cohort = factor(cohort, levels = c("Healthy", "Insomnia",
                                              "Diabetes", "Prostate"))) |>
    filter(!is.na(mean), !is.na(sd)) |>
    ggplot(aes(x = mean, y = sd, color = cohort)) +
    geom_point(size = 0.4, alpha = 0.1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
    facet_wrap(~metric, scales = "free", ncol = 3) +
    scale_color_manual(values = cohort_colors) +
    labs(title  = "Patient Mean vs. SD by Cohort (Heteroscedasticity Check)",
         x = "Patient mean", y = "Patient SD", color = "Cohort") +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          legend.position = "bottom")
  
  ggsave(p_mean_sd, filename = paste0(path.out, "activity_mean_sd.png"), width = 10, height = 7, dpi = 150)
  
  return()
  }


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# ac_active_plot_bysex_func
# Violin plots by sex calories out
# ───────────────────────────────────────────────────────────────────────────────────────────────────
ac_active_plot_bysex_func <- function(df, path.out) {
  # ── Calories out by cohort and sex ────────────────────────────────────────────
  p_v3_data <- df |>
    mutate(
      sex = case_when(
        sex_at_birth == "Female" ~ "Female",
        sex_at_birth == "Male"   ~ "Male",
        TRUE                     ~ "Other"
      ),
      sex    = factor(sex, levels = c("Female", "Male", "Other")),
      cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy", "Prostate"))
    )
  summary_tbl <- p_v3_data |>
    group_by(cohort, sex) |>
    summarise(
      N      = n(),
      Mean   = round(mean(calories_out__mean,             na.rm = TRUE), 1),
      SD     = round(sd(calories_out__mean,               na.rm = TRUE), 1),
      Median = round(median(calories_out__mean,           na.rm = TRUE), 1),
      Q1     = round(quantile(calories_out__mean, 0.25,   na.rm = TRUE), 1),
      Q3     = round(quantile(calories_out__mean, 0.75,   na.rm = TRUE), 1),
      .groups = "drop"
    ) |>
    arrange(cohort, sex)
  print(summary_tbl)
  kw_v3 <- p_v3_data |>
    group_by(sex) |>
    filter(n_distinct(cohort) > 1) |>
    kruskal_test(calories_out__mean ~ cohort)
  
  dodge <- position_dodge(width = 0.9)
  
  p_v3 <- ggplot(p_v3_data, aes(x = cohort, y = calories_out__mean, fill = sex)) +
    geom_violin(color = NA, alpha = 0.8, position = dodge) +
    geom_boxplot(width = 0.08, outlier.size = 0.3, outlier.alpha = 0.6,
                 color = "gray30", position = dodge) +
    geom_text(data = kw_v3,
              aes(x = 3.8, y = Inf,
                  label = paste0("KW p", ifelse(p < 0.001, "<0.001",
                                                paste0("=", round(p, 3))),
                                 " (", sex, ")")),
              inherit.aes = FALSE,
              vjust = seq(1.5, 4, length.out = nrow(kw_v3)),
              size = 2.8, fontface = "italic") +
    scale_fill_manual(values = c("Female" = "#E91E8C",
                                 "Male"   = "#3498DB",
                                 "Other"  = "#95A5A6")) +
    labs(title    = "Patient-Level Mean Daily Calories Out by Cohort and Sex",
         subtitle = "Kruskal-Wallis test between cohorts within each sex group",
         x = NULL, y = "Mean calories/day", fill = "Sex at birth") +
    theme_bw(base_size = 11) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  
    ggsave(p_v3, filename = paste0(path.out, "activity_caloriesout_violin_bysex.png"), width = 7, height = 6, dpi = 150)
return()
}


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# ac_sleep_plot_func
# get comparing plots for sleep for each cohort
# Violin plots; CV distribution; Mean vs SD scatter of each variable
# ───────────────────────────────────────────────────────────────────────────────────────────────────
ac_sleep_plot_func <- function(df, path.out) {
  cohort_colors <- c(
    "Diabetes"       = "#E74C3C",
    "Prostate" = "#3498DB",
    "Insomnia"       = "#F39C12",
    "Healthy"        = "#27AE60"
  )
  
  # ── Plot 1: Activity minutes by cohort ───────────────────────────────────────
  p_v1_data <- df |>
    select(person_id, cohort,
           minute_in_bed__mean, minute_asleep__mean, minute_awake__mean, minute_deep__mean, minute_light__mean, minute_rem__mean, minute_wake__mean) |>
    pivot_longer(ends_with("__mean"), names_to = "metric", values_to = "mean_val") |>
    mutate(
      metric = str_remove(metric, "__mean"),
      metric = str_replace_all(metric, "_", " ") |> str_to_title(),
      cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy", "Prostate")) # "Healthy", 
    )
  summary(p_v1_data)
  # Kruskal-Wallis p-value per metric
  kw_v1 <- p_v1_data |>
    group_by(metric) |>
    kruskal_test(mean_val ~ cohort) |>
    mutate(label = paste0("KW p", ifelse(p < 0.001, "<0.001", paste0("=", round(p, 3)))))
  
  
  p_v1 <- ggplot(p_v1_data, aes(x = cohort, y = mean_val, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.6,
                 fill = "white", color = "gray30") +
    geom_text(data = kw_v1, aes(x = 2.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.5, size = 3, fontface = "italic") +
    facet_wrap(~metric, scales = "free_y", ncol = 7) +
    scale_fill_manual(values = cohort_colors) +
    labs(title = "Patient-Level Mean Sleep Minutes by Cohort",
         x = NULL, y = "Mean minutes/day", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  ggsave(p_v1,      filename = paste0(path.out, "sleep_stage_violin.png"), width = 12, height = 6, dpi = 150)
  
  ### sleep efficiency
  p_v2_data <- df |>
    filter(!is.na(cohort)) |>
    select(person_id, cohort, sleep_efficiency__mean) |>
    rename(mean_val = sleep_efficiency__mean) |>
    mutate(cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy","Prostate")))
  
  sleep_efficiency_upper <- p_v2_data |>
    group_by(cohort) |>
    summarise(upper = quantile(mean_val, 0.75, na.rm = TRUE) +
                3 * IQR(mean_val, na.rm = TRUE)) |>
    pull(upper) |> max()
  print(sleep_efficiency_upper)
  
  outlier_counts <- p_v2_data |>
    group_by(cohort) |>
    summarise(
      n_out  = sum(mean_val > sleep_efficiency_upper, na.rm = TRUE),
      label  = ifelse(n_out > 0, paste0("n=", n_out, " above"), ""),
      y_pos  = sleep_efficiency_upper * 0.98
    )
  print(outlier_counts)
  # ── KW test ───────────────────────────────────────────────────────────────────
  kw_v2 <- p_v2_data |> kruskal_test(mean_val ~ cohort)
  kw_label <- paste0("Kruskal-Wallis p",
                     ifelse(kw_v2$p < 0.001, "<0.001",
                            paste0("=", round(kw_v2$p, 3))))
  dunn_v2 <- p_v2_data |>
    dunn_test(mean_val ~ cohort, p.adjust.method = "bonferroni") |>
    add_xy_position(x = "cohort")
  # ── Plot ──────────────────────────────────────────────────────────────────────
  p_v2 <- ggplot(p_v2_data, aes(x = cohort, y = mean_val, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8,
                data = p_v2_data |> filter(mean_val <= sleep_efficiency_upper)) +
    geom_boxplot(width = 0.1, outlier.shape = NA,
                 fill = "white", color = "gray30") +
    geom_text(data = outlier_counts,
              aes(x = cohort, y = y_pos, label = label),
              inherit.aes = FALSE, size = 3, color = "#C0392B",
              vjust = 1, fontface = "italic") +
    annotate("text", x = 3.9, y = sleep_efficiency_upper * 1.04,
             hjust = 0, vjust = 0, size = 2.2, fontface = "italic",
             label = kw_label) +
    stat_pvalue_manual(dunn_v2, label = "p.adj.signif", tip.length = 0.01,
                       hide.ns = TRUE, size = 3, y.position = c(100, 103, 102, 101, 104)) +
    coord_cartesian(ylim = c(65, sleep_efficiency_upper * 1.05)) +
    scale_fill_manual(values = cohort_colors) +
    labs(title    = "Patient-Level Mean Daily Sleep Efficiency by Cohort",
         x = NULL, y = "Mean sleep efficiency/day", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  
  ggsave(
    p_v2,
    filename = paste0(path.out, "sleep_efficiency_violin.png"),
    width = 7, height = 9, dpi = 150)
  
  return()
}

# ───────────────────────────────────────────────────────────────────────────────────────────────────
# final_demographic_plot_condition_func
# plots for conditions for final df
# curves, barplot, 
# ───────────────────────────────────────────────────────────────────────────────────────────────────
final_demographic_plot_cii_func <- function(df, path.out) {
  # ── CCI category labels ───────────────────────────────────────────────────────
  cci_categorized <- df %>%
    mutate(cci_category = case_when(
      cci_score == 0            ~ "None (0)",
      cci_score <= 2            ~ "Mild (1-2)",
      cci_score <= 4            ~ "Moderate (3-4)",
      cci_score >= 5            ~ "Severe (≥5)"
    )) %>%
    mutate(cci_category = factor(cci_category,
                                 levels = c("None (0)", "Mild (1-2)",
                                            "Moderate (3-4)", "Severe (≥5)")))
  
  # ── Summary table: median CCI + n and % per category within each cohort ───────
  cci_summary <- cci_categorized %>%
    group_by(cohort) %>%
    mutate(
      n_cohort   = n(),
      median_cci = median(cci_score)
    ) %>%
    group_by(cohort, cci_category, n_cohort, median_cci) %>%
    summarise(n_category = n(), .groups = "drop") %>%
    mutate(pct = round(n_category / n_cohort * 100, 1))
  
  print(cci_summary)
  
  cci_categorized2 <- df %>%
    mutate(cci_category = case_when(
      cci_score <= 7            ~ "<=7",
      cci_score == 8            ~ "8",
      cci_score == 9            ~ "9",
      cci_score == 10            ~ "10",
      cci_score >= 11            ~ ">=11"
    )) %>%
    mutate(cci_category = factor(cci_category,
                                 levels = c("<=7", "8",
                                            "9", "10", ">=11")))
  
  # ── Summary table: median CCI + n and % per category within each cohort ───────
  cci_summary2 <- cci_categorized2 %>%
    group_by(cohort) %>%
    mutate(
      n_cohort   = n(),
      median_cci = median(cci_score)
    ) %>%
    group_by(cohort, cci_category, n_cohort, median_cci) %>%
    summarise(n_category = n(), .groups = "drop") %>%
    mutate(pct = round(n_category / n_cohort * 100, 1))
  
  print(cci_summary2)
  
  # ── Median CCI per cohort ─────────────────────────────────────────────────────
  df %>%
    group_by(cohort) %>%
    summarise(
      n          = n(),
      median_cci = median(cci_score),
      mean_cci   = round(mean(cci_score), 2),
      sd_cci     = round(sd(cci_score), 2),
      .groups    = "drop"
    )
  
  # ── PLOT 1: Stacked bar — % of patients in each CCI category by cohort ────────
  category_colors <- c(
    "None (0)"       = "#27AE60",
    "Mild (1-2)"     = "#F39C12",
    "Moderate (3-4)" = "#E67E22",
    "Severe (≥5)"    = "#E74C3C"
  )
  
  p_bar <- cci_summary %>%
    ggplot(aes(x = cohort, y = pct, fill = cci_category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(pct, "%\n(n=", n_category, ")")),
              position = position_stack(vjust = 0.5),
              size = 3.2, color = "white", fontface = "bold") +
    scale_fill_manual(values = category_colors) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "CCI Category Distribution by Cohort",
      x     = NULL,
      y     = "% of Patients",
      fill  = "CCI Category"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(p_bar, filename = paste0(path.out, "final_cci_stackbarplot.png"), width = 14, height = 8, dpi = 150)
  
  
  # ── PLOT 2: Density curve of CCI score by cohort ──────────────────────────────
  # Add median lines per cohort
  cohort_medians <- df %>%
    group_by(cohort) %>%
    summarise(median_cci = median(cci_score), .groups = "drop")
  
  p_density <- df %>%
    ggplot(aes(x = cci_score, fill = cohort, color = cohort)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(data = cohort_medians,
               aes(xintercept = median_cci, color = cohort),
               linetype = "dashed", linewidth = 0.8) +
    geom_text(data = cohort_medians,
              aes(x = median_cci, y = Inf,
                  label = paste0("median=", median_cci),
                  color = cohort),
              vjust = 1.5, hjust = -0.1, size = 3.2, show.legend = FALSE) +
    scale_fill_manual(values  = cohort_colors) +
    scale_color_manual(values = cohort_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
      title = "CCI Score Distribution by Cohort",
      x     = "CCI Score (with age adjustment)",
      y     = "Density",
      fill  = "Cohort",
      color = "Cohort"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(p_density, filename = paste0(path.out, "final_cci_density.png"), width = 14, height = 8, dpi = 150)
  
  
  # ── PLOT 3: Grouped bar — n patients per category per cohort (raw counts) ──────
  p_grouped <- cci_summary %>%
    ggplot(aes(x = cci_category, y = n_category, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(aes(label = scales::comma(n_category)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 3) +
    scale_fill_manual(values = cohort_colors) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "Number of Patients per CCI Category by Cohort",
      x     = "CCI Category",
      y     = "Number of Patients",
      fill  = "Cohort"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(p_grouped, filename = paste0(path.out, "final_cci_barplot.png"), width = 14, height = 8, dpi = 150)
  
  
  return()
}

# ───────────────────────────────────────────────────────────────────────────────────────────────────
# final_demographic_plot_age_func
# plots for age for final df
# curves and barplot for older age
# ───────────────────────────────────────────────────────────────────────────────────────────────────
final_demographic_plot_age_func <- function(df, path.out) {
  cohort_colors <- c(
    "Diabetes"  = "#E74C3C",
    "Prostate"  = "#3498DB",
    "Insomnia"  = "#F39C12",
    "Healthy"   = "#27AE60"
  )

  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 3 — age distribution within each cohort
  # ════════════════════════════════════════════════════════════════════════════
  p3 <- df %>%
    distinct(cohort, person_id, age) %>%
    ggplot(aes(x = age, fill = cohort, color = cohort)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    scale_fill_manual(values  = cohort_colors) +
    scale_color_manual(values = cohort_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
      title = "Age Distribution by Cohort",
      x     = "Age (years)",
      y     = "Density",
      fill  = "Cohort",
      color = "Cohort"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  ggsave(p3, filename = paste0(path.out, "final_age_curve.png"), width = 14, height = 8, dpi = 150)
  
  p4 <- df %>%
    distinct(cohort, person_id, age) %>%
    mutate(age_group = case_when(
      age < 65              ~ "<65",
      age >= 65 & age < 70  ~ "65-70",
      age >= 70 & age < 75  ~ "70-75",
      age >= 75 & age < 80  ~ "75-80",
      age >= 80 & age < 85  ~ "80-85",
      age >= 85 & age < 90  ~ "85-90",
      age >= 90             ~ ">90"
    )) %>%
    mutate(age_group = factor(age_group, levels = c("<65", "65-70", "70-75", 
                                                    "75-80", "80-85", "85-90", ">90"))) %>%
    group_by(cohort, age_group) %>%
    summarise(n_patients = n_distinct(person_id), .groups = "drop") %>%
    ggplot(aes(x = age_group, y = n_patients, fill = cohort)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_fill_manual(values = cohort_colors) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = "Number of Patients by Age Group and Cohort",
      x     = "Age Group",
      y     = "Number of Patients",
      fill  = "Cohort"
    ) +
    theme_minimal(base_size = 13) + facet_wrap(~ cohort, scales = "free_y", nrow = 2) +
    theme(
      legend.position = "none",   # redundant with facet labels
      strip.text      = element_text(face = "bold")
    )
  ggsave(p4, filename = paste0(path.out, "final_age_barplot.png"), width = 14, height = 8, dpi = 150)
  
  return()
}



