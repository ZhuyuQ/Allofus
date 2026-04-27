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
# get summary statistics for commobidity conditions for each cohort
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
  
  
  # ── Plot 3: Calories out by cohort ────────────────────────────────────────────
  p_v3_data <- df |>
    select(person_id, cohort, calories_out__mean) |>
    rename(mean_val = calories_out__mean) |>
    mutate(cohort = factor(cohort, levels = c("Insomnia", "Diabetes", "Healthy", "Prostate"))) # healthy
  
  
  kw_v3 <- p_v3_data |>
    kruskal_test(mean_val ~ cohort)
  
  dunn_v3 <- p_v3_data |>
    dunn_test(mean_val ~ cohort, p.adjust.method = "bonferroni") |>
    add_xy_position(x = "cohort")
  
  p_v3 <- ggplot(p_v3_data, aes(x = cohort, y = mean_val, fill = cohort)) +
    geom_violin(color = NA, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.size = 0.3, outlier.alpha = 0.6,
                 fill = "white", color = "gray30") +
    stat_pvalue_manual(dunn_v3, label = "p.adj.signif", tip.length = 0.01,
                       hide.ns = TRUE, size = 3) +
    scale_fill_manual(values = cohort_colors) +
    labs(title    = "Patient-Level Mean Daily Calories Out by Cohort",
         x = NULL, y = "Mean calories/day", fill = "Cohort") +
    theme_bw(base_size = 11) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "bottom")
  ggsave(p_v3, filename = paste0(path.out, "activity_caloriesout_violin.png"), width = 7, height = 6, dpi = 150)
  
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
  
  ggsave(p_mean_sd, filename = paste0(path.out, "activity_mean_sd.png"),            width = 10, height = 7, dpi = 150)
  
  return()
  }





