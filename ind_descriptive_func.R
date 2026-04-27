library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(patchwork)   
#install.packages("corrplot")
library(scales)
library(corrplot)


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# condition_func 
# get summary statistics for commobidity conditions for each cohort
# No condition, Hyperlipidemia 55822004, obesity 414916001, Essential hypertension59621000, coronary arteriosclerosis 53741008
# ───────────────────────────────────────────────────────────────────────────────────────────────────
condition_func <- function(df) {
    n <- length(unique(df$person_id));
    print(n)
    # non of teh conditions
    print('non of the conditions')
    print(nrow(df[which(df$standard_concept_code %in% c('55822004', '414916001', '59621000', '53741008', '19829001')),])) 
    atleast1_df <- df[which(df$standard_concept_code %in% c('55822004', '414916001', '59621000', '53741008', '19829001')),]
    print('non of the conditions')
    print(length(unique(df$person_id)) - length(unique(atleast1_df$person_id)))
    print((length(unique(df$person_id)) - length(unique(atleast1_df$person_id)))/n)
    
    # Hyperlipidemia 55822004
    print('Hyperlipidemia')
    print(nrow(df[which(df$standard_concept_code == '55822004'),])) 
    hyperlip_df <- df[which(df$standard_concept_code == '55822004'),]
    print(length(unique(hyperlip_df$person_id)))
    print(length(unique(hyperlip_df$person_id))/n)
    print(hyperlip_df[1:3,])
          
    # obesity 414916001
    print('obesity')
    print(nrow(df[which(df$standard_concept_code == '414916001'),]) )
    obese_df <- df[which(df$standard_concept_code == '414916001'),]
    print(length(unique(obese_df$person_id)))
    print(length(unique(obese_df$person_id))/n)
    print(obese_df[1:3,])

    # # Essential hypertension59621000
    print('hypertension')
    print(nrow(df[which(df$standard_concept_code == '59621000'),]))
    hypertension_df <-df[which(df$standard_concept_code == '59621000'),]
    print(length(unique(hypertension_df$person_id)))
    print(length(unique(hypertension_df$person_id))/n)
    print(hypertension_df[1:3,])
    
    # pulmonary disease is defined within the structural hierarchy, commonly under the concept Disorder of lung (19829001). 
    print('pulmonary')
    print(nrow(df[which(df$standard_concept_code == '19829001'),]))
    pulmonary_df <- df[which(df$standard_concept_code == '19829001'),]
    print(length(unique(pulmonary_df$person_id)))
    print(length(unique(pulmonary_df$person_id))/n)
    print(pulmonary_df[1:3,])

    #coronary arteriosclerosis 53741008
    print('coronary arteriosclerosis')
    print(nrow(df[which(df$standard_concept_code == '53741008'),]))
    artery_df <- df[which(df$standard_concept_code == '53741008 '),]
    print(length(unique(artery_df$person_id)))
    print(length(unique(artery_df$person_id))/n)
    print(artery_df[1:3,])

    return();
}


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# wear_compliance_func
# get summary statistics for wear compliance 
# barplots for length of FU & total days' of wearing, scatterplot; percentage of wear; hetamap grid of wear
# ───────────────────────────────────────────────────────────────────────────────────────────────────
wear_compliance_func <- function(path.out, df) {
    wear_summary <- df |>
      group_by(person_id) |>
      summarise(
        n_days       = n(),
        first_date   = min(date),
        last_date    = max(date),
        followup_days = as.numeric(last_date - first_date) + 1,
        .groups = "drop"
      )

    cat("── Wear compliance ──\n")
    print(summary(wear_summary[, c("n_days", "followup_days")]))
    
    ##### histogram ####################################################  
    p_wear_ndays <- ggplot(wear_summary, aes(x = n_days)) +
      geom_histogram(bins = 50, fill = "#4A90D9", color = "white") +
      labs(title = "Observations per Patient",
           x = "Number of recorded days", y = "Count of patients") +
      theme_bw()
    ggsave(p_wear_ndays, filename = paste0(path.out, 'activity_nday_dist_histogram.png'), width = 5, height = 3);

    p_wear_followup <- ggplot(wear_summary, aes(x = followup_days)) +
      geom_histogram(bins = 50, fill = "#4A90D9", color = "white") +
      labs(title = "Observations per Patient",
           x = "Follow-up duration (days)", y = "Count of patients") +
      theme_bw()
    ggsave(p_wear_followup, filename = paste0(path.out, 'activity_followup_length_dist_histogram.png'), width = 5, height = 3);

    scatter_wear_compliance <- ggplot(wear_summary, aes(x = followup_days, y = n_days)) +
      geom_point(color = "#4A90D9", size = 1, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "#E67E22", linewidth = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        title = "Recorded Days vs. Follow-up Duration",
        x = "Follow-up duration (days)",
        y = "Number of recorded days"
      ) +
      theme_bw()
    ggsave(scatter_wear_compliance, filename = paste0(path.out, 'activity_followup_vs_nday_scatterplot.png'), width = 5, height = 5);

    wear_summary_percent <- df |>
      group_by(person_id) |>
      summarise(
        n_days        = n(),
        first_date    = min(date),
        last_date     = max(date),
        followup_days = as.numeric(last_date - first_date) + 1,
        .groups       = "drop"
      ) |>
      mutate(
        missing_days    = followup_days - n_days,
        pct_missing     = (missing_days / followup_days) * 100,
        compliance_tier = cut(pct_missing,
                              breaks = c(-Inf, 10, 20, 30, 40, 50, Inf),
                              labels = c("High (>90%)", "80-90%",
                                         "70-80%", "60-70%", "50-60%", "Poor (<50%)"))
      )
    print(summary(wear_summary_percent))
    
    ### Marginal Plot: % Missing distribution ####################################################
    p_main <- ggplot(wear_summary_percent, aes(x = followup_days, y = pct_missing)) +
      geom_point(size = 1, alpha = 0.75) +
      labs(x = "Follow-up duration (days)",
           y = "% Days missing") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none")

    p_top <- ggplot(wear_summary_percent, aes(x = followup_days)) +
      geom_histogram(fill = "#4A90D9", color = "white", bins = 50) +
      labs(y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())

    p_right <- ggplot(wear_summary_percent, aes(x = pct_missing)) +
      geom_histogram(fill = "#8E44AD", color = "white", bins = 50) +
      coord_flip() +
      labs(y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())

    p_marplot <- (p_top + plot_spacer()) / (p_main + p_right) +
      plot_layout(heights = c(1, 6), widths = c(10, 1)) +
      plot_annotation(
        title    = "Marginal Plot: % Missing Days vs. Follow-up Duration",
        theme    = theme(plot.title = element_text(face = "bold"))
      )
    ggsave(p_marplot, filename = paste0(path.out, 'activity_missing_percent_barplot.png'), width = 7, height = 5);

    ##### heatmap - grid ####################################################
    wear_grid <- df |>
      group_by(person_id) |>
      mutate(study_day = as.numeric(date - min(date)) + 1) |>
      ungroup() |>
      select(person_id, study_day) |>
      mutate(has_data = 1)

    # Full grid: every person x every day up to their max study day
    full_grid <- wear_summary |>
      select(person_id, followup_days) |>
      rowwise() |>
      mutate(study_day = list(1:followup_days)) |>
      unnest(study_day) |>
      ungroup()

    # Join to flag observed vs missing
    heatmap_df <- full_grid |>
      left_join(wear_grid, by = c("person_id", "study_day")) |>
      mutate(has_data = replace_na(has_data, 0))

    # Order patients by follow-up length (longest on top)
    patient_order <- wear_summary |>
      arrange(desc(followup_days)) |>
      pull(person_id)

    heatmap_df <- heatmap_df |>
      mutate(person_id = factor(person_id, levels = patient_order))

    heatmap <- ggplot(heatmap_df, aes(x = study_day, y = person_id, fill = factor(has_data))) +
      geom_tile(linewidth = 0) +
      scale_fill_manual(values = c("0" = "white", "1" = "black"),
                        labels = c("0" = "Missing", "1" = "Observed"),
                        name   = NULL) +
      labs(
        title    = "Patient Wear Compliance Heatmap",
        x        = "Study day",
        y        = "Patient ID"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid       = element_blank(),
        axis.text.y = element_blank(),
        legend.position  = "bottom",
        legend.background = element_rect(color = "black", linewidth = 0.5)
      )
    ggsave(heatmap, filename = paste0(path.out, 'activity_missing_heatmap.png'), width = 7, height = 5);
    
    return()
}


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# activity_prep_stats
# get cohort level stats for each vairable
# output the dataset into path, and output summary table
# ───────────────────────────────────────────────────────────────────────────────────────────────────
activity_prep_stats <- function(path.out, df, cohort) {
    # ── 1. Compute patient-level stats ─────────────────────────────────────────
    key_metrics <- c("steps", "calories_out", "sedentary_minutes",
                    "lightly_active_minutes", "fairly_active_minutes",
                    "very_active_minutes", "sex_at_birth")

    patient_stats <- df |>
    group_by(person_id) |>
    summarise(
        n_obs = n(),
        across(all_of(key_metrics), list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x,   na.rm = TRUE),
        cv   = ~ifelse(mean(.x, na.rm = TRUE) > 0,
                        sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE) * 100, NA)
        ), .names = "{.col}__{.fn}"),
        .groups = "drop"
    )
    patient_stats$cohort <- cohort
    write.table(patient_stats, paste0(path.out, cohort,'_activity_summary_stats.txt'), sep = '\t', quote = F, row.names = F)

    # ── 2. Summary table: distribution of patient means ────────────────────────
    mean_cols <- paste0(key_metrics, "__mean")

    cohort_summary <- patient_stats |>
    select(all_of(mean_cols)) |>
    pivot_longer(everything(), names_to = "metric", values_to = "value") |>
    mutate(metric = str_remove(metric, "__mean")) |>
    group_by(metric) |>
    summarise(
        N        = sum(!is.na(value)),
        Mean     = mean(value,              na.rm = TRUE),
        SD       = sd(value,                na.rm = TRUE),
        Median   = median(value,            na.rm = TRUE),
        Q1       = quantile(value, 0.25,    na.rm = TRUE),
        Q3       = quantile(value, 0.75,    na.rm = TRUE),
        P5       = quantile(value, 0.05,    na.rm = TRUE),
        P95      = quantile(value, 0.95,    na.rm = TRUE),
        .groups  = "drop"
    )
    print(cohort_summary, n = Inf)
    return(patient_stats);
}


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# miss_gap_func
# get the missing cohort statistics and curve
# ───────────────────────────────────────────────────────────────────────────────────────────────────
miss_gap_func <- function(path.out, df) {
  gap_summary <- df |>
    group_by(person_id) |>
    summarise(
      first_date = min(date),
      last_date  = max(date),
      obs_dates  = list(sort(unique(date))),
      .groups    = "drop"
    ) |>
    rowwise() |>
    mutate(
      # Full expected sequence from first to last date
      full_seq      = list(seq(first_date, last_date, by = "day")),
      # Missing days = expected minus observed
      missing_seq   = list(setdiff(as.numeric(full_seq),
                                   as.numeric(obs_dates))),
      total_missing = length(missing_seq),
      followup_days = as.numeric(last_date - first_date) + 1,
      
      # Identify consecutive runs in missing days
      gap_lengths = list({
        m <- missing_seq
        if (length(m) == 0) {
          integer(0)
        } else {
          # Label each missing day with a run ID
          run_id  <- cumsum(c(1L, diff(m) != 1L))
          as.integer(table(run_id))
        }
      }),
      
      n_gaps   = length(gap_lengths),
      max_gap  = ifelse(n_gaps == 0, 0L, max(gap_lengths)),
      mean_gap = ifelse(n_gaps == 0, 0,  mean(gap_lengths))
    ) |>
    ungroup() |>
    select(person_id, first_date, last_date,
           followup_days, total_missing,
           n_gaps, max_gap, mean_gap, gap_lengths)
  print(summary(gap_summary$max_gap))
  # ── 3. Threshold sensitivity table ───────────────────────────────────────────
  # How many patients would be excluded at each candidate threshold?
  thresholds <- c(0, 7, 14, 30, 60)
  
  threshold_tbl <- map_dfr(thresholds, ~tibble(
    threshold      = .x,
    n_excluded     = sum(gap_summary$max_gap > .x, na.rm = TRUE),
    pct_excluded   = mean(gap_summary$max_gap > .x, na.rm = TRUE) * 100,
    n_retained     = sum(gap_summary$max_gap <= .x, na.rm = TRUE),
    pct_retained   = mean(gap_summary$max_gap <= .x, na.rm = TRUE) * 100
  ))
  
  cat("\n── Exclusion sensitivity by threshold ──\n")
  print(threshold_tbl, n = Inf)
  p_ecdf <- gap_summary |>
    mutate(max_gap = pmin(max_gap, 60)) |>
    ggplot(aes(x = max_gap)) +
    stat_ecdf(color = "#2C3E50", linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 60),
                       breaks = c(0, 14, 30, 60)) +
    #scale_y_continuous(labels = scales::percent) +
    labs(title    = "ECDF of Max Continuous Missing Days",
         subtitle = "Y = proportion of patients with max gap ≤ x days",
         x = "Max consecutive missing days",
         y = "Cumulative proportion retained") +
    theme_bw(base_size = 11)

  ggsave(p_ecdf, filename = paste0(path.out, 'activity_missing_gap.png'), width = 9, height = 5);
  return();
}

