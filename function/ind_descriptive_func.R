library(tidyverse)
library(readxl)
library(lubridate)
library(ggplot2)
library(patchwork)   
#install.packages("corrplot")
library(scales)
library(corrplot)
library(data.table)
library(stringr)

# ── Helper: format as "n (x%)" ───────────────────────────────────────────────
fmt <- function(n, total) {
  pct <- round(n / total * 100)
  if (pct < 1) paste0(format(n, big.mark = ","), " (<1%)")
  else         paste0(format(n, big.mark = ","), " (", pct, "%)")
}

# ───────────────────────────────────────────────────────────────────────────────────────────────────
# condition_func 
# get summary statistics for commobidity conditions for each cohort
# No condition, Hyperlipidemia 55822004, obesity 414916001, Essential hypertension59621000, 
# coronary arteriosclerosis 53741008 diabetes mellitus, 73211009
# ───────────────────────────────────────────────────────────────────────────────────────────────────
condition_func <- function(df) {
    n <- length(unique(df$person_id));
    print(n)
    # non of teh conditions
    print('non of the conditions')
    print(nrow(df[which(df$standard_concept_code %in% c('55822004', '414916001', '59621000', '53741008', '19829001', '73211009')),])) 
    atleast1_df <- df[which(df$standard_concept_code %in% c('55822004', '414916001', '59621000', '53741008', '19829001', '73211009')),]
    print('non of the conditions')
    print(length(unique(df$person_id)) - length(unique(atleast1_df$person_id)))
    print((length(unique(df$person_id)) - length(unique(atleast1_df$person_id)))/n)
    
    # Hyperlipidemia 55822004
    print('Hyperlipidemia')
    print(nrow(df[which(df$standard_concept_code == '55822004'),])) 
    hyperlip_df <- df[which(df$standard_concept_code == '55822004'),]
    print(length(unique(hyperlip_df$person_id)))
    print(length(unique(hyperlip_df$person_id))/n)
    print(hyperlip_df[1:3,1:3])
          
    # obesity 414916001
    print('obesity')
    print(nrow(df[which(df$standard_concept_code == '414916001'),]) )
    obese_df <- df[which(df$standard_concept_code == '414916001'),]
    print(length(unique(obese_df$person_id)))
    print(length(unique(obese_df$person_id))/n)
    print(obese_df[1:3,1:3])

    # # Essential hypertension59621000
    print('hypertension')
    print(nrow(df[which(df$standard_concept_code == '59621000'),]))
    hypertension_df <-df[which(df$standard_concept_code == '59621000'),]
    print(length(unique(hypertension_df$person_id)))
    print(length(unique(hypertension_df$person_id))/n)
    print(hypertension_df[1:3,1:3])
    
    # pulmonary disease is defined within the structural hierarchy, commonly under the concept Disorder of lung (19829001). 
    print('pulmonary')
    print(nrow(df[which(df$standard_concept_code == '19829001'),]))
    pulmonary_df <- df[which(df$standard_concept_code == '19829001'),]
    print(length(unique(pulmonary_df$person_id)))
    print(length(unique(pulmonary_df$person_id))/n)
    print(pulmonary_df[1:3,1:3])

    #coronary arteriosclerosis 53741008
    print('coronary arteriosclerosis')
    print(nrow(df[which(df$standard_concept_code == '53741008'),]))
    artery_df <- df[which(df$standard_concept_code == '53741008'),]
    print(length(unique(artery_df$person_id)))
    print(length(unique(artery_df$person_id))/n)
    print(artery_df[1:3,])
    
    #diabetes mellitus, 73211009
    print('diabetes mellitus')
    print(nrow(df[which(df$standard_concept_code == '73211009'),]))
    diabetes_df <- df[which(df$standard_concept_code == '73211009'),]
    print(length(unique(diabetes_df$person_id)))
    print(length(unique(diabetes_df$person_id))/n)
    print(diabetes_df[1:3,1:3])
    
    return();
}




# ───────────────────────────────────────────────────────────────────────────────────────────────────
# demographic_tbl_func 
# get summary statistics for demographic meta data
# ───────────────────────────────────────────────────────────────────────────────────────────────────
demographic_tbl_func <- function(df) {
    # ── 1. Derive age from date_of_birth ─────────────────────────────────────────
    age <- df %>%
      mutate(
        dob       = as.Date(substr(date_of_birth, 1, 10)),
        age       = as.integer(floor(interval(dob, Sys.Date()) / years(1))),
        age_group = case_when(
          age >= 18 & age <= 44 ~ "18-44",
          age >= 45 & age <= 64 ~ "45-64",
          age >= 65              ~ ">65",
          TRUE                   ~ NA_character_
        )
      )
    
    N <- nrow(df)
    
    # ── 2. Age ────────────────────────────────────────────────────────────────────
    age_tbl <- age %>%
      count(age_group) %>%
      mutate(
        Demographic = age_group,
        Value       = map_chr(n, ~ fmt(.x, N))
      ) %>%
      filter(!is.na(Demographic)) %>%
      arrange(Demographic) %>%
      select(Demographic, Value)
    
    # ── 3. Ethnicity ──────────────────────────────────────────────────────────────
    # Map ethnicity_concept_id:
    # everything else = Others
    eth_tbl <- df %>%
      mutate(eth_group = case_when(
        ethnicity_concept_id == 38003563 ~ "Hispanic or Latino",
        ethnicity_concept_id == 38003564 ~ "Not Hispanic or Latino",
        TRUE                             ~ "Others"
      )) %>%
      count(eth_group) %>%
      mutate(
        Demographic = eth_group,
        Value       = map_chr(n, ~ fmt(.x, N))
      ) %>%
      select(Demographic, Value)
    
    # ── 4. Race ───────────────────────────────────────────────────────────────────
    # Standard All of Us race concept IDs
    race_map <- c(
      "8657" = "American Indian or Alaska Native",
      "8515" = "Asian",
      "8516" = "Black or African American",
      "38003615" = "Middle Eastern or North African",
      "8657" = "Native Hawaiian or Other Pacific Islander",
      "8527" = "White",
      "2000000008" = "More than one population"
    )
    
    race_tbl <- df %>%
      mutate(race_group = case_when(
        as.character(race_concept_id) %in% names(race_map) ~
          race_map[as.character(race_concept_id)],
        TRUE ~ "Others"
      )) %>%
      count(race_group) %>%
      mutate(
        Demographic = race_group,
        Value       = map_chr(n, ~ fmt(.x, N))
      ) %>%
      arrange(Demographic) %>%
      select(Demographic, Value)
    
    # ── 5. Sex at birth ───────────────────────────────────────────────────────────
    # 45878463 = Female, 45873798 = Male, else Others
    sex_tbl <- df %>%
      mutate(sex_group = case_when(
        sex_at_birth_concept_id == 45878463 ~ "Female",
        sex_at_birth_concept_id == 45880669 ~ "Male",
        TRUE                                 ~ "Others"
      )) %>%
      count(sex_group) %>%
      mutate(
        Demographic = sex_group,
        Value       = map_chr(n, ~ fmt(.x, N))
      ) %>%
      select(Demographic, Value)
    
    # ── 6. Assemble the full table ────────────────────────────────────────────────
    demo_table <- bind_rows(
      tibble(Demographic = "N",              Value = format(N, big.mark = ",")),
      tibble(Demographic = "Age",            Value = ""),
      age_tbl,
      tibble(Demographic = "Ethnicity",      Value = ""),
      eth_tbl,
      tibble(Demographic = "Race",           Value = ""),
      race_tbl,
      tibble(Demographic = "Sex at birth",   Value = ""),
      sex_tbl
    )
    
    print(demo_table, n = Inf)

  return()
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
activity_prep_stats <- function(path.out, df, cohort, meta.data) {
    # ── 1. Compute patient-level stats ─────────────────────────────────────────
    key_metrics <- c("steps", "calories_out", "sedentary_minutes",
                    "lightly_active_minutes", "fairly_active_minutes",
                    "very_active_minutes")

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
    patient_stats <-  merge(patient_stats, meta.data[,c('person_id', 'sex_at_birth')], by = 'person_id')
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
# sleep_sum_main_nap_func
# get cohort level stats for each variable
# output the dataset into path, and output summary table
# ───────────────────────────────────────────────────────────────────────────────────────────────────
sleep_sum_main_nap_func <- function(df) {
  sleep_sum <- df |>
    group_by(person_id, date) |>
    summarise(
      # Flags — did this day have main sleep or nap?
      has_main_sleep = any(is_main_sleep == "true"),
      has_nap = any(is_main_sleep == "false"),
      n_records = n(),
      
      # Sum across main + nap for that day
      minute_in_bed = sum(minute_in_bed,       na.rm = TRUE),  # Total time in bed
      minute_asleep = sum(minute_asleep,        na.rm = TRUE),  # Total sleep time
      minute_awake = sum(minute_awake,         na.rm = TRUE),  # Wake after sleep onset
      minute_restless = sum(minute_restless,      na.rm = TRUE),
      minute_deep = sum(minute_deep,      na.rm = TRUE),
      minute_light = sum(minute_light,      na.rm = TRUE),
      minute_rem = sum(minute_rem,      na.rm = TRUE),
      minute_wake = sum(minute_wake,      na.rm = TRUE),
      minute_after_wakeup = sum(minute_after_wakeup,  na.rm = TRUE),
      main_minute_in_bed = sum(minute_in_bed      * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_asleep = sum(minute_asleep      * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_awake = sum(minute_awake       * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_deep = sum(minute_deep       * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_light = sum(minute_light       * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_rem = sum(minute_rem       * (is_main_sleep == "true"), na.rm = TRUE),
      main_minute_wake = sum(minute_wake       * (is_main_sleep == "true"), na.rm = TRUE),
      
      # Daily sleep efficiency based on totals
      sleep_efficiency  = minute_asleep / minute_in_bed * 100,
      
      .groups = "drop"
  )
  return(sleep_sum)

}

# ───────────────────────────────────────────────────────────────────────────────────────────────────
# sleep_prep_stats
# get cohort level stats for each variable
# output the dataset into path, and output summary table; main = "main"/"nap"/"sum"
# ───────────────────────────────────────────────────────────────────────────────────────────────────
library(data.table)
library(stringr)

sleep_prep_stats <- function(path.out, df, cohort, main) {
  # ── 1. Compute patient-level stats ─────────────────────────────────────────
  if (main == "sum") {
    key_metrics <- c("minute_in_bed", "minute_asleep", "minute_awake",
                     "minute_deep", "minute_light", "minute_rem", "minute_wake", "sleep_efficiency")
  } else {
    key_metrics <- c("minute_in_bed", "minute_asleep", "minute_awake",
                     "minute_deep", "minute_light", "minute_rem", "minute_wake")
    
  }
  
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
  write.table(patient_stats, paste0(path.out, cohort,'_', main, '_sleep_summary_stats.txt'), sep = '\t', quote = F, row.names = F)
  
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


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# filter_pid_func
# output person_id with  longest consecutive run of missing days ≤ 14;≥ 70% compliance
# ───────────────────────────────────────────────────────────────────────────────────────────────────
filter_pid_func <- function(path.out, df, cohort, type) {

  wear_summary <- df |>
    group_by(person_id) |>
    summarise(
      n_days        = n_distinct(date),        # ✅ unique dates
      first_date    = min(date),
      last_date     = max(date),
      followup_days = as.numeric(last_date - first_date) + 1,
      obs_dates     = list(sort(unique(date))),
      .groups       = "drop"
    ) |>
    mutate(
      pct_wear = (n_days / followup_days) * 100
    ) |>
    rowwise() |>
    mutate(
      full_seq    = list(seq(first_date, last_date, by = "day")),
      missing_seq = list(setdiff(as.numeric(full_seq), as.numeric(obs_dates))),
      max_gap     = {
        m <- missing_seq
        if (length(m) == 0) 0L
        else max(as.integer(table(cumsum(c(1L, diff(m) != 1L)))))
      }
    ) |>
    ungroup()
  
  # Inspect before filtering
  cat("Before filter:", nrow(wear_summary), "patients\n")
  cat("pct_wear >= 70:", sum(wear_summary$pct_wear >= 70), "\n")
  cat("max_gap <= 14: ", sum(wear_summary$max_gap <= 14), "\n")
  
  eligible_ids <- wear_summary |>
    filter(
      pct_wear >= 70,
      max_gap  <= 14
    ) |>
    pull(person_id)
  
  cat("After filter:", length(eligible_ids), "patients\n")
  
  write.table(
    eligible_ids,
    paste0(path.out, cohort, "_", type, "_pid_filtered.txt"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  return(invisible(eligible_ids))

}



# ───────────────────────────────────────────────────────────────────────────────────────────────────
# sleep_prep_stats_from_gcs
# read large files for healthy cohort
# ───────────────────────────────────────────────────────────────────────────────────────────────────
sum_main_nap_dt <- function(dt) {
  dt_agg <- dt[, .(
    has_main_sleep      = any(is_main_sleep == "true"),
    has_nap             = any(is_main_sleep == "false"),
    n_records           = .N,
    minute_in_bed       = sum(minute_in_bed,       na.rm = TRUE),
    minute_asleep       = sum(minute_asleep,        na.rm = TRUE),
    minute_awake        = sum(minute_awake,         na.rm = TRUE),
    minute_restless     = sum(minute_restless,      na.rm = TRUE),
    minute_deep         = sum(minute_deep,          na.rm = TRUE),
    minute_light        = sum(minute_light,         na.rm = TRUE),
    minute_rem          = sum(minute_rem,           na.rm = TRUE),
    minute_wake         = sum(minute_wake,          na.rm = TRUE),
    minute_after_wakeup = sum(minute_after_wakeup,  na.rm = TRUE),
    main_minute_in_bed  = sum(minute_in_bed  * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_asleep  = sum(minute_asleep  * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_awake   = sum(minute_awake   * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_deep    = sum(minute_deep    * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_light   = sum(minute_light   * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_rem     = sum(minute_rem     * (is_main_sleep == "true"), na.rm = TRUE),
    main_minute_wake    = sum(minute_wake    * (is_main_sleep == "true"), na.rm = TRUE)
    # sleep_efficiency intentionally excluded — computed after full dedup
  ), by = .(person_id, sleep_date)]
}

sleep_pipeline_from_gcs <- function(path.out, export_path, cohort, main) {
  
  if (main == "sum") {
    # sleep_efficiency excluded from key_metrics — computed separately below
    key_metrics <- c("minute_in_bed", "minute_asleep", "minute_awake",
                     "minute_deep", "minute_light", "minute_rem", "minute_wake")
    include_efficiency <- TRUE
  } else {
    key_metrics        <- c("minute_in_bed", "minute_asleep", "minute_awake",
                            "minute_deep", "minute_light", "minute_rem", "minute_wake")
    include_efficiency <- FALSE
  }
  
  raw_cols_needed <- c("person_id", "sleep_date", "is_main_sleep",
                       "minute_in_bed", "minute_asleep", "minute_awake",
                       "minute_restless", "minute_deep", "minute_light",
                       "minute_rem", "minute_wake", "minute_after_wakeup")
  
  csv_files <- system2('gsutil', args = c('ls', export_path),
                       stdout = TRUE, stderr = FALSE)
  cat("Found", length(csv_files), "CSV files.\n")
  dir.create(path.out, showWarnings = FALSE, recursive = TRUE)
  
  daily_file   <- paste0(path.out, cohort, "_", main, "_daily_summary.txt")
  stats_file   <- paste0(path.out, cohort, "_", main, "_sleep_summary_stats.txt")
  log_file     <- paste0(path.out, cohort, "_", main, "_processed_files.log")
  first_write  <- TRUE
  
  # Track already-processed files to allow safe rerun
  done_files <- if (file.exists(log_file)) readLines(log_file) else character(0)
  
  for (i in seq_along(csv_files)) {
    csv      <- csv_files[[i]]
    csv_base <- basename(csv)
    if (csv_base %in% done_files) {
      cat("Skipping (already processed):", csv_base, "\n")
      first_write <- FALSE   # daily_file already has a header
      next
    }
    
    cat("File", i, "/", length(csv_files), ":", csv_base, "\n")
    
    tryCatch({
      # Safely select only columns that actually exist in this shard
      all_cols  <- names(fread(cmd = paste("gsutil cat", shQuote(csv)),
                               nrows = 0))
      use_cols  <- intersect(raw_cols_needed, all_cols)
      missing   <- setdiff(raw_cols_needed, all_cols)
      if (length(missing) > 0)
        cat("  WARNING: missing columns:", paste(missing, collapse = ", "), "\n")
      
      dt_raw   <- fread(cmd = paste("gsutil cat", shQuote(csv)),
                        select = use_cols, na.strings = c("", "NA", "NULL"))
      
      # Add missing columns as NA so downstream code doesn't break
      for (col in missing) dt_raw[, (col) := NA_real_]
      
      dt_daily <- sum_main_nap_dt(dt_raw)
      rm(dt_raw); gc()
      
      fwrite(dt_daily, daily_file, sep = "\t", quote = FALSE,
             append = !first_write, col.names = first_write)
      first_write <- FALSE
      rm(dt_daily); gc()
      
      # Log success — enables safe rerun
      cat(csv_base, file = log_file, append = TRUE, sep = "\n")
      
    }, error = function(e) {
      cat("  ERROR on", csv_base, ":", conditionMessage(e), "\n")
    })
  }
  
  # ── Dedup: collapse any (person_id, sleep_date) split across shards ─────────
  cat("Deduplicating daily summary...\n")
  dt_daily_all <- fread(daily_file)
  
  minute_cols <- c("minute_in_bed", "minute_asleep", "minute_awake",
                   "minute_restless", "minute_deep", "minute_light",
                   "minute_rem", "minute_wake", "minute_after_wakeup",
                   "main_minute_in_bed", "main_minute_asleep", "main_minute_awake",
                   "main_minute_deep", "main_minute_light", "main_minute_rem",
                   "main_minute_wake")
  
  dt_daily_all <- dt_daily_all[, c(
    list(has_main_sleep = any(has_main_sleep),
         has_nap        = any(has_nap),
         n_records      = sum(n_records)),
    lapply(.SD, sum, na.rm = TRUE)
  ), by = .(person_id, sleep_date), .SDcols = minute_cols]
  
  # Compute sleep_efficiency AFTER dedup, from correct totals
  dt_daily_all[, sleep_efficiency := minute_asleep / minute_in_bed * 100]
  
  # Check for implausible values
  n_bad_eff <- dt_daily_all[sleep_efficiency > 100 | sleep_efficiency < 0, .N]
  if (n_bad_eff > 0)
    cat("  WARNING:", n_bad_eff, "rows with sleep_efficiency outside [0,100]\n")
  
  fwrite(dt_daily_all, daily_file, sep = "\t", quote = FALSE)
  
  # ── Patient-level stats ───────────────────────────────────────────────────
  cat("Computing patient stats...\n")
  stats_metrics <- if (include_efficiency) c(key_metrics, "sleep_efficiency") else key_metrics
  
  patient_stats <- dt_daily_all[, {
    result <- list(n_days = .N,
                   pct_days_with_nap       = mean(has_nap,        na.rm = TRUE) * 100,
                   pct_days_with_mainsleep = mean(has_main_sleep,  na.rm = TRUE) * 100)
    for (col in stats_metrics) {
      x <- get(col)
      m <- mean(x, na.rm = TRUE)
      s <- sd(x,   na.rm = TRUE)
      result[[paste0(col, "__mean")]] <- m
      result[[paste0(col, "__sd")]]   <- s
      result[[paste0(col, "__cv")]]   <- ifelse(m > 0, s / m * 100, NA_real_)
    }
    # Patient-level efficiency from totals — more accurate than mean of daily ratios
    result[["sleep_efficiency_overall"]] <-
      sum(minute_asleep, na.rm = TRUE) / sum(minute_in_bed, na.rm = TRUE) * 100
    result
  }, by = person_id]
  
  patient_stats[, cohort := cohort]
  fwrite(patient_stats, stats_file, sep = "\t", quote = FALSE)
  rm(dt_daily_all); gc()
  
  # ── Cohort summary ────────────────────────────────────────────────────────
  mean_cols <- paste0(stats_metrics, "__mean")
  cohort_summary <- melt(
    patient_stats[, c("person_id", mean_cols), with = FALSE],
    id.vars = "person_id", variable.name = "metric", value.name = "value"
  )[, metric := str_remove(as.character(metric), "__mean")
  ][, .(N      = sum(!is.na(value)),
        Mean   = mean(value,           na.rm = TRUE),
        SD     = sd(value,             na.rm = TRUE),
        Median = median(value,         na.rm = TRUE),
        Q1     = quantile(value, 0.25, na.rm = TRUE),
        Q3     = quantile(value, 0.75, na.rm = TRUE),
        P5     = quantile(value, 0.05, na.rm = TRUE),
        P95    = quantile(value, 0.95, na.rm = TRUE)),
    by = metric]
  
  print(as.data.frame(cohort_summary))
  return(invisible(patient_stats))
}


# ───────────────────────────────────────────────────────────────────────────────────────────────────
# hr_prep_stats
# get cohort level stats for each variable
# output the dataset into path, and output summary table
# ───────────────────────────────────────────────────────────────────────────────────────────────────
hr_prep_stats <- function(path.out, df, cohort) {
  # ── 1. Compute patient-level stats ─────────────────────────────────────────
  key_metrics <- c("min_heart_rate", "max_heart_rate", "minute_in_zone", "calorie_count")
  
  patient_stats <- df |>
    group_by(person_id, zone_name) |>
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
  write.table(patient_stats, paste0(path.out, cohort,'_hr_summary_stats.txt'), sep = '\t', quote = F, row.names = F)
  
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