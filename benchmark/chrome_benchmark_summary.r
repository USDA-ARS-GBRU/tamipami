library(dplyr)
library(gt)
library(readr)

# Read chrome_app_results.csv
chrome_app <- read_csv("chrome_app_results.csv")

# Summarize by experiment and enzyme
summary_df <- chrome_app %>%
  group_by(experiment, enzyme) %>%
  summarise(
    elapsed_time_sec_mean = mean(elapsed_time_sec, na.rm = TRUE),
    elapsed_time_sec_sd = sd(elapsed_time_sec, na.rm = TRUE),
    max_mem_mb_mean = mean(max_mem_mb, na.rm = TRUE),
    max_mem_mb_sd = sd(max_mem_mb, na.rm = TRUE),
    avg_cpu_percent_mean = mean(avg_cpu_percent, na.rm = TRUE),
    avg_cpu_percent_sd = sd(avg_cpu_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    elapsed_time_sec = sprintf("%.2f ± %.2f", elapsed_time_sec_mean, elapsed_time_sec_sd),
    max_mem_mb = sprintf("%.2f ± %.2f", max_mem_mb_mean, max_mem_mb_sd),
    avg_cpu_percent = sprintf("%.4f ± %.4f", avg_cpu_percent_mean, avg_cpu_percent_sd)
  )

# Read benchmark_results.csv for read counts
benchmark_results <- read_csv("benchmark_results.csv")

# Summarize read counts by Experiment
reads_summary <- benchmark_results %>%
  group_by(Experiment) %>%
  summarise(
    Control_reads = as.integer(mean(Control_reads, na.rm = TRUE)),
    Exp_reads = as.integer(mean(Exp_reads, na.rm = TRUE)),
    .groups = "drop"
  )

# Merge summaries
final_summary <- summary_df %>%
  left_join(reads_summary, by = c("experiment" = "Experiment")) %>%
  select(
    experiment, enzyme,
    elapsed_time_sec, max_mem_mb, avg_cpu_percent,
    Control_reads, Exp_reads
  )

# Create the gt table
gt_tbl <- final_summary %>%
  gt() %>%
  tab_header(
    title = "Chrome App Benchmark Results"
  ) %>%
  cols_label(
    experiment = "Experiment",
    enzyme = "Enzyme",
    elapsed_time_sec = "Elapsed Time (sec.)",
    max_mem_mb = "Max Mem (MB)",
    avg_cpu_percent = "Avg CPU (%)",
    Control_reads = "Control Read pairs",
    Exp_reads = "Experimental Read pairs"
  ) %>%
  fmt_number(
    columns = c("Control_reads", "Exp_reads"),
    decimals = 0,
    use_seps = TRUE
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 18,
    heading.title.font.weight = "bold"
  )

# Print the table
print(gt_tbl)

# Save the table as a PDF and HTML
gtsave(gt_tbl, "chrome_app_benchmark.pdf", zoom = 0.7)
gtsave(gt_tbl, "chrome_app_benchmark.html", inline_css=TRUE)