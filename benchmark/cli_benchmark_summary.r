library(dplyr)
library(gt)
library(readr)

# Read the data
df <- read_csv("benchmark_results.csv")

# Columns to summarize with mean and sd
cols_to_summarize <- c("process_time", "predict_time", "process_max_mem", "predict_max_mem", "process_avg_cpu", "predict_avg_cpu")

# Summarise and combine mean ± sd into one column
summary_df <- df %>%
  group_by(Experiment, Enzyme) %>%
  summarise(
    across(
      all_of(cols_to_summarize),
      list(mean = ~mean(.), sd = ~sd(.)),
      .names = "{.col}_{.fn}"
    ),
    Control_reads = as.integer(mean(Control_reads)),
    Exp_reads = as.integer(mean(Exp_reads))
  ) %>%
  ungroup() %>%
  mutate(
    process_time = sprintf("%.2f ± %.2f", process_time_mean, process_time_sd),
    predict_time = sprintf("%.2f ± %.2f", predict_time_mean, predict_time_sd),
    process_max_mem = sprintf("%.2f ± %.2f", process_max_mem_mean, process_max_mem_sd),
    predict_max_mem = sprintf("%.2f ± %.2f", predict_max_mem_mean, predict_max_mem_sd),
    process_avg_cpu = sprintf("%.2f ± %.2f", process_avg_cpu_mean, process_avg_cpu_sd),
    predict_avg_cpu = sprintf("%.2f ± %.2f", predict_avg_cpu_mean, predict_avg_cpu_sd)
  ) %>%
  select(
    Experiment, Enzyme,
    process_time, predict_time, process_max_mem, predict_max_mem,
    process_avg_cpu, predict_avg_cpu,
    Control_reads, Exp_reads
  )

# Create the gt table
gt_tbl <- summary_df %>%
  gt() %>%
  tab_header(
    title = "Command Line Benchmark Results"
  ) %>%
  cols_label(
    Experiment = "Experiment",
    Enzyme = "Enzyme",
    process_time = "Process Time (sec.)",
    predict_time = "Predict Time (sec.)",
    process_max_mem = "Process Max Mem (MB)",
    predict_max_mem = "Predict Max Mem (MB)",
    process_avg_cpu = "Process Avg CPU (%)",
    predict_avg_cpu = "Predict Avg CPU (%)",
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

# Save the table as a PDF
gtsave(gt_tbl, "cli_benchmark.pdf",  zoom = 0.7)
gtsave(gt_tbl, "cli_benchmark2.html") #, inline_css=TRUE)