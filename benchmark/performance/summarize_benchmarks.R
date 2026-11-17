#!/usr/bin/env Rscript
# Publication-style figures and tables for the TamiPami benchmarks, in R.
#
# Reads every results/cli_*.csv and results/app_*.csv (labels come from git, nothing is hard coded)
# and writes to results/R/:
#   cli_benchmark.{png,pdf}  app_benchmark.{png,pdf}  tables.html
#
# Needs: ggplot2, dplyr, tidyr, readr, patchwork, kableExtra, jsonlite, scales
# (all in the tilapiastats_env conda env). Run from anywhere:
#   ~/miniforge3/envs/tilapiastats_env/bin/Rscript summarize_benchmarks.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(kableExtra); library(jsonlite); library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
res  <- file.path(here, "results")
out  <- file.path(res, "R")
dir.create(out, showWarnings = FALSE)

# ---- load -------------------------------------------------------------------------------------
read_set <- function(prefix) {
  files <- list.files(res, pattern = paste0("^", prefix, "_.*\\.csv$"), full.names = TRUE)
  bind_rows(lapply(files, function(f) {
    meta <- fromJSON(sub("\\.csv$", ".meta.json", f))
    read_csv(f, show_col_types = FALSE) %>%
      mutate(label = meta$label, dirty = meta$dirty, commit = substr(meta$commit, 1, 9),
             stamp = as.Date(meta$timestamp_utc), meta_ok = !is.null(meta$params$cpu_delta))
  }))
}
cli <- read_set("cli") %>% rename_with(tolower)
app <- read_set("app")
stopifnot(nrow(cli) > 0)
sysinfo <- fromJSON(list.files(res, pattern = "\\.meta\\.json$", full.names = TRUE)[1])$system

lab <- function(d, exp = "experiment", enz = "enzyme") {
  d$exp_lab <- factor(paste0(d[[exp]], ": ", d[[enz]]))
  d$exp_lab <- reorder(d$exp_lab, as.numeric(d[[exp]]))
  d
}
cli <- lab(cli); app <- lab(app)
version <- unique(cli$label)[1]
tag <- paste0(version, if (any(cli$dirty)) " (uncommitted changes)" else "", " | n = ",
              max(cli$run), " runs per experiment | ", sysinfo$Chip, ", ", sysinfo$Memory)

# ---- plotting ---------------------------------------------------------------------------------
BLUE <- "#4C78A8"; ORANGE <- "#F58518"
theme_bench <- theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.x = element_text(size = 9, angle = 40, hjust = 1), plot.title = element_text(face = "bold", size = 11),
        axis.title.x = element_blank())

# Bar = mean, error bar = SD, dots = individual runs.
panel <- function(d, y, title, ylab, fill = BLUE) {
  s <- d %>% group_by(exp_lab) %>% summarise(m = mean(.data[[y]]), sd = sd(.data[[y]]), .groups = "drop")
  ggplot(s, aes(exp_lab, m)) +
    geom_col(fill = fill, width = 0.65, alpha = 0.9) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.18, linewidth = 0.4) +
    geom_jitter(data = d, aes(y = .data[[y]]), width = 0.12, size = 0.9, colour = "grey25", alpha = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = label_comma()) +
    labs(title = title, y = ylab) + theme_bench
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h)
}

p_cli <- (panel(cli, "process_time", "Process: wall time", "Seconds") +
          panel(cli, "predict_time", "Predict: wall time", "Seconds", ORANGE)) /
         (panel(cli, "process_max_mem", "Process: peak memory", "MB") +
          panel(cli, "predict_max_mem", "Predict: peak memory", "MB", ORANGE)) +
  plot_annotation(title = "TamiPami command line benchmark", subtitle = tag,
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.subtitle = element_text(colour = "grey35", size = 9)))
save_fig(p_cli, "cli_benchmark", 9, 7)

has_app <- nrow(app) > 0
if (has_app) {
  app_ok_cpu <- all(app$meta_ok)
  panels <- list(panel(app, "elapsed_time_sec", "End-to-end time (submit to results)", "Seconds"),
                 panel(app, "server_max_mem_mb", "Server peak memory", "MB", ORANGE),
                 panel(app, "browser_max_mem_mb", "Browser peak memory", "MB", "#54A24B"))
  if (app_ok_cpu) panels[[4]] <- panel(app, "server_avg_cpu_percent", "Server average CPU", "% of one core", "#B279A2")
  p_app <- wrap_plots(panels, ncol = 2) +
    plot_annotation(title = "TamiPami web application benchmark (Selenium + Chrome)",
                    subtitle = paste0(unique(app$label)[1], " | n = ", max(app$run),
                                      " runs per experiment, fresh server per run"),
                    theme = theme(plot.title = element_text(face = "bold", size = 14),
                                  plot.subtitle = element_text(colour = "grey35", size = 9)))
  save_fig(p_app, "app_benchmark", 9, if (app_ok_cpu) 7 else 7)
}

# ---- tables -----------------------------------------------------------------------------------
msd <- function(x, d = 1) sprintf(paste0("%.", d, "f ± %.", d, "f"), mean(x), sd(x))
comma <- function(x) format(round(x), big.mark = ",")

cli_tab <- cli %>% group_by(Experiment = experiment, Enzyme = enzyme) %>%
  summarise(`Process time (s)` = msd(process_time), `Predict time (s)` = msd(predict_time),
            `Process peak mem (MB)` = msd(process_max_mem, 0), `Predict peak mem (MB)` = msd(predict_max_mem, 0),
            `Process CPU (%)` = msd(process_avg_cpu, 0), `Predict CPU (%)` = msd(predict_avg_cpu, 0),
            `Control read pairs` = comma(mean(control_reads)), `Exp. read pairs` = comma(mean(exp_reads)),
            .groups = "drop") %>% arrange(Experiment)

style_tab <- function(t, caption, note) {
  kbl(t, caption = caption, align = c("c", "l", rep("r", ncol(t) - 2)), escape = FALSE) %>%
    kable_classic(full_width = FALSE, html_font = "Helvetica") %>%
    row_spec(0, bold = TRUE, extra_css = "border-bottom: 2px solid #333;") %>%
    column_spec(1:2, bold = TRUE) %>%
    footnote(general = note, general_title = "", footnote_as_chunk = TRUE)
}

html <- c("<html><head><meta charset='utf-8'><title>TamiPami benchmark tables</title></head><body style='margin:2em'>",
          as.character(style_tab(cli_tab, "Table 1. Command line benchmark (mean ± SD).",
            paste0(version, "; ", tag, ". CPU is total CPU time / wall time (100% = one core). ",
                   "Memory is peak resident set size summed over the process tree, including bbmerge.sh (Java)."))))

if (has_app) {
  app_tab <- app %>% group_by(Experiment = experiment, Enzyme = enzyme) %>%
    summarise(`Time (s)` = msd(elapsed_time_sec),
              `Server baseline (MB)` = msd(server_baseline_mem_mb, 0),
              `Server peak (MB)` = msd(server_max_mem_mb, 0),
              `Server CPU (%)` = if (app_ok_cpu) msd(server_avg_cpu_percent, 0) else NA_character_,
              `Browser peak (MB)` = msd(browser_max_mem_mb, 0), .groups = "drop") %>% arrange(Experiment)
  if (!app_ok_cpu) app_tab$`Server CPU (%)` <- NULL
  html <- c(html, "<br>", as.character(style_tab(app_tab, "Table 2. Web application benchmark (mean ± SD).",
            "Each run used a fresh Streamlit server; baseline is server memory before submitting. Time is submit click to results shown.")))
}

sys_tab <- data.frame(Attribute = c("Model", "Chip", "Cores", "Memory", "macOS", "Browser", "Python"),
                      Value = c(sysinfo$`Model Name`, sysinfo$Chip, sysinfo$`Total Number of Cores`, sysinfo$Memory,
                                sysinfo$macos, sysinfo$chrome, sysinfo$python))
html <- c(html, "<br>", as.character(style_tab(sys_tab, "Table 3. Test system.", "Collected automatically when the benchmarks ran.")),
          "</body></html>")
writeLines(html, file.path(out, "tables.html"))
cat("Wrote", paste(list.files(out), collapse = ", "), "to", out, "\n")
