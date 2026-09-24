# Summary plots and paper tables for the original (non-JTPA) simulations.
# Produces four PDFs and four LaTeX tables: one complier/overall pair and one
# nudge pair for each of XGBoost and MGCV.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

quiet_source <- function(file, local = .GlobalEnv) {
  invisible(capture.output(
    suppressPackageStartupMessages(suppressMessages(source(file, local = local)))
  ))
}

quiet_source("LIV_xgboost.R")

SD_quantile <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  as.numeric((q[2] - q[1]) / (qnorm(0.75) - qnorm(0.25)))
}

set.seed(2025)
methods <- c("xgb", "mgcv")
sample_sizes <- c(2000L, 5000L)
policy_labels <- c("00", "01", "10", "11")
actions <- list(
  "00" = c(0, 0),
  "01" = c(0, 1),
  "10" = c(1, 0),
  "11" = c(1, 1)
)
dynamic_policy_source <- c("00" = "(A_0,0)", "01" = "(A_0,1)")
intervention_labels <- c(
  "(0,0)", "(0,1)", "(1,0)", "(1,1)", "(A_0,0)", "(A_0,1)"
)
truth_accuracy <- as.integer(Sys.getenv("SUMMARY_TRUTH_ACCURACY", "1000000"))
if (is.na(truth_accuracy) || truth_accuracy < 1L) {
  stop("SUMMARY_TRUTH_ACCURACY must be a positive integer.", call. = FALSE)
}
output_dir <- file.path("outputs", "simulation")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

model_spec <- list(
  overall = list(
    parameter_file = "Script/par/par_overall.R",
    result_prefix = "Overall",
    truth_function = ATE_overall
  ),
  complier = list(
    parameter_file = "Script/par/par_complier.R",
    result_prefix = "Complier",
    truth_function = ATE_local
  ),
  nudge = list(
    parameter_file = "Script/par/par_nudge.R",
    result_prefix = "Nudge",
    truth_function = ATE_local
  )
)

load_result <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("RES", envir = env)) {
    stop("Missing RES in ", path, call. = FALSE)
  }
  env$RES
}

source_model <- function(estimand) {
  spec <- model_spec[[estimand]]
  if (is.null(spec)) stop("Unknown estimand: ", estimand, call. = FALSE)
  quiet_source(spec$parameter_file)
  list(D = D, TM = TM, truth_function = spec$truth_function)
}

make_truth <- function(estimand) {
  model <- source_model(estimand)
  truth <- lapply(actions, function(action) {
    suppressMessages(model$truth_function(
      model$D, model$TM, PROB = action, accuracy = truth_accuracy
    ))
  })
  list(TM = model$TM, truth = truth)
}

summarize_one <- function(method, estimand, n, policy, truth) {
  spec <- model_spec[[estimand]]
  path <- file.path(
    "data", paste0(method, n), paste0(spec$result_prefix, "_", policy, ".RData")
  )
  if (!file.exists(path)) stop("Result file not found: ", path, call. = FALSE)

  RES <- load_result(path)
  times <- seq_along(truth)
  if (length(dim(RES)) != 3L || dim(RES)[1] != length(times) || dim(RES)[2] < 5L) {
    stop(
      "Unexpected RES dimensions in ", path, "; expected ", length(times),
      " x at least 5 x repetitions, got ", paste(dim(RES), collapse = " x "),
      call. = FALSE
    )
  }

  estimate <- apply(RES[, 2, ], 1, median, na.rm = TRUE)
  data.frame(
    method = method,
    estimand = estimand,
    n = n,
    policy = policy,
    Time = times - 1L,
    Truth = truth,
    Estimate = estimate,
    Bias = abs(estimate - truth),
    SD = apply(RES[, 2, ], 1, SD_quantile),
    SE = apply(RES[, 3, ], 1, median, na.rm = TRUE),
    Lower = apply(RES[, 4, ], 1, median, na.rm = TRUE),
    Upper = apply(RES[, 5, ], 1, median, na.rm = TRUE),
    CR95 = vapply(times, function(t) {
      mean(RES[t, 4, ] <= truth[t] & RES[t, 5, ] >= truth[t], na.rm = TRUE)
    }, numeric(1))
  )
}

make_summary <- function(method, estimand, truth_obj) {
  bind_rows(lapply(sample_sizes, function(n) {
    bind_rows(lapply(policy_labels, function(policy) {
      summarize_one(method, estimand, n, policy, truth_obj$truth[[policy]])
    }))
  }))
}

# Time 0 fixes both treatments.  Time 1 leaves A_0 natural, which gives the
# original table's Y(A_0, 0) and Y(A_0, 1) entries.
make_paper_table <- function(summary_df) {
  fixed <- summary_df %>%
    filter(Time == 0) %>%
    mutate(intervention = paste0("(", substr(policy, 1, 1), ",", substr(policy, 2, 2), ")"))
  natural_a0 <- summary_df %>%
    filter(Time == 1, policy %in% names(dynamic_policy_source)) %>%
    mutate(intervention = unname(dynamic_policy_source[policy]))

  bind_rows(fixed, natural_a0) %>%
    select(method, estimand, n, Time, policy, intervention, Bias, SD, SE, CR95)
}

metric_names <- c("Bias", "SD", "SE", "CR95")

extract_metric_values <- function(paper_df, estimand_name, sample_size, metric_name) {
  values <- paper_df %>%
    filter(.data$estimand == estimand_name, n == sample_size) %>%
    arrange(match(intervention, intervention_labels))
  if (nrow(values) != length(intervention_labels) ||
      anyDuplicated(values$intervention)) {
    stop("Expected six unique interventions for ", estimand_name,
         ", n = ", sample_size, call. = FALSE)
  }
  metric <- values[[metric_name]]
  names(metric) <- values$intervention
  metric[intervention_labels]
}

make_integrated_table <- function(paper_df) {
  rows <- list()
  for (sample_size in sample_sizes) {
    for (metric_name in metric_names) {
      values <- c(
        extract_metric_values(paper_df, "complier", sample_size, metric_name),
        extract_metric_values(paper_df, "overall", sample_size, metric_name)
      )
      rows[[paste(sample_size, metric_name, sep = "_")]] <- data.frame(
        Size = sample_size,
        Metric = ifelse(metric_name == "CR95", "CR", metric_name),
        t(values),
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  names(out) <- c(
    "Size", "Metric",
    paste0("Complier_", intervention_labels),
    paste0("Overall_", intervention_labels)
  )
  row.names(out) <- NULL
  out
}

make_single_estimand_table <- function(paper_df, estimand_name) {
  rows <- list()
  for (sample_size in sample_sizes) {
    for (metric_name in metric_names) {
      values <- extract_metric_values(paper_df, estimand_name, sample_size, metric_name)
      rows[[paste(sample_size, metric_name, sep = "_")]] <- data.frame(
        Size = sample_size,
        Metric = ifelse(metric_name == "CR95", "CR", metric_name),
        t(values),
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  names(out) <- c("Size", "Metric", intervention_labels)
  row.names(out) <- NULL
  out
}

format_table_value <- function(x, digits = 4L) {
  value <- as.numeric(x)
  if (abs(value) < 0.5 * 10^(-digits)) value <- 0
  sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", digits, "f"), value))
}

make_latex_table <- function(integrated_table) {
  labels <- c(
    "$(0,0)$", "$(0,1)$", "$(1,0)$", "$(1,1)$", "$(A_0,0)$", "$(A_0,1)$"
  )
  lines <- c(
    "\\begin{tabular}{cc|cccccc|cccccc}",
    "\\toprule",
    "\\multirow{2}{*}{Size}",
    "& \\multirow{2}{*}{Metric}",
    "& \\multicolumn{6}{c|}{Intervention for DGP~1: Complier mean}",
    "& \\multicolumn{6}{c}{Intervention for DGP~2: Overall mean} \\\\",
    "&",
    paste0("& ", paste(c(labels, labels), collapse = " & "), " \\\\"),
    "\\midrule"
  )
  for (sample_size in sample_sizes) {
    size_df <- integrated_table[integrated_table$Size == sample_size, ]
    for (i in seq_len(nrow(size_df))) {
      row <- size_df[i, ]
      digits <- if (row$Metric == "CR") 3L else 4L
      values <- vapply(row[3:ncol(row)], format_table_value, character(1), digits = digits)
      size_cell <- if (i == 1L) paste0("\\multirow{4}{*}{", sample_size, "}") else ""
      lines <- c(lines, paste0(size_cell, " & ", row$Metric, " & ",
                               paste(values, collapse = " & "), " \\\\"))
    }
    if (sample_size != tail(sample_sizes, 1L)) lines <- c(lines, "\\addlinespace")
  }
  c(lines, "\\bottomrule", "\\end{tabular}")
}

make_single_latex_table <- function(single_table) {
  labels <- c(
    "$(0,0)$", "$(0,1)$", "$(1,0)$", "$(1,1)$", "$(A_0,0)$", "$(A_0,1)$"
  )
  lines <- c(
    "\\begin{tabular}{cc|cccccc}",
    "\\toprule",
    "\\multirow{2}{*}{Size}",
    "& \\multirow{2}{*}{Metric}",
    "& \\multicolumn{6}{c}{Intervention} \\\\",
    paste0("& & ", paste(labels, collapse = " & "), " \\\\"),
    "\\midrule"
  )
  for (sample_size in sample_sizes) {
    size_df <- single_table[single_table$Size == sample_size, ]
    for (i in seq_len(nrow(size_df))) {
      row <- size_df[i, ]
      digits <- if (row$Metric == "CR") 3L else 4L
      values <- vapply(row[3:ncol(row)], format_table_value, character(1), digits = digits)
      size_cell <- if (i == 1L) paste0("\\multirow{4}{*}{", sample_size, "}") else ""
      lines <- c(lines, paste0(size_cell, " & ", row$Metric, " & ",
                               paste(values, collapse = " & "), " \\\\"))
    }
    if (sample_size != tail(sample_sizes, 1L)) lines <- c(lines, "\\addlinespace")
  }
  c(lines, "\\bottomrule", "\\end{tabular}")
}

plot_combined_summary <- function(summary_df, method_label) {
  ggplot(summary_df, aes(x = Time)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#8ab6d6", alpha = 0.18) +
    geom_line(aes(y = Estimate, color = "Median estimate"), linewidth = 0.65) +
    geom_point(aes(y = Estimate, color = "Median estimate"), size = 1.3) +
    geom_line(aes(y = Truth, color = "Truth"), linewidth = 0.65, linetype = "dashed") +
    geom_point(aes(y = Truth, color = "Truth"), size = 1.3) +
    facet_grid(
      rows = vars(estimand, n), cols = vars(policy), scales = "free_y",
      labeller = labeller(n = function(x) paste0("n = ", x))
    ) +
    scale_x_continuous(breaks = sort(unique(summary_df$Time))) +
    scale_color_manual(values = c("Median estimate" = "#2f6f95", "Truth" = "#b23a48")) +
    labs(title = paste(method_label, "complier and overall simulation summary"),
         x = "Time", y = "ATE", color = NULL) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f2f2f2", color = "#bdbdbd"),
          legend.text = element_text(size = 16),
          strip.text = element_text(face = "bold", size = 17),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(face = "bold", size = 18))
}

plot_nudge_summary <- function(summary_df, method_label) {
  ggplot(summary_df, aes(x = Time)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#8ab6d6", alpha = 0.18) +
    geom_line(aes(y = Estimate, color = "Median estimate"), linewidth = 0.65) +
    geom_point(aes(y = Estimate, color = "Median estimate"), size = 1.3) +
    geom_line(aes(y = Truth, color = "Truth"), linewidth = 0.65, linetype = "dashed") +
    geom_point(aes(y = Truth, color = "Truth"), size = 1.3) +
    facet_grid(rows = vars(n), cols = vars(policy), scales = "free_y",
               labeller = labeller(n = function(x) paste0("n = ", x))) +
    scale_x_continuous(breaks = sort(unique(summary_df$Time))) +
    scale_color_manual(values = c("Median estimate" = "#2f6f95", "Truth" = "#b23a48")) +
    labs(title = paste(method_label, "nudge simulation summary"),
         x = "Time", y = "Nudge mean", color = NULL) +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f2f2f2", color = "#bdbdbd"),
          legend.text = element_text(size = 16),
          strip.text = element_text(face = "bold", size = 17),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(face = "bold", size = 18))
}

truth_by_estimand <- lapply(names(model_spec), make_truth)
names(truth_by_estimand) <- names(model_spec)

write_method_summary <- function(method) {
  method_label <- if (method == "xgb") "XGBoost" else toupper(method)
  complier <- make_summary(method, "complier", truth_by_estimand$complier)
  overall <- make_summary(method, "overall", truth_by_estimand$overall)
  nudge <- make_summary(method, "nudge", truth_by_estimand$nudge)

  combined <- bind_rows(complier, overall) %>%
    mutate(estimand = factor(estimand, levels = c("complier", "overall")))
  combined_table <- make_integrated_table(make_paper_table(combined))
  nudge_table <- make_single_estimand_table(make_paper_table(nudge), "nudge")

  files <- c(
    file.path(output_dir, paste0("summary_", method, "_complier_overall_table_integrated.tex")),
    file.path(output_dir, paste0("plot_", method, "_complier_overall_ggplot.pdf")),
    file.path(output_dir, paste0("summary_", method, "_nudge_table.tex")),
    file.path(output_dir, paste0("plot_", method, "_nudge_ggplot.pdf"))
  )
  writeLines(make_latex_table(combined_table), files[1])
  ggsave(files[2], plot_combined_summary(combined, method_label), width = 10, height = 8)
  writeLines(make_single_latex_table(nudge_table), files[3])
  ggsave(files[4], plot_nudge_summary(nudge, method_label), width = 10, height = 4.8)
  files
}

generated_files <- unlist(lapply(methods, write_method_summary), use.names = FALSE)
cat("Generated:\n", paste(generated_files, collapse = "\n"), "\n", sep = "")
