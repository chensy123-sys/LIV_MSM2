# JTPA summary figures and tables for all 12 target outcomes.
#
# The first eight targets are full two-period interventions.  x0, x1, xp and
# xn are the second-period entries from the corresponding 00, 01, 0p and 0n
# fits, respectively.  At that point A_0 remains natural, so they estimate
# Y(A_0, 0), Y(A_0, 1), Y(A_0, p) and Y(A_0, n).

# Do not clear .GlobalEnv: doing so can interrupt non-interactive Rscript
# execution before the remaining expressions are evaluated.

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

sample_sizes <- 10000L
# Override only for a quick diagnostic run, e.g. JTPA_TRUTH_ACCURACY=100000.
truth_accuracy <- as.integer(Sys.getenv("JTPA_TRUTH_ACCURACY", "1000000"))
if (is.na(truth_accuracy) || truth_accuracy < 1L) {
  stop("JTPA_TRUTH_ACCURACY must be a positive integer.", call. = FALSE)
}
set.seed(2025)
threshold <- 0.8368623
policy_labels <- c("00", "01", "0p", "0n", "10", "11", "1p", "1n")
dynamic_policy_labels <- c("x0", "x1", "xp", "xn")
dgp_labels <- c("00", "0n", "0p", "01", "10", "1n", "1p", "11", "x0", "xn", "xp", "x1")

# The result vector has one entry for each intervention start time.  Time 0
# fixes both A_0 and A_1; Time 1 fixes only A_1 and leaves A_0 natural.
dynamic_policy_source <- c("00" = "x0", "01" = "x1", "0p" = "xp", "0n" = "xn")
dgp_math_labels <- c(
  "00" = "00", "01" = "01",
  "0p" = "0d+", "0n" = "0d-",
  "10" = "10", "11" = "11",
  "1p" = "1d+", "1n" = "1d-",
  "x0" = "x0", "x1" = "x1",
  "xp" = "xd+", "xn" = "xd-"
)
dgp_plotmath_labels <- c(
  "00" = "paste(0,0)",
  "01" = "paste(0,1)",
  "0p" = "paste(0,d^\"+\")",
  "0n" = "paste(0,d^\"-\")",
  "10" = "paste(1,0)",
  "11" = "paste(1,1)",
  "1p" = "paste(1,d^\"+\")",
  "1n" = "paste(1,d^\"-\")",
  "x0" = "paste(x,0)",
  "x1" = "paste(x,1)",
  "xp" = "paste(x,d^\"+\")",
  "xn" = "paste(x,d^\"-\")"
)
dgp_latex_labels <- c(
  "00" = "00", "01" = "01", "0p" = "0d^+", "0n" = "0d^-",
  "10" = "10", "11" = "11", "1p" = "1d^+", "1n" = "1d^-",
  "x0" = "x0", "x1" = "x1", "xp" = "xd^+", "xn" = "xd^-"
)

acts <- list(
  "00" = c(0, 0),
  "01" = c(0, 1),
  "0p" = list(function(fY) 0, function(fY) as.numeric(fY$L < threshold)),
  "0n" = list(function(fY) 0, function(fY) as.numeric(fY$L >= threshold)),
  "10" = c(1, 0),
  "11" = c(1, 1),
  "1p" = list(function(fY) 1, function(fY) as.numeric(fY$L < threshold)),
  "1n" = list(function(fY) 1, function(fY) as.numeric(fY$L >= threshold))
)

output_dir <- file.path("outputs", "jtpa")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

result_dir <- function(n) {
  preferred <- file.path("data", paste0("jtpa", n))
  alternate <- paste0("data_jtpa", n)
  if (dir.exists(preferred)) preferred else alternate
}

load_result <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("RES", envir = env)) stop("Missing RES in ", path, call. = FALSE)
  env$RES
}

source_model <- function(estimand) {
  if (estimand == "overall") {
    quiet_source("Script/par/par_overall_continuous_jtpa.R")
    truth_fun <- ATE_overall
  } else if (estimand == "complier") {
    quiet_source("Script/par/par_complier_continuous_jtpa.R")
    truth_fun <- ATE_local
  } else {
    stop("Unknown estimand: ", estimand, call. = FALSE)
  }
  list(D = D, TM = TM, truth_fun = truth_fun)
}

make_truth <- function(estimand) {
  model <- source_model(estimand)
  truth <- lapply(acts, function(act) {
    suppressMessages(model$truth_fun(
      model$D, model$TM, PROB = act, accuracy = truth_accuracy
    ))
  })
  list(TM = model$TM, truth = truth)
}

summarize_one <- function(estimand, n, policy, truth) {
  file_prefix <- if (estimand == "overall") "Overall" else "Complier"
  path <- file.path(result_dir(n), paste0(file_prefix, "_", policy, ".RData"))
  if (!file.exists(path)) stop("Result file not found: ", path, call. = FALSE)

  RES <- load_result(path)
  times <- seq_along(truth)
  if (length(dim(RES)) != 3L || dim(RES)[2] < 5L || dim(RES)[1] != length(times)) {
    stop(
      "Unexpected RES dimensions in ", path, "; expected ", length(times),
      " x at least 5 x repetitions, got ", paste(dim(RES), collapse = " x "),
      call. = FALSE
    )
  }

  estimate <- apply(RES[, 2, ], 1, median, na.rm = TRUE)
  data.frame(
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

make_summary <- function(estimand, truth_obj) {
  bind_rows(lapply(sample_sizes, function(n) {
    bind_rows(lapply(policy_labels, function(policy) {
      summarize_one(estimand, n, policy, truth_obj$truth[[policy]])
    }))
  }))
}

# This mirrors the former table's (A_0, 0) and (A_0, 1) rows, now retaining
# all four natural-A_0 targets and all eight full interventions.
make_dgp_summary <- function(summary_df) {
  fixed <- summary_df %>%
    filter(Time == 0) %>%
    mutate(DGP = policy)
  natural_a0 <- summary_df %>%
    filter(Time == 1, policy %in% names(dynamic_policy_source)) %>%
    mutate(DGP = unname(dynamic_policy_source[policy]))

  bind_rows(fixed, natural_a0) %>%
    mutate(
      DGP = factor(DGP, levels = dgp_labels),
      Outcome = unname(dgp_math_labels[as.character(DGP)])
    ) %>%
    arrange(estimand, n, DGP)
}

metric_names <- c("Bias", "SD", "SE", "CR95")

extract_metric_values <- function(dgp_df, estimand_name, sample_size, metric) {
  values <- dgp_df %>%
    filter(.data$estimand == estimand_name, n == sample_size) %>%
    arrange(DGP)
  if (nrow(values) != length(dgp_labels) || anyDuplicated(values$DGP)) {
    stop("Expected one row for each of the 12 DGPs for ", estimand_name,
         ", n = ", sample_size, call. = FALSE)
  }
  result <- values[[metric]]
  names(result) <- as.character(values$DGP)
  result[dgp_labels]
}

make_integrated_table <- function(dgp_df) {
  rows <- lapply(sample_sizes, function(sample_size) {
    bind_rows(lapply(metric_names, function(metric) {
      values <- c(
        extract_metric_values(dgp_df, "complier", sample_size, metric),
        extract_metric_values(dgp_df, "overall", sample_size, metric)
      )
      value_frame <- as.data.frame(as.list(values), check.names = FALSE)
      names(value_frame) <- paste0("value_", seq_along(values))
      data.frame(
        Size = sample_size,
        Metric = ifelse(metric == "CR95", "CR", metric),
        value_frame,
        check.names = FALSE
      )
    }))
  })
  table <- bind_rows(rows)
  names(table) <- c(
    "Size", "Metric", paste0("Complier_", dgp_labels),
    paste0("Overall_", dgp_labels)
  )
  table
}

format_table_value <- function(x, digits = 4L) {
  value <- as.numeric(x)
  if (abs(value) < 0.5 * 10^(-digits)) value <- 0
  sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", digits, "f"), value))
}

make_latex_table <- function(dgp_df) {
  if (length(sample_sizes) != 1L) {
    stop("The JTPA table layout expects exactly one sample size.", call. = FALSE)
  }
  sample_size <- sample_sizes[[1L]]
  labels <- paste0("$", dgp_latex_labels[dgp_labels], "$")
  columns <- paste(c("cc|", rep("c", length(dgp_labels))), collapse = "")
  lines <- c(
    paste0("\\begin{tabular}{", columns, "}"),
    "\\toprule",
    paste0("\\multirow{2}{*}{Estimand} & \\multirow{2}{*}{Metric} & \\multicolumn{", length(dgp_labels),
           "}{c}{Intervention} \\\\"),
    paste0("& & ", paste(labels, collapse = " & "), " \\\\"),
    "\\midrule"
  )
  for (estimand_name in c("complier", "overall")) {
    if (estimand_name == "overall") lines <- c(lines, "\\midrule")
    estimand_label <- if (estimand_name == "overall") "Overall" else "Complier"
    for (metric_name in metric_names) {
      display_metric <- if (metric_name == "CR95") "CR" else metric_name
      digits <- if (display_metric %in% c("Bias", "SD", "SE")) {
        2L
      } else if (display_metric == "CR") {
        3L
      } else {
        4L
      }
      values <- extract_metric_values(
        dgp_df, estimand_name, sample_size, metric_name
      )
      values <- vapply(values, format_table_value, character(1), digits = digits)
      estimand_cell <- if (metric_name == metric_names[[1L]]) {
        paste0("\\multirow{4}{*}{", estimand_label, "}")
      } else {
        ""
      }
      lines <- c(lines, paste0(estimand_cell, " & ", display_metric, " & ",
                               paste(values, collapse = " & "), " \\\\"))
    }
  }
  c(lines, "\\bottomrule", "\\end{tabular}")
}

plot_jtpa_summary <- function(dgp_df) {
  ggplot(dgp_df, aes(DGP, Estimate)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.35) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.18,
                  color = "#2f6f95") +
    geom_point(color = "#2f6f95", size = 2.1) +
    geom_point(aes(y = Truth, shape = "Truth"), color = "#b23a48", size = 2.3) +
    facet_grid(estimand ~ ., scales = "free_y",
               labeller = labeller(estimand = c(
                 overall = "Overall mean", complier = "Complier mean"
               ))) +
    scale_shape_manual(values = c(Truth = 4), name = NULL) +
    scale_x_discrete(labels = function(x) parse(text = dgp_plotmath_labels[x])) +
    labs(
      x = NULL, y = "Mean earnings (USD)"
    ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = c(0.96, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 16),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#f2f2f2", color = "#bdbdbd"),
      strip.text = element_text(face = "bold", size = 17),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 18)
    )
}

truth_by_estimand <- lapply(c("overall", "complier"), make_truth)
names(truth_by_estimand) <- c("overall", "complier")

summary_df <- bind_rows(
  make_summary("complier", truth_by_estimand$complier),
  make_summary("overall", truth_by_estimand$overall)
) %>%
  mutate(estimand = factor(estimand, levels = c("complier", "overall")))
dgp_summary <- make_dgp_summary(summary_df)

latex_path <- file.path(output_dir, "jtpa_12_dgp_table.tex")
plot_path <- file.path(output_dir, "jtpa_12_dgp_ggplot.pdf")
writeLines(make_latex_table(dgp_summary), latex_path)
ggsave(plot_path, plot_jtpa_summary(dgp_summary), width = 14, height = 6.5)

cat("Generated:\n", paste(c(latex_path, plot_path), collapse = "\n"), "\n", sep = "")
