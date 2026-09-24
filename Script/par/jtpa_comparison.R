compare_jtpa_dgp <- function(dag, model, n = 50000L, seed = 20250919L,
                             real_path = "JTPA/jpta_han.csv",
                             output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to draw the JTPA comparison plots.")
  }
  if (!file.exists(real_path)) {
    stop("JTPA data not found: ", real_path)
  }
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }

  columns <- c("L0", "Z0", "A0", "L1", "Z1", "A1", "Y2")
  binary_columns <- columns[columns != "L1"]
  real_raw <- read.csv(real_path)
  if (!all(columns %in% names(real_raw))) {
    stop("JTPA data must contain: ", paste(columns, collapse = ", "))
  }
  real <- real_raw[columns]

  set.seed(seed)
  simulated_raw <- simcausal::sim(dag, n = as.integer(n), wide = TRUE)
  simulated <- setNames(
    simulated_raw[c("L_0", "Z_0", "A_0", "L_1", "Z_1", "A_1", "Y_2")],
    columns
  )
  if (anyNA(real) || anyNA(simulated)) {
    stop("Observed variables contain missing values.")
  }

  sources <- c("JTPA", "Simulation")
  colors <- c(JTPA = "#177E89", Simulation = "#D56A4C")
  binary_data <- do.call(rbind, lapply(binary_columns, function(variable) {
    data.frame(
      variable = variable,
      value = rep(c("0", "1"), times = 2),
      source = rep(sources, each = 2),
      proportion = c(mean(real[[variable]] == 0),
                     mean(real[[variable]] == 1),
                     mean(simulated[[variable]] == 0),
                     mean(simulated[[variable]] == 1))
    )
  }))
  binary_data$variable <- factor(binary_data$variable,
                                 levels = binary_columns)
  binary_data$source <- factor(binary_data$source, levels = sources)
  binary_plot <- ggplot2::ggplot(
    binary_data,
    ggplot2::aes(x = value, y = proportion, fill = source)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72),
                      width = 0.64) +
    ggplot2::facet_wrap(~variable, ncol = 3) +
    ggplot2::scale_fill_manual(values = colors, name = NULL) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::labs(title = paste("Binary marginals:", model),
                  x = "Value", y = "Share") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", size = 17),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 15),
                   axis.title = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(size = 18))

  l1_grid <- sort(unique(c(seq(min(c(real$L1, simulated$L1)),
                               max(c(real$L1, simulated$L1)),
                               length.out = 600),
                           min(real$L1), min(simulated$L1))))
  l1_data <- data.frame(
    L1 = rep(l1_grid, times = 2),
    cumulative = c(stats::ecdf(real$L1)(l1_grid),
                   stats::ecdf(simulated$L1)(l1_grid)),
    source = factor(rep(sources, each = length(l1_grid)), levels = sources)
  )
  l1_plot <- ggplot2::ggplot(
    l1_data, ggplot2::aes(x = L1, y = cumulative, color = source)
  ) +
    ggplot2::geom_step(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::labs(title = paste("L1 marginal distribution:", model),
                  x = "Standardized log pre-program earnings",
                  y = "Cumulative share") +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 15),
                   axis.title = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(size = 18))

  observed_cor <- stats::cor(real[columns])
  simulated_cor <- stats::cor(simulated[columns])
  correlation_gap <- simulated_cor - observed_cor
  pairs <- which(upper.tri(observed_cor), arr.ind = TRUE)
  correlations <- data.frame(
    variable_1 = columns[pairs[, 1]],
    variable_2 = columns[pairs[, 2]],
    JTPA = observed_cor[pairs],
    simulation = simulated_cor[pairs],
    difference = correlation_gap[pairs]
  )
  correlations <- correlations[order(-abs(correlations$difference)), ]
  rownames(correlations) <- NULL

  correlation_data <- expand.grid(x = columns, y = columns)
  correlation_data$difference <- as.vector(correlation_gap)
  correlation_data$x <- factor(correlation_data$x, levels = columns)
  correlation_data$y <- factor(correlation_data$y, levels = rev(columns))
  correlation_data$label <- ifelse(
    as.character(correlation_data$x) == as.character(correlation_data$y),
    "",
    ifelse(abs(correlation_data$difference) < 0.005, "0.00",
           sprintf("%+.2f", correlation_data$difference))
  )
  limit <- max(0.05, max(abs(correlation_gap)))
  correlation_plot <- ggplot2::ggplot(
    correlation_data,
    ggplot2::aes(x = x, y = y, fill = difference)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 4.5) +
    ggplot2::scale_fill_gradient2(
      low = "#B35C51", mid = "#F7F7F7", high = "#177E89",
      midpoint = 0, limits = c(-limit, limit),
      name = "Simulation - JTPA"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = paste("Pearson correlation gaps:", model),
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_text(size = 15, angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 15),
                   plot.title = ggplot2::element_text(size = 18))

  marginal_summary <- data.frame(
    metric = c(paste0("P(", binary_columns, " = 1)"),
               "L1 mean", "L1 SD", "L1 at minimum", "L1 > 1"),
    JTPA = c(vapply(real[binary_columns], mean, numeric(1)),
             mean(real$L1), stats::sd(real$L1),
             mean(real$L1 == min(real$L1)), mean(real$L1 > 1)),
    simulation = c(vapply(simulated[binary_columns], mean, numeric(1)),
                   mean(simulated$L1), stats::sd(simulated$L1),
                   mean(simulated$L1 == min(real$L1)),
                   mean(simulated$L1 > 1))
  )
  marginal_summary$difference <- marginal_summary$simulation -
    marginal_summary$JTPA

  contrast <- function(data, response, instrument) {
    group_means <- tapply(data[[response]], data[[instrument]], mean)
    unname(group_means["1"] - group_means["0"])
  }
  contrast_pairs <- list(c("A0", "Z0"), c("A1", "Z1"),
                         c("L1", "Z0"), c("Y2", "Z0"), c("Y2", "Z1"))
  contrasts <- data.frame(
    contrast = vapply(contrast_pairs, function(pair) {
      paste0(pair[1], " | ", pair[2], " = 1 versus 0")
    }, character(1)),
    JTPA = vapply(contrast_pairs, function(pair) {
      contrast(real, pair[1], pair[2])
    }, numeric(1)),
    simulation = vapply(contrast_pairs, function(pair) {
      contrast(simulated, pair[1], pair[2])
    }, numeric(1))
  )
  contrasts$difference <- contrasts$simulation - contrasts$JTPA

  cat("\n", model, "comparison: JTPA n =", nrow(real),
      ", simulation n =", nrow(simulated), "\n")
  for (table in list(marginal_summary, contrasts, correlations)) {
    display_table <- table
    numeric_columns <- vapply(display_table, is.numeric, logical(1))
    display_table[numeric_columns] <- lapply(
      display_table[numeric_columns], round, digits = 3
    )
    print(display_table, row.names = FALSE)
  }

  if (is.null(output_dir)) {
    print(binary_plot)
    print(l1_plot)
    print(correlation_plot)
  } else {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(output_dir, paste0("jtpa_", model))
    ggplot2::ggsave(paste0(prefix, "_binary_marginals.pdf"), binary_plot,
                    width = 11, height = 7, dpi = 180, bg = "white")
    ggplot2::ggsave(paste0(prefix, "_l1_ecdf.pdf"), l1_plot,
                    width = 9, height = 5, dpi = 180, bg = "white")
    ggplot2::ggsave(paste0(prefix, "_correlation_gaps.pdf"), correlation_plot,
                    width = 9, height = 7, dpi = 180, bg = "white")
    utils::write.csv(correlations, paste0(prefix, "_correlations.csv"),
                     row.names = FALSE)
    utils::write.csv(marginal_summary,
                     paste0(prefix, "_marginal_summary.csv"), row.names = FALSE)
    utils::write.csv(contrasts, paste0(prefix, "_iv_contrasts.csv"),
                     row.names = FALSE)
  }

  invisible(list(marginals = marginal_summary, correlations = correlations,
                 contrasts = contrasts,
                 plots = list(binary = binary_plot, L1 = l1_plot,
                              correlations = correlation_plot)))
}
