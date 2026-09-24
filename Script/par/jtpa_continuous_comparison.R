# Compare a continuous-outcome JTPA DAG with observed JTPA data.
# Each run produces four plots and does not write CSV tables.
compare_jtpa_continuous_dgp <- function(dag, model, n = 50000L,
                                        seed = 20250919L,
                                        real_path = "JTPA/jpta_han.csv",
                                        output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to draw the JTPA comparison plots.")
  }
  if (!file.exists(real_path)) stop("JTPA data not found: ", real_path)
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n < 2 || n != as.integer(n)) {
    stop("n must be an integer greater than one.")
  }

  columns <- c("L0", "Z0", "A0", "L1", "Z1", "A1", "Y2")
  binary_columns <- c("L0", "Z0", "A0", "Z1", "A1")
  raw <- if (grepl("\\.tab$", real_path, ignore.case = TRUE)) {
    utils::read.delim(real_path)
  } else {
    utils::read.csv(real_path)
  }
  if (all(columns %in% names(raw))) {
    observed <- raw[columns]
  } else if (all(c("D2", "Z2", "earnings", "prevearn", "edu", "n_hs2", "sex") %in%
                 names(raw))) {
    observed <- data.frame(
      L0 = raw$sex,
      Z0 = as.numeric(raw$n_hs2 >= 13),
      A0 = as.numeric(raw$edu >= 12),
      L1 = as.numeric(scale(log1p(raw$prevearn))),
      Z1 = raw$Z2,
      A1 = raw$D2,
      Y2 = raw$earnings
    )
  } else {
    stop("JTPA data must contain processed L0,Z0,A0,L1,Z1,A1,Y2 or raw JTPA columns.")
  }

  set.seed(seed)
  simulated_raw <- simcausal::sim(dag, n = as.integer(n), wide = TRUE)
  simulated <- setNames(
    simulated_raw[c("L_0", "Z_0", "A_0", "L_1", "Z_1", "A_1", "Y_2")],
    columns
  )
  if (!all(vapply(observed, is.numeric, logical(1))) ||
      !all(vapply(simulated, is.numeric, logical(1))) ||
      anyNA(observed) || anyNA(simulated)) {
    stop("Observed and simulated variables must be numeric and complete.")
  }
  if (any(observed$Y2 < 0) || any(simulated$Y2 < 0)) {
    stop("Earnings must be nonnegative.")
  }
  if (length(unique(observed$Y2)) <= 2L) {
    stop("Observed Y2 is discrete; supply continuous earnings data.")
  }
  for (variable in binary_columns) {
    if (any(!observed[[variable]] %in% c(0, 1)) ||
        any(!simulated[[variable]] %in% c(0, 1))) {
      stop(variable, " must contain only zero and one.")
    }
  }

  sources <- c("JTPA", "Simulation")
  colors <- c(JTPA = "#177E89", Simulation = "#D56A4C")
  style <- ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 15),
                   axis.title = ggplot2::element_text(size = 18),
                   plot.title = ggplot2::element_text(size = 18),
                   plot.subtitle = ggplot2::element_text(size = 16))

  # A binary variable's mean is P(variable = 1); one bar per variable/source.
  binary_means <- do.call(rbind, lapply(binary_columns, function(variable) {
    data.frame(variable = variable, source = sources,
               mean = c(mean(observed[[variable]]),
                        mean(simulated[[variable]])))
  }))
  binary_means$variable <- factor(binary_means$variable, levels = binary_columns)
  binary_means$source <- factor(binary_means$source, levels = sources)
  binary_plot <- ggplot2::ggplot(
    binary_means, ggplot2::aes(variable, mean, fill = source)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = .72),
                      width = .64) +
    ggplot2::scale_fill_manual(values = colors, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(title = paste("Means of 0/1 variables:", model),
                  x = NULL, y = "Mean / share equal to 1") + style

  # L1 and Y2 each have a point mass. Show that mass in the subtitle and
  # compare the density of the remaining continuous values on the same scale.
  density_data <- function(jtpa, simulation) {
    if (length(jtpa) < 2L || length(simulation) < 2L) {
      stop("Each density needs at least two non-atom observations per source.")
    }
    span <- range(c(jtpa, simulation))
    bandwidth <- stats::bw.nrd0(jtpa)
    curves <- lapply(list(jtpa, simulation), function(x) {
      stats::density(x, bw = bandwidth, from = span[1], to = span[2], n = 1024)
    })
    data.frame(value = rep(curves[[1L]]$x, 2L),
               density = c(curves[[1L]]$y, curves[[2L]]$y),
               source = factor(rep(sources, each = length(curves[[1L]]$x)),
                               levels = sources))
  }
  l1_min <- min(observed$L1)
  observed_floor <- abs(observed$L1 - l1_min) < 1e-8
  simulated_floor <- abs(simulated$L1 - l1_min) < 1e-8
  l1_density <- density_data(observed$L1[!observed_floor],
                             simulated$L1[!simulated_floor])
  l1_plot <- ggplot2::ggplot(
    l1_density, ggplot2::aes(value, density, color = source)
  ) +
    ggplot2::geom_line(linewidth = .9) +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::labs(
      title = paste("L1 density above its minimum:", model),
      subtitle = sprintf("Share at minimum: JTPA %.1f%%; simulation %.1f%%",
                         100 * mean(observed_floor), 100 * mean(simulated_floor)),
      x = "Standardized log pre-program earnings", y = "Conditional density"
    ) + style

  observed_positive <- observed$Y2[observed$Y2 > 0]
  simulated_positive <- simulated$Y2[simulated$Y2 > 0]
  y2_density <- density_data(log1p(observed_positive),
                             log1p(simulated_positive))
  y2_plot <- ggplot2::ggplot(
    y2_density, ggplot2::aes(value, density, color = source)
  ) +
    ggplot2::geom_line(linewidth = .9) +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::labs(
      title = paste("Y2 positive-earnings density:", model),
      subtitle = sprintf("Zero share: JTPA %.1f%%; simulation %.1f%%",
                         100 * mean(observed$Y2 == 0),
                         100 * mean(simulated$Y2 == 0)),
      x = "log(1 + earnings in USD)", y = "Conditional density"
    ) + style

  gap <- stats::cor(simulated[columns]) - stats::cor(observed[columns])
  correlation_data <- expand.grid(x = columns, y = columns)
  correlation_data$difference <- as.vector(gap)
  correlation_data$x <- factor(correlation_data$x, levels = columns)
  correlation_data$y <- factor(correlation_data$y, levels = rev(columns))
  correlation_data$label <- ifelse(
    as.character(correlation_data$x) == as.character(correlation_data$y), "",
    ifelse(abs(correlation_data$difference) < .005, "0.00",
           sprintf("%+.2f", correlation_data$difference))
  )
  limit <- max(.05, max(abs(gap)))
  correlation_plot <- ggplot2::ggplot(
    correlation_data, ggplot2::aes(x, y, fill = difference)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = .5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 4.5) +
    ggplot2::scale_fill_gradient2(
      low = "#B35C51", mid = "#F7F7F7", high = "#177E89",
      midpoint = 0, limits = c(-limit, limit),
      breaks = c(-limit, 0, limit),
      labels = function(x) sprintf("%+.2f", x),
      name = "Simulation - JTPA",
      guide = ggplot2::guide_colorbar(
        title.position = "top", direction = "vertical",
        barwidth = grid::unit(.4, "cm"), barheight = grid::unit(4, "cm")
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = paste("Pearson correlation gaps:", model),
                  x = NULL, y = NULL) + style +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  plots <- list(correlation_gaps = correlation_plot,
                binary_means = binary_plot,
                l1_density = l1_plot,
                y2_density = y2_plot)
  if (is.null(output_dir)) {
    for (plot in plots) print(plot)
  } else {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- paste0("jtpa_", model, "_continuous_")
    filenames <- paste0(stem, names(plots), ".pdf")
    paths <- file.path(output_dir, filenames)
    for (name in names(plots)) {
      ggplot2::ggsave(
        file.path(output_dir, paste0(stem, name, ".pdf")), plots[[name]],
        width = 9,
        height = if (name == "correlation_gaps") 7 else 5,
        dpi = 180, bg = "white"
      )
    }

    # Remove only files created by earlier JTPA comparison versions.
    old_base <- c("binary_marginals.png", "l1_ecdf.png",
                  "correlation_gaps.png", "correlations.csv",
                  "marginal_summary.csv", "iv_contrasts.csv")
    old_continuous <- c(
      "binary_marginals.png", "l1_ecdf.png", "earnings_ecdf.png",
      "positive_earnings_density.png", "earnings_by_treatment.png",
      "correlation_gaps.png", "correlations.csv", "marginal_summary.csv",
      "iv_contrasts.csv", "earnings_by_treatment.csv",
      "earnings_fit_metrics.csv", "binary_means.png", "l1_density.png",
      "y2_density.png"
    )
    candidates <- c(
      file.path(output_dir, paste0("jtpa_", model, "_", old_base)),
      file.path(output_dir, paste0(stem, old_continuous))
    )
    stale <- setdiff(candidates[file.exists(candidates)], paths)
    if (length(stale) && !all(file.remove(stale))) {
      warning("Some old JTPA comparison files could not be removed.")
    }
  }
  invisible(list(plots = plots, binary_means = binary_means,
                 correlation_gaps = gap))
}
