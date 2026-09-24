# Run from the repository root to reproduce the fixed earnings quantile map.
# The raw JTPA earnings determine the target marginal distribution; the map
# does not independently validate the fitted outcome distribution.
calibrate_jtpa_earnings_quantile_map <- function(
    output_path = "Script/par/jtpa_earnings_quantile_map.csv",
    n = 250000L, seed = 20261001L) {
  observed <- utils::read.delim("JTPA/jtpa_han.tab")$earnings
  if (!is.numeric(observed) || anyNA(observed) || any(observed < 0)) {
    stop("JTPA earnings must be complete, nonnegative numeric values.")
  }
  p <- seq(0, 1, length.out = 2001L)
  map <- data.frame(
    probability = p,
    observed_positive = as.numeric(stats::quantile(observed[observed > 0],
                                                    p, type = 8))
  )
  old_options <- options(jtpa.earnings.apply_quantile_map = FALSE)
  on.exit(options(old_options), add = TRUE)

  for (model in c("overall", "complier")) {
    source(paste0("Script/par/par_", model, "_continuous_jtpa.R"), local = TRUE)
    set.seed(seed)
    generated <- simcausal::sim(D, n = as.integer(n), wide = TRUE)$Y_2
    map[[paste0(model, "_positive")]] <-
      as.numeric(stats::quantile(generated[generated > 0], p, type = 8))
  }
  utils::write.csv(map, output_path, row.names = FALSE)
  invisible(map)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_path <- if (length(args)) args[[1L]] else
    "Script/par/jtpa_earnings_quantile_map.csv"
  calibrate_jtpa_earnings_quantile_map(output_path)
}
