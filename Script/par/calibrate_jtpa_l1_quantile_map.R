# Recreate the fixed L1 transport from the original JTPA pre-program earnings.
# Run from the repository root before regenerating the Y2 earnings map.
calibrate_jtpa_l1_quantile_map <- function(
    output_path = "Script/par/jtpa_l1_quantile_map.csv",
    n = 250000L, seed = 20261006L) {
  raw <- utils::read.delim("JTPA/jtpa_han.tab")
  observed <- as.numeric(scale(log1p(raw$prevearn)))
  floor_l1 <- min(observed)
  p <- seq(0, 1, length.out = 2001L)
  map <- data.frame(
    probability = p,
    observed_nonfloor = as.numeric(stats::quantile(
      observed[observed > floor_l1 + 1e-8], p, type = 8
    ))
  )
  old_options <- options(jtpa.l1.apply_quantile_map = FALSE,
                         jtpa.earnings.apply_quantile_map = FALSE)
  on.exit(options(old_options), add = TRUE)

  for (model in c("overall", "complier")) {
    source(paste0("Script/par/par_", model, "_continuous_jtpa.R"), local = TRUE)
    set.seed(seed)
    generated <- simcausal::sim(D, n = as.integer(n), wide = TRUE)$L_1
    map[[paste0(model, "_raw_nonfloor")]] <-
      as.numeric(stats::quantile(generated[generated > floor_l1 + 1e-8],
                                 p, type = 8))
  }
  utils::write.csv(map, output_path, row.names = FALSE)
  invisible(map)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_path <- if (length(args)) args[[1L]] else
    "Script/par/jtpa_l1_quantile_map.csv"
  calibrate_jtpa_l1_quantile_map(output_path)
}
