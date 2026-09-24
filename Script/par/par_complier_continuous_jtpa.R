TM <- 1L
library(simcausal); options(simcausal.verbose = FALSE)
# Observed targets are calibrated to JTPA; latent associations are scenarios.

# JTPA earnings are in dollars, including zero earnings. The hurdle probability
# and positive-earnings distribution were calibrated to jtpa_han.tab. Treatment
# and latent-U coefficients remain simulation choices, not causal estimates.
r_jtpa_complier_earnings <- local({
  weight <- 0.3477221
  shape1 <- 0.8201464; scale1 <- 0.08127246
  shape2 <- 2.100793; scale2 <- 1
  power <- 0.610408
  moment <- weight * exp(power * log(scale1) + lgamma(shape1 + power) - lgamma(shape1)) +
    (1 - weight) * exp(power * log(scale2) + lgamma(shape2 + power) - lgamma(shape2))
  # A fixed monotone transport matches the observed positive-earnings CDF.
  # It was fitted on 250,000 observational draws (seed 20261001) and the
  # original JTPA earnings. Recalibrate after changing upstream DAG parameters.
  apply_transport <- isTRUE(getOption("jtpa.earnings.apply_quantile_map", TRUE))
  if (apply_transport) {
    map_path <- "Script/par/jtpa_earnings_quantile_map.csv"
    if (!file.exists(map_path)) stop("Earnings calibration file not found: ", map_path)
    map <- utils::read.csv(map_path)
    raw_q <- map$complier_positive
    observed_q <- map$observed_positive
    if (is.null(raw_q) || is.null(observed_q) ||
        anyNA(raw_q) || anyNA(observed_q) ||
        any(diff(raw_q) <= 0) || any(diff(observed_q) < 0)) {
      stop("Invalid earnings quantile calibration file.")
    }
    transport <- stats::approxfun(raw_q, observed_q, rule = 2, ties = "ordered")
    last <- length(raw_q)
    tail_power <- log(observed_q[last] / observed_q[last - 1L]) /
      log(raw_q[last] / raw_q[last - 1L])
  }
  function(n, zero_prob, positive_mean) {
    zero_prob <- rep_len(zero_prob, n)
    positive_mean <- rep_len(positive_mean, n)
    if (anyNA(zero_prob) || any(zero_prob < 0 | zero_prob > 1) ||
        anyNA(positive_mean) || any(positive_mean <= 0)) {
      stop("Invalid earnings distribution parameters.")
    }
    earnings <- numeric(n)
    paid <- which(runif(n) >= zero_prob)
    component1 <- runif(length(paid)) < weight
    base <- numeric(length(paid))
    base[component1] <- rgamma(sum(component1), shape = shape1, scale = scale1)
    base[!component1] <- rgamma(sum(!component1), shape = shape2, scale = scale2)
    raw_earnings <- positive_mean[paid] * base^power / moment
    if (!apply_transport) {
      earnings[paid] <- raw_earnings
      return(earnings)
    }
    calibrated <- transport(raw_earnings)
    upper <- raw_earnings > raw_q[last]
    calibrated[upper] <- observed_q[last] *
      (raw_earnings[upper] / raw_q[last])^tail_power
    earnings[paid] <- calibrated
    earnings
  }
})
floor_l1 <- -1.52409977806294

# The floor probability stays structural. A monotone map calibrates only the
# continuous part of L1; the observed quantiles reproduce earnings heaping.
calibrate_l1_complier <- local({
  apply_transport <- isTRUE(getOption("jtpa.l1.apply_quantile_map", TRUE))
  if (apply_transport) {
    path <- "Script/par/jtpa_l1_quantile_map.csv"
    if (!file.exists(path)) stop("L1 calibration file not found: ", path)
    map <- utils::read.csv(path)
    raw_q <- map$complier_raw_nonfloor
    observed_q <- map$observed_nonfloor
    if (is.null(raw_q) || is.null(observed_q) || anyNA(raw_q) ||
        anyNA(observed_q) || any(diff(raw_q) < 0) ||
        any(diff(observed_q) < 0)) stop("Invalid L1 calibration file.")
    transport <- stats::approxfun(raw_q, observed_q, rule = 2, ties = "ordered")
    last <- length(raw_q)
    upper_slope <- (observed_q[last] - observed_q[last - 1L]) /
      (raw_q[last] - raw_q[last - 1L])
  }
  function(value) {
    if (!apply_transport) return(value)
    result <- value
    nonfloor <- value > floor_l1 + 1e-8
    result[nonfloor] <- transport(value[nonfloor])
    upper <- nonfloor & value > raw_q[last]
    result[upper] <- observed_q[last] + upper_slope * (value[upper] - raw_q[last])
    result
  }
})

compliance_map <- function(Z, R) {
  if (length(Z) != length(R)) {
    stop("Instrument and compliance vectors must have equal length.")
  }
  A <- rep(NA_integer_, length(Z))
  A[R == "nt"] <- 0L
  A[R == "at"] <- 1L
  A[R == "co"] <- Z[R == "co"]
  if (anyNA(A)) {
    stop("Unknown compliance type.")
  }
  A
}

D <- DAG.empty() +
  node("L", t = 0, distr = "rbern", prob = 0.465) +
  node("U", t = 0, distr = "rnorm", mean = 0, sd = 1) +
  node("Z", t = 0, distr = "rbern", prob = plogis(0.984 + 0.040 * L[t])) +
  node("scoreCo", t = 0, distr = "rconst",
       const = plogis(-2.100 - 0.200 * L[t] + 0.750 * U[t])) +
  node("scoreAt", t = 0, distr = "rconst",
       const = (1 - scoreCo[t]) * plogis(0.430 + 0.450 * U[t] + 0.020 * L[t])) +
  node("Rtype", t = 0, distr = "rcat.b1",
       probs = cbind(scoreCo[t], scoreAt[t], 1 - scoreCo[t] - scoreAt[t])) +
  node("R", t = 0, distr = "rconst",
       const = ifelse(Rtype[t] == 1, "co",
                      ifelse(Rtype[t] == 2, "at", "nt"))) +
  node("A", t = 0, distr = "rconst", const = compliance_map(Z[t], R[t])) +

  node("Ufloor", t = 1:TM, distr = "rbern",
       prob = plogis(-0.349 - 0.481 * A[t-1] - 0.712 * L[t-1] - 0.250 * U[t-1])) +
  node("Udraw", t = 1:TM, distr = "rlnorm",
       meanlog = log(0.595), sdlog = 0.530) +
  node("L", t = 1:TM, distr = "rconst",
       const = calibrate_l1_complier(
         ifelse(Ufloor[t] == 1, floor_l1,
                pmin(1.600, pmax(floor_l1 + 0.001,
                                 1.150 + 0.090 * A[t-1] + 0.100 * L[t-1] +
                                   0.050 * U[t-1] - Udraw[t]))))) +
  node("U", t = 1:TM, distr = "rnorm", mean = 0.400 * U[t-1], sd = 0.900) +
  node("Z", t = 1:TM, distr = "rbern",
       prob = plogis(0.630 + 0.080 * A[t-1] - 0.020 * tanh(L[t]))) +
  node("scoreCo", t = 1:TM, distr = "rconst",
       const = plogis(0.500 + 0.250 * A[t-1] + 0.050 * L[t] + 0.500 * U[t])) +
  node("scoreAt", t = 1:TM, distr = "rconst",
       const = (1 - scoreCo[t]) * plogis(-3.200 + 0.400 * U[t])) +
  node("Rtype", t = 1:TM, distr = "rcat.b1",
       probs = cbind(scoreCo[t], scoreAt[t], 1 - scoreCo[t] - scoreAt[t])) +
  node("R", t = 1:TM, distr = "rconst",
       const = ifelse(Rtype[t] == 1, "co",
                      ifelse(Rtype[t] == 2, "at", "nt"))) +
  node("A", t = 1:TM, distr = "rconst", const = compliance_map(Z[t], R[t])) +
  node("Y", t = TM + 1L, distr = "r_jtpa_complier_earnings",
       zero_prob = plogis(-1.879364 - 0.18 * A[t-2] - 0.35 * A[t-1] -
                           0.10 * L[t-2] - 0.55 * L[t-1] - 0.20 * U[t-1]),
       positive_mean = exp(9.418864 + 0.20 * A[t-2] + 0.12 * A[t-1] +
                             0.30 * L[t-2] + 0.20 * L[t-1] -
                             0.08 * A[t-2] * A[t-1] + 0.25 * U[t-1]))

D <- set.DAG(D, vecfun = c("cbind", "compliance_map", "calibrate_l1_complier",
                           "pmax", "pmin"))

# Run from the repository root after sourcing this file.
# plot_complier_continuous_jtpa_fit(output_dir = "outputs/jtpa_complier")
plot_complier_continuous_jtpa_fit <- local({
  model_dag <- D
  function(n = 100000L, seed = 20250919L,
           real_path = "JTPA/jpta_han.csv", output_dir = NULL) {
    source("Script/par/jtpa_continuous_comparison.R", local = TRUE)
    compare_jtpa_continuous_dgp(model_dag, "complier", n, seed,
                                real_path, output_dir)
  }
})

# Keep the call used by the original complier parameter script.
plot_complier_jtpa_fit <- plot_complier_continuous_jtpa_fit
