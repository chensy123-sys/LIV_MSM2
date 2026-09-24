# Run from the repository root: Rscript Script/par/check_jtpa_dgps.R
library(simcausal)

real <- read.csv("JTPA/jpta_han.csv")
columns <- c("L0", "Z0", "A0", "L1", "Z1", "A1", "Y2")
binary_columns <- c("L0", "Z0", "A0", "Z1", "A1")
real <- real[columns]

observed <- function(data) {
  setNames(data[c("L_0", "Z_0", "A_0", "L_1", "Z_1", "A_1", "Y_2")],
           columns)
}

joint_probabilities <- function(data) {
  code <- as.integer(as.matrix(data[binary_columns]) %*% (2^(seq_along(binary_columns) - 1L)))
  tabulate(code + 1L, nbins = 2^length(binary_columns)) / nrow(data)
}

conditional_means <- function(data) {
  c(
    as.vector(with(data, tapply(A0, list(Z0, L0), mean))),
    as.vector(with(data, tapply(A1, list(Z1, A0), mean))),
    as.vector(with(data, tapply(as.integer(L1 == min(real$L1)),
                                list(A0, L0), mean)))
  )
}

specs <- list(
  list(file = "Script/par/par_complier_continuous_jtpa.R", type = "complier"),
  list(file = "Script/par/par_overall_continuous_jtpa.R", type = "overall")
)
selected_types <- commandArgs(trailingOnly = TRUE)
if (length(selected_types) > 0L) {
  stopifnot(all(selected_types %in% c("complier", "overall")))
  specs <- Filter(function(spec) spec$type %in% selected_types, specs)
}

for (spec in specs) {
  source(spec$file)
  set.seed(20250919)
  data <- sim(D, n = 100000, wide = TRUE)
  synthetic <- observed(data)

  margin_gap <- max(abs(colMeans(synthetic[binary_columns]) - colMeans(real[binary_columns])))
  conditional_gap <- max(abs(conditional_means(synthetic) -
                             conditional_means(real)))
  joint_tv <- sum(abs(joint_probabilities(synthetic) -
                        joint_probabilities(real))) / 2
  l1_ks <- suppressWarnings(ks.test(round(synthetic$L1, 6),
                                     round(real$L1, 6))$statistic)
  floor_gap <- abs(mean(synthetic$L1 == min(real$L1)) -
                     mean(real$L1 == min(real$L1)))
  high_gap <- abs(mean(synthetic$L1 > 1) - mean(real$L1 > 1))

  cat("\n", spec$type, "\n", sep = "")
  print(round(c(max_margin_gap = margin_gap,
                max_conditional_gap = conditional_gap,
                binary_joint_tv = joint_tv,
                l1_ks_rounded = l1_ks,
                l1_floor_gap = floor_gap,
                l1_high_gap = high_gap), 4))
  stopifnot(margin_gap < 0.03, conditional_gap < 0.04,
            joint_tv < 0.05, l1_ks < 0.04,
            floor_gap < 0.02, high_gap < 0.02)

  latent_correlations <- c(U0_A0 = cor(data$U_0, data$A_0),
                           U1_A1 = cor(data$U_1, data$A_1))
  print(round(latent_correlations, 4))
  stopifnot(latent_correlations["U0_A0"] > 0.15)
  if (spec$type == "complier") {
    stopifnot(latent_correlations["U1_A1"] > 0.10)
  }

  pz0 <- plogis(0.984 + 0.04 * data$L_0)
  pz1 <- plogis(0.63 + 0.08 * data$A_0 - 0.02 * tanh(data$L_1))
  stopifnot(all(pz0 > 0.2 & pz0 < 0.8),
            all(pz1 > 0.2 & pz1 < 0.8))
  residual_correlations <- c(
    cor(data$Z_0 - pz0, data$U_0),
    cor(data$Z_1 - pz1, data$U_0),
    cor(data$Z_1 - pz1, data$U_1)
  )
  print(round(setNames(residual_correlations,
                       c("Z0_U0", "Z1_U0", "Z1_U1")), 4))
  stopifnot(max(abs(residual_correlations)) < 0.02)

  if (spec$type == "complier") {
    stopifnot(all(data$R_0 %in% c("nt", "at", "co")),
              all(data$R_1 %in% c("nt", "at", "co")),
              all(compliance_map(rep(1L, nrow(data)), data$R_0) >=
                    compliance_map(rep(0L, nrow(data)), data$R_0)),
              all(compliance_map(rep(1L, nrow(data)), data$R_1) >=
                    compliance_map(rep(0L, nrow(data)), data$R_1)))
    for (a0 in 0:1) {
      act <- node("A", t = 0, distr = "rconst", const = a0)
      counterfactual <- sim(D + action(name = "do", nodes = act),
                            n = 100000, actions = "do", wide = TRUE)$do
      joint_compliers <- mean(counterfactual$R_0 == "co" &
                                counterfactual$R_1 == "co")
      cat("Joint complier share under A0 =", a0, ":",
          round(joint_compliers, 4), "\n")
      stopifnot(joint_compliers > 0.05)
    }
  } else {
    delta0 <- deltaA(data$L_0)
    delta1 <- deltaA(data$L_1, data$A_0, stage = 1L)
    odds_product0 <- OPA(data$L_0, data$U_0)
    odds_product1 <- OPA(data$L_1, data$U_1, data$A_0, stage = 1L)
    baseline0 <- pi(data$OPA_0, data$deltaA_0)
    baseline1 <- pi(data$OPA_1, data$deltaA_1)
    stopifnot(max(abs(data$deltaA_0 - delta0)) < 1e-12,
              max(abs(data$deltaA_1 - delta1)) < 1e-12,
              max(abs(data$OPA_0 - odds_product0)) < 1e-12,
              max(abs(data$OPA_1 - odds_product1)) < 1e-12,
              all(baseline0 > 0),
              all(baseline1 > 0),
              all(baseline0 + data$deltaA_0 < 1),
              all(baseline1 + data$deltaA_1 < 1),
              min(data$deltaA_0) > 0.1,
              min(data$deltaA_1) > 0.58)
    recovered_odds0 <- baseline0 / (1 - baseline0) *
      (baseline0 + data$deltaA_0) / (1 - baseline0 - data$deltaA_0)
    recovered_odds1 <- baseline1 / (1 - baseline1) *
      (baseline1 + data$deltaA_1) / (1 - baseline1 - data$deltaA_1)
    stopifnot(max(abs(log(recovered_odds0) - log(data$OPA_0))) < 1e-10,
              max(abs(log(recovered_odds1) - log(data$OPA_1))) < 1e-10)
    cat("Conditional first-stage ranges:",
        paste(round(range(data$deltaA_0), 4), collapse = " to "), ";",
        paste(round(range(data$deltaA_1), 4), collapse = " to "), "\n")

    baseline_mean <- mean(baseline1)
    baseline_cap <- 1 - min(data$deltaA_1)
    history_component <- pz1 * data$deltaA_1
    correlation_bound <-
      (sqrt(baseline_mean * (baseline_cap - baseline_mean)) +
         sd(history_component)) / sd(data$A_1)
    cat("Current overall U1-A1 correlation bound:",
        round(correlation_bound, 4), "\n")
  }

  exclusion_mean <- numeric(2)
  exclusion_variance <- numeric(2)
  for (z in 0:1) {
    act <- c(node("Z", t = 0, distr = "rconst", const = z),
             node("Z", t = 1, distr = "rconst", const = z),
             node("A", t = 0, distr = "rconst", const = 0),
             node("A", t = 1, distr = "rconst", const = 1))
    set.seed(20250919 + z)
    counterfactual <- sim(D + action(name = "do", nodes = act),
                          n = 100000, actions = "do", wide = TRUE)$do
    exclusion_mean[z + 1] <- mean(counterfactual$Y_2)
    exclusion_variance[z + 1] <- stats::var(counterfactual$Y_2)
  }
  exclusion_tolerance <- 4 * sqrt(sum(exclusion_variance) / 100000)
  cat("Fixed-treatment IV exclusion difference:",
      round(diff(exclusion_mean), 4),
      "; Monte Carlo tolerance:", round(exclusion_tolerance, 4), "\n")
  stopifnot(abs(diff(exclusion_mean)) < exclusion_tolerance)

  l1_reduced_form <- diff(tapply(real$L1, real$Z0, mean))
  l1_reduced_form_sim <- diff(tapply(synthetic$L1, synthetic$Z0, mean))
  y_reduced_form <- diff(tapply(real$Y2, real$Z1, mean))
  y_reduced_form_sim <- diff(tapply(synthetic$Y2, synthetic$Z1, mean))
  cat("Observed L1 contrast by Z0:", round(l1_reduced_form, 4),
      "; simulated:", round(l1_reduced_form_sim, 4), "\n")
  cat("Observed Y2 contrast by Z1:", round(y_reduced_form, 4),
      "; simulated:", round(y_reduced_form_sim, 4), "\n")
  if (abs(l1_reduced_form_sim - l1_reduced_form) > 0.05 ||
      abs(y_reduced_form_sim - y_reduced_form) > 0.02) {
    cat("Reduced-form contrasts remain noticeably different from JTPA.\n")
  }
}

cat("\nStructural and prespecified similarity checks passed; see reduced-form gaps above.\n")
