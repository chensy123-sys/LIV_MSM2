source("Script/par/par_overall.R")
source("LIV_xgboost.R")

set.seed(2025)
Looptime <- 500
sample_sizes <- c(2000, 5000)

for (n in sample_sizes) {
  output_dir <- file.path("data", paste0("xgb", n))
  dir.create(output_dir, showWarnings = FALSE)

  RES <- array(0, dim = c(TM + 2, 5, Looptime))
  for (loop in seq_len(Looptime)) {
    print(paste("n =", n, "loop =", loop))
    data_LIV <- sim_LIV_overall(n, D)
    fit_LIV <- DML_LIV(mydat = data_LIV, action_list = c(1, 0), estimand = "overall")
    RES[, , loop] <- as.matrix(predict_LIV(fit_LIV))
  }
  save(RES, file = file.path(output_dir, "Overall_10.RData"))
}
