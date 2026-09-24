library(simcausal)
library(tibble)
library(dplyr)
library(mgcv)
library(purrr)
library(ggplot2)

# Shared utilities -------------------------------------------------------------

.liv_sim_observed <- function(n, D, drop_patterns) {
  mydat <- tibble(sim(D, n = n, wide = FALSE))
  for (pattern in drop_patterns) {
    mydat <- mydat %>% select(-contains(pattern))
  }
  return(mydat)
}

.liv_make_id_folds <- function(mydat, id_var = "ID", K = 5, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ids <- sort(unique(mydat[[id_var]]))
  folds <- sample(rep(1:K, length.out = length(ids)))
  id_folds <- split(ids, folds)

  lapply(seq_along(id_folds), function(k) {
    test_ids <- id_folds[[k]]
    train_ids <- setdiff(ids, test_ids)
    list(train_ids = train_ids, test_ids = test_ids)
  })
}

.liv_build_fY_history <- function(mydat, time, H) {
  fY <- mydat %>% filter(t == time) %>% arrange(ID) %>% select(L, Z, A)

  h <- 1
  while (h <= H - 1) {
    past_time <- time - h
    if (past_time < 0) {
      break
    }
    temp <- mydat %>%
      filter(t == past_time) %>%
      arrange(ID) %>%
      select(L, Z, A) %>%
      rename_with(~ paste0(., "_lag", h))
    fY <- bind_cols(fY, temp)
    h <- h + 1
  }

  return(fY)
}

.liv_build_fA_history <- function(mydat, time, H) {
  fA <- mydat %>% filter(t == time) %>% arrange(ID) %>% select(L, Z)

  h <- 1
  while (h <= H - 1) {
    past_time <- time - h
    if (past_time < 0) {
      break
    }
    temp <- mydat %>%
      filter(t == past_time) %>%
      arrange(ID) %>%
      select(L, Z) %>%
      rename_with(~ paste0(., "_lag", h))
    fA <- bind_cols(fA, temp)
    h <- h + 1
  }

  h <- 1
  while (h <= H) {
    past_time <- time - h
    if (past_time < 0) {
      break
    }
    temp <- mydat %>%
      filter(t == past_time) %>%
      arrange(ID) %>%
      select(A) %>%
      rename_with(~ paste0(., "_lag", h))
    fA <- bind_cols(fA, temp)
    h <- h + 1
  }

  return(fA)
}

.liv_build_fZ_history <- function(mydat, time, H) {
  fZ <- mydat %>% filter(t == time) %>% arrange(ID) %>% select(L)

  h <- 1
  while (h <= H - 1) {
    past_time <- time - h
    if (past_time < 0) {
      break
    }
    temp <- mydat %>%
      filter(t == past_time) %>%
      arrange(ID) %>%
      select(L) %>%
      rename_with(~ paste0(., "_lag", h))
    fZ <- bind_cols(fZ, temp)
    h <- h + 1
  }

  h <- 1
  while (h <= H) {
    past_time <- time - h
    if (past_time < 0) {
      break
    }
    temp <- mydat %>%
      filter(t == past_time) %>%
      arrange(ID) %>%
      select(Z, A) %>%
      rename_with(~ paste0(., "_lag", h))
    fZ <- bind_cols(fZ, temp)
    h <- h + 1
  }

  return(fZ)
}

.liv_get_treatment <- function(action_list, time, fY) {
  action <- action_list[[time + 1]]
  if (is.function(action)) {
    return(action(fY))
  }
  action_list[time + 1]
}

.liv_mgcv_variable_types <- function(X, cts_num = 4) {
  continuous <- vapply(X, function(x) is.numeric(x) && length(unique(x)) > cts_num, logical(1))
  list(
    continuous = names(X)[continuous],
    discrete = names(X)[!continuous]
  )
}

.liv_mgcv_family_from_response <- function(Y) {
  y <- Y[!is.na(Y)]
  if (length(y) == 0) {
    return(gaussian())
  }

  if (is.logical(y)) {
    return(binomial())
  }

  if (is.numeric(y) || is.integer(y)) {
    values <- unique(as.numeric(y))
    if (length(values) <= 2 && all(values %in% c(0, 1))) {
      return(binomial())
    }
    return(gaussian())
  }

  values <- unique(tolower(as.character(y)))
  if (length(values) <= 2 && all(values %in% c("0", "1", "false", "true"))) {
    return(binomial())
  }

  gaussian()
}

.liv_mgcv_constant_fit <- function(Y, family) {
  list(
    type = "constant",
    value = mean(as.numeric(Y)),
    family = family$family
  )
}

.liv_mgcv_prepare_discrete <- function(data, discrete_vars, levels = NULL) {
  if (length(discrete_vars) == 0) {
    return(list(data = data, levels = levels))
  }

  if (is.null(levels)) {
    levels <- lapply(discrete_vars, function(var) {
      vals <- unique(data[[var]])
      if (is.numeric(data[[var]]) || is.integer(data[[var]]) || is.logical(data[[var]])) {
        vals <- unique(c(vals, 0, 1))
      }
      sort(as.character(vals))
    })
    names(levels) <- discrete_vars
  }

  for (var in discrete_vars) {
    data[[var]] <- factor(as.character(data[[var]]), levels = levels[[var]])
  }

  data$.disc_stratum <- interaction(data[, discrete_vars, drop = FALSE], drop = FALSE, sep = ":")
  list(data = data, levels = levels)
}

.liv_mgcv_stratum_fit <- function(Y, X, family, discrete_vars, levels) {
  prepared <- .liv_mgcv_prepare_discrete(as.data.frame(X), discrete_vars, levels)
  Xp <- prepared$data
  y <- as.numeric(Y)

  if (family$family == "binomial") {
    global <- mean(y)
    means <- tapply(y, Xp$.disc_stratum, mean)
  } else {
    global <- mean(y)
    means <- tapply(y, Xp$.disc_stratum, mean)
  }

  list(
    type = "stratum",
    means = means,
    global = global,
    family = family$family,
    discrete_vars = discrete_vars,
    discrete_levels = prepared$levels
  )
}

.liv_mgcv_basis_k <- function(train_data, continuous_vars, discrete_vars, use_by) {
  d <- length(continuous_vars)
  min_k <- d + 2
  unique_cont <- nrow(unique(train_data[, continuous_vars, drop = FALSE]))
  local_n <- nrow(train_data)

  if (use_by && length(discrete_vars) > 0) {
    observed <- table(train_data$.disc_stratum)
    observed <- observed[observed > 0]
    local_n <- min(as.integer(observed))
  }

  upper_k <- min(10, unique_cont - 1, local_n - 1)
  if (upper_k < min_k) {
    return(NULL)
  }
  max(min_k, upper_k)
}

.liv_mgcv_smooth_term <- function(continuous_vars, k = NULL, use_by = FALSE) {
  args <- continuous_vars
  if (use_by) {
    args <- c(args, "by = .disc_stratum")
  }
  if (!is.null(k)) {
    args <- c(args, paste0("k = ", k))
  }
  paste0("s(", paste(args, collapse = ", "), ")")
}

.liv_mgcv_build_formula <- function(continuous_vars, discrete_vars, use_by = TRUE, k = NULL) {
  if (length(continuous_vars) == 0 && length(discrete_vars) == 0) {
    return(Y ~ 1)
  }

  if (length(continuous_vars) == 0) {
    return(Y ~ .disc_stratum)
  }

  smooth_term <- .liv_mgcv_smooth_term(continuous_vars, k = k, use_by = FALSE)
  if (length(discrete_vars) == 0) {
    return(as.formula(paste("Y ~", smooth_term)))
  }

  if (use_by) {
    smooth_term <- .liv_mgcv_smooth_term(continuous_vars, k = k, use_by = TRUE)
  }
  as.formula(paste("Y ~ .disc_stratum +", smooth_term))
}

.liv_mgcv_parametric_formula <- function(continuous_vars, discrete_vars) {
  terms <- c()
  if (length(discrete_vars) > 0) {
    terms <- c(terms, ".disc_stratum")
  }
  terms <- c(terms, continuous_vars)
  as.formula(paste("Y ~", paste(terms, collapse = " + ")))
}

.liv_mgcv_try_gam <- function(gam_formula, train_data, family) {
  tryCatch(
    mgcv::gam(
      gam_formula,
      data = train_data,
      family = family,
      method = "REML",
      select = TRUE,
      drop.unused.levels = FALSE
    ),
    error = function(e) NULL
  )
}

.liv_mgcv_fit <- function(Y, X, family = NULL, cts_num = 4) {
  if (is.null(family)) {
    family <- .liv_mgcv_family_from_response(Y)
  }

  X <- as.data.frame(X)
  if (length(unique(Y)) == 1) {
    return(.liv_mgcv_constant_fit(Y, family))
  }

  types <- .liv_mgcv_variable_types(X, cts_num = cts_num)
  continuous_vars <- types$continuous
  discrete_vars <- types$discrete

  if (length(continuous_vars) == 0) {
    return(.liv_mgcv_stratum_fit(Y, X, family, discrete_vars, levels = NULL))
  }

  prepared <- .liv_mgcv_prepare_discrete(X, discrete_vars)
  train_data <- prepared$data
  train_data$Y <- Y
  k_by <- .liv_mgcv_basis_k(train_data, continuous_vars, discrete_vars, use_by = TRUE)
  k_pooled <- .liv_mgcv_basis_k(train_data, continuous_vars, discrete_vars, use_by = FALSE)

  fit <- NULL
  if (length(discrete_vars) == 0 || !is.null(k_by)) {
    gam_formula <- .liv_mgcv_build_formula(continuous_vars, discrete_vars, use_by = TRUE, k = k_by)
    fit <- .liv_mgcv_try_gam(gam_formula, train_data, family)
  }

  if (is.null(fit) && length(discrete_vars) > 0 && !is.null(k_pooled)) {
    gam_formula <- .liv_mgcv_build_formula(continuous_vars, discrete_vars, use_by = FALSE, k = k_pooled)
    fit <- .liv_mgcv_try_gam(gam_formula, train_data, family)
  }

  if (is.null(fit)) {
    gam_formula <- .liv_mgcv_parametric_formula(continuous_vars, discrete_vars)
    fit <- .liv_mgcv_try_gam(gam_formula, train_data, family)
  }

  if (is.null(fit)) {
    return(.liv_mgcv_constant_fit(Y, family))
  }

  list(
    type = "mgcv",
    object = fit,
    family = family$family,
    fallback = mean(as.numeric(Y)),
    continuous_vars = continuous_vars,
    discrete_vars = discrete_vars,
    discrete_levels = prepared$levels
  )
}

.liv_mgcv_predict <- function(fit, newdata) {
  newdata <- as.data.frame(newdata)

  if (fit$type == "constant") {
    pred <- rep(fit$value, nrow(newdata))
  } else if (fit$type == "stratum") {
    prepared <- .liv_mgcv_prepare_discrete(newdata, fit$discrete_vars, fit$discrete_levels)
    pred <- as.numeric(fit$means[as.character(prepared$data$.disc_stratum)])
    pred[is.na(pred)] <- fit$global
  } else {
    prepared <- .liv_mgcv_prepare_discrete(newdata, fit$discrete_vars, fit$discrete_levels)
    pred <- tryCatch(
      as.numeric(predict(fit$object, newdata = prepared$data, type = "response")),
      error = function(e) rep(fit$fallback, nrow(newdata))
    )
    pred[!is.finite(pred)] <- fit$fallback
  }

  if (fit$family == "binomial") {
    pred <- pmin(pmax(pred, 1e-6), 1 - 1e-6)
  }
  pred
}

.liv_fold_error <- function(k, error) {
  structure(
    list(fold = k, message = conditionMessage(error)),
    class = "liv_fold_error"
  )
}

.liv_stop_on_fold_errors <- function(results) {
  fold_failed <- vapply(results, inherits, logical(1), what = "liv_fold_error")
  fork_failed <- vapply(results, inherits, logical(1), what = "try-error")
  failed <- fold_failed | fork_failed

  if (any(failed)) {
    messages <- vapply(seq_along(results)[failed], function(i) {
      x <- results[[i]]
      if (inherits(x, "liv_fold_error")) {
        return(paste0("fold ", x$fold, ": ", x$message))
      }
      condition <- attr(x, "condition")
      msg <- if (is.null(condition)) as.character(x)[1] else conditionMessage(condition)
      paste0("fold ", i, ": ", msg)
    }, character(1))
    stop(paste(messages, collapse = "\n"), call. = FALSE)
  }
  results
}

.liv_describe_object <- function(x) {
  object_class <- paste(class(x), collapse = "/")
  object_names <- names(x)
  if (is.null(object_names)) {
    object_names <- "<none>"
  } else {
    object_names <- paste(object_names, collapse = ", ")
  }
  paste0("class=", object_class, "; names=", object_names)
}

.liv_validate_fold_matrix <- function(x, n, p, name, fold) {
  expected_dim <- c(n, p)
  if (!is.matrix(x) || length(dim(x)) != 2 || any(dim(x) != expected_dim)) {
    stop(
      sprintf(
        "%s from fold %s has dimension %s; expected %s x %s.",
        name,
        fold,
        paste(dim(x), collapse = " x "),
        n,
        p
      ),
      call. = FALSE
    )
  }
}

.liv_validate_local_fold_result <- function(x, n, p, fold) {
  if (inherits(x, "liv_fold_error")) {
    stop(paste0("fold ", x$fold, ": ", x$message), call. = FALSE)
  }
  if (inherits(x, "try-error")) {
    condition <- attr(x, "condition")
    msg <- if (is.null(condition)) as.character(x)[1] else conditionMessage(condition)
    stop(paste0("fold ", fold, ": ", msg), call. = FALSE)
  }
  if (!is.list(x) || is.null(x$ETA) || is.null(x$ETA_den)) {
    stop(
      paste0(
        "fold ", fold, " returned an invalid result; expected list(ETA, ETA_den), got ",
        .liv_describe_object(x), "."
      ),
      call. = FALSE
    )
  }
  .liv_validate_fold_matrix(x$ETA, n, p, "ETA", fold)
  .liv_validate_fold_matrix(x$ETA_den, n, p, "ETA_den", fold)
}

.liv_predict_point <- function(fit_LIV, time, local = FALSE) {
  id_order <- if (!is.null(fit_LIV$id_order)) {
    fit_LIV$id_order
  } else {
    sort(unique(fit_LIV$input_data$ID))
  }
  n <- length(id_order)
  K <- length(fit_LIV$folds)
  ETA_merge <- rep(0, n)

  if (local) {
    ETA_den_merge <- rep(0, n)
  }

  for (k in 1:K) {
    index <- match(fit_LIV$folds[[k]]$test_ids, id_order)
    if (anyNA(index)) {
      stop("Fold IDs do not match the fitted ID order.")
    }
    ETA_merge[index] <- fit_LIV$ETA[index, k, time + 1]
    if (local) {
      ETA_den_merge[index] <- fit_LIV$ETA_den[index, k, time + 1]
    }
  }

  if (local) {
    ATE <- mean(ETA_merge) / mean(ETA_den_merge)
    SD <- sd(ETA_merge - ATE * ETA_den_merge) / mean(ETA_den_merge) / sqrt(n)
  } else {
    ATE <- mean(ETA_merge)
    SD <- sd(ETA_merge) / sqrt(n)
  }

  list(time = time, est = ATE, std = SD)
}

.liv_predict <- function(fit_LIV, local = FALSE) {
  res <- unlist(.liv_predict_point(fit_LIV, time = 0, local = local))
  for (time in 1:(dim(fit_LIV$ETA)[3] - 1)) {
    res <- rbind(res, unlist(.liv_predict_point(fit_LIV, time, local = local)))
  }
  res <- as.data.frame(res)
  rownames(res) <- NULL
  res <- res %>% mutate(lower = est - 1.96 * std, upper = est + 1.96 * std)
  return(res)
}

.liv_generate_action_from_vector <- function(treatment_vec, varname = "A") {
  nodes <- list()
  for (i in 1:length(treatment_vec)) {
    temp <- treatment_vec[i]
    if (temp == 1) {
      nodes <- c(nodes, node(varname, t = i - 1, distr = "rconst", const = 1))
    } else {
      nodes <- c(nodes, node(varname, t = i - 1, distr = "rconst", const = 0))
    }
  }
  return(nodes)
}

.liv_plot <- function(res) {
  ggplot(res, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI"), alpha = 0.2) +
    geom_line(aes(y = est, color = "Estimate"), linewidth = 1) +
    geom_point(aes(y = est, color = "Estimate")) +
    geom_line(aes(y = ate, color = "True ATE"), linewidth = 1) +
    geom_point(aes(y = ate, color = "True ATE")) +
    scale_color_manual(
      name = "Lines",
      values = c("Estimate" = "blue", "True ATE" = "red")
    ) +
    scale_fill_manual(
      name = "Confidence Band",
      values = c("95% CI" = "blue")
    ) +
    labs(
      x = "Time", y = "Estimate",
      title = "Estimate with 95% Confidence Band"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

.liv_lapply_folds <- function(K, FUN) {
  run_fold <- function(k) {
    tryCatch(FUN(k), error = function(e) .liv_fold_error(k, e))
  }

  if (K <= 1) {
    return(.liv_stop_on_fold_errors(lapply(seq_len(K), run_fold)))
  }

  in_rstudio_or_knitr <- identical(Sys.getenv("RSTUDIO"), "1") ||
    nzchar(Sys.getenv("RSTUDIO_SESSION_ID")) ||
    isTRUE(getOption("knitr.in.progress"))
  use_parallel <- isTRUE(getOption("liv.parallel", !in_rstudio_or_knitr))

  if (!use_parallel) {
    return(.liv_stop_on_fold_errors(lapply(seq_len(K), run_fold)))
  }

  if (.Platform$OS.type == "windows") {
    warning("Parallel fold fitting uses forked workers and is not available on Windows; falling back to sequential fitting.")
    return(.liv_stop_on_fold_errors(lapply(seq_len(K), run_fold)))
  }

  results <- parallel::mclapply(seq_len(K), run_fold, mc.cores = K, mc.preschedule = FALSE)
  .liv_stop_on_fold_errors(results)
}

# Estimation -------------------------------------------------------------------

.liv_dml_overall <- function(mydat,
                             action_list,
                             K = 5,
                             H = 2,
                             superLearners = "mgcv",
                             print_state = FALSE) {
  TM <- max(mydat$t) - 1
  id_order <- sort(unique(mydat$ID))
  n <- length(id_order)

  ETA <- array(0, dim = c(n, K, TM + 2))
  terminal_Y <- mydat %>% filter(t == TM + 1) %>% arrange(ID) %>% pull(Y)
  ETA[, , TM + 2] <- matrix(rep(terminal_Y, K), n, K)
  folds <- .liv_make_id_folds(mydat, id_var = "ID", K = K)

  fold_results <- .liv_lapply_folds(K, function(k) {
    ETA_k <- matrix(0, nrow = n, ncol = TM + 2)
    ETA_k[, TM + 2] <- terminal_Y

    time <- TM
    while (time >= 0) {
      train_ids <- folds[[k]]$train_ids
      test_ids <- folds[[k]]$test_ids
      train <- match(train_ids, id_order)
      test <- match(test_ids, id_order)
      if (anyNA(train) || anyNA(test)) {
        stop("Fold IDs do not match the observed IDs.")
      }

      oY <- ETA_k[, time + 2]
      fY <- .liv_build_fY_history(mydat, time, H)
      oA <- mydat %>% filter(t == time) %>% arrange(ID) %>% pull(A)
      fA <- .liv_build_fA_history(mydat, time, H)
      oZ <- mydat %>% filter(t == time) %>% arrange(ID) %>% pull(Z)
      fZ <- .liv_build_fZ_history(mydat, time, H)

      fit_Y <- .liv_mgcv_fit(oY[train], fY[train, ])
      fit_A <- .liv_mgcv_fit(oA[train], fA[train, ], family = binomial())
      fit_Z <- .liv_mgcv_fit(oZ[train], fZ[train, ], family = binomial())

      Y11 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 1, A = 1))
      Y10 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 1, A = 0))
      Y01 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 0, A = 1))
      Y00 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 0, A = 0))
      pi1 <- .liv_mgcv_predict(fit_A, fA %>% mutate(Z = 1))
      pi0 <- .liv_mgcv_predict(fit_A, fA %>% mutate(Z = 0))
      omega <- pi0
      kappa <- pi1 - pi0
      f <- .liv_mgcv_predict(fit_Z, fZ)
      trt <- .liv_get_treatment(action_list, time, fY)

      mu1 <- trt * Y11 * pi1 + (trt - 1) * Y10 * (1 - pi1)
      mu0 <- trt * Y01 * pi0 + (trt - 1) * Y00 * (1 - pi0)
      gamma <- (mu1 - mu0) / kappa
      xi <- mu0
      weight <- oZ / f - (1 - oZ) / (1 - f)
      eta <- weight * ((oA - 1 + trt) / kappa * oY - xi / kappa) +
        (1 - weight * (oA - omega) / kappa) * gamma
      ETA_k[, time + 1] <- eta

      if (print_state) {
        message("fold ", k, ", time ", time)
      }
      time <- time - 1
    }
    ETA_k
  })

  for (k in seq_len(K)) {
    .liv_validate_fold_matrix(fold_results[[k]], n, TM + 2, "ETA", k)
    ETA[, k, ] <- fold_results[[k]]
  }

  list(
    estimand = "overall",
    ETA = ETA,
    folds = folds,
    id_order = id_order,
    input_data = mydat,
    Action = action_list,
    H = H,
    SupersuperLearners = superLearners
  )
}

.liv_dml_local <- function(mydat,
                           action_list,
                           K = 5,
                           H = 2,
                           superLearners = "mgcv",
                           print_state = FALSE) {
  TM <- max(mydat$t) - 1
  id_order <- sort(unique(mydat$ID))
  n <- length(id_order)

  ETA <- array(0, dim = c(n, K, TM + 2))
  terminal_Y <- mydat %>% filter(t == TM + 1) %>% arrange(ID) %>% pull(Y)
  ETA[, , TM + 2] <- matrix(rep(terminal_Y, K), n, K)

  ETA_den <- array(0, dim = c(n, K, TM + 2))
  ETA_den[, , TM + 2] <- 1
  folds <- .liv_make_id_folds(mydat, id_var = "ID", K = K)

  fold_results <- .liv_lapply_folds(K, function(k) {
    ETA_k <- matrix(0, nrow = n, ncol = TM + 2)
    ETA_den_k <- matrix(0, nrow = n, ncol = TM + 2)
    ETA_k[, TM + 2] <- terminal_Y
    ETA_den_k[, TM + 2] <- 1

    time <- TM
    while (time >= 0) {
      train_ids <- folds[[k]]$train_ids
      test_ids <- folds[[k]]$test_ids
      train <- match(train_ids, id_order)
      test <- match(test_ids, id_order)
      if (anyNA(train) || anyNA(test)) {
        stop("Fold IDs do not match the observed IDs.")
      }

      oY <- ETA_k[, time + 2]
      oden <- ETA_den_k[, time + 2]
      fY <- .liv_build_fY_history(mydat, time, H)
      oA <- mydat %>% filter(t == time) %>% arrange(ID) %>% pull(A)
      fA <- .liv_build_fA_history(mydat, time, H)
      oZ <- mydat %>% filter(t == time) %>% arrange(ID) %>% pull(Z)
      fZ <- .liv_build_fZ_history(mydat, time, H)

      fit_Y <- .liv_mgcv_fit(oY[train], fY[train, ])
      if (length(unique(oden[train])) == 1) {
        fit_den <- NULL
      } else {
        fit_den <- .liv_mgcv_fit(oden[train], fY[train, ], family = gaussian())
      }
      fit_A <- .liv_mgcv_fit(oA[train], fA[train, ], family = binomial())
      fit_Z <- .liv_mgcv_fit(oZ[train], fZ[train, ], family = binomial())

      Y11 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 1, A = 1))
      Y10 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 1, A = 0))
      Y01 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 0, A = 1))
      Y00 <- .liv_mgcv_predict(fit_Y, fY %>% mutate(Z = 0, A = 0))

      if (is.null(fit_den)) {
        den11 <- den10 <- den01 <- den00 <- rep(oden[train][1], n)
      } else {
        den11 <- .liv_mgcv_predict(fit_den, fY %>% mutate(Z = 1, A = 1))
        den10 <- .liv_mgcv_predict(fit_den, fY %>% mutate(Z = 1, A = 0))
        den01 <- .liv_mgcv_predict(fit_den, fY %>% mutate(Z = 0, A = 1))
        den00 <- .liv_mgcv_predict(fit_den, fY %>% mutate(Z = 0, A = 0))
      }

      pA1 <- .liv_mgcv_predict(fit_A, fA %>% mutate(Z = 1))
      pA0 <- .liv_mgcv_predict(fit_A, fA %>% mutate(Z = 0))
      f <- .liv_mgcv_predict(fit_Z, fZ)
      trt <- .liv_get_treatment(action_list, time, fY)

      mu1 <- trt * Y11 * pA1 + (trt - 1) * Y10 * (1 - pA1)
      mu0 <- trt * Y01 * pA0 + (trt - 1) * Y00 * (1 - pA0)
      pi1 <- trt * den11 * pA1 + (trt - 1) * den10 * (1 - pA1)
      pi0 <- trt * den01 * pA0 + (trt - 1) * den00 * (1 - pA0)

      weight <- oZ / f - (1 - oZ) / (1 - f)
      eta <- weight * (oA - 1 + trt) * oY +
        (1 - oZ / f) * mu1 -
        (1 - (1 - oZ) / (1 - f)) * mu0
      eta_den <- weight * (oA - 1 + trt) * oden +
        (1 - oZ / f) * pi1 -
        (1 - (1 - oZ) / (1 - f)) * pi0
      ETA_k[, time + 1] <- eta
      ETA_den_k[, time + 1] <- eta_den

      if (print_state) {
        message("fold ", k, ", time ", time)
      }
      time <- time - 1
    }
    list(ETA = ETA_k, ETA_den = ETA_den_k)
  })

  for (k in seq_len(K)) {
    .liv_validate_local_fold_result(fold_results[[k]], n, TM + 2, k)
    ETA[, k, ] <- fold_results[[k]]$ETA
    ETA_den[, k, ] <- fold_results[[k]]$ETA_den
  }

  list(
    estimand = "local",
    ETA = ETA,
    ETA_den = ETA_den,
    folds = folds,
    id_order = id_order,
    input_data = mydat,
    Action = action_list,
    H = H,
    SupersuperLearners = superLearners
  )
}

# Public simulation wrappers ---------------------------------------------------

sim_LIV_overall <- function(n, D) {
  .liv_sim_observed(n, D, drop_patterns = c("OPA", "deltaA", "pA0", "U"))
}

sim_LIV_local <- function(n, D) {
  .liv_sim_observed(n, D, drop_patterns = c("score", "R", "U"))
}

# Public DML wrappers ----------------------------------------------------------

DML_LIV <- function(mydat,
                    action_list,
                    estimand = c("overall", "local"),
                    K = 5,
                    H = 2,
                    superLearners = "mgcv",
                    print_state = FALSE) {
  estimand <- match.arg(estimand)
  if (estimand == "overall") {
    return(.liv_dml_overall(mydat, action_list, K, H, superLearners, print_state))
  }
  .liv_dml_local(mydat, action_list, K, H, superLearners, print_state)
}

# Public prediction wrappers ---------------------------------------------------

predict_LIV <- function(fit_LIV) {
  if (is.null(fit_LIV$estimand)) {
    stop("The fitted object does not contain an estimand field. Refit it with DML_LIV().", call. = FALSE)
  }
  .liv_predict(fit_LIV, local = fit_LIV$estimand == "local")
}

# Public interventional mean wrappers -----------------------------------------

ATE_overall <- function(D, TM, PROB, accuracy = 1e5) {
  PROB <- .liv_generate_action_from_vector(PROB)
  ate <- rep(0, TM + 2)

  for (s in 0:TM) {
    prob <- PROB[(s + 1):length(PROB)]
    Dact <- D + action(name = "do", nodes = prob)
    dat_do <- sim(Dact, n = accuracy, actions = "do", wide = TRUE)$do
    ate[s + 1] <- mean(dat_do[[paste("Y", TM + 1, sep = "_")]])
  }

  ate[TM + 2] <- mean(sim_LIV_overall(accuracy, D)$Y, na.rm = TRUE)
  return(ate)
}

ATE_local <- function(D, TM, PROB, accuracy = 1e5, state = 1) {
  PROB <- .liv_generate_action_from_vector(PROB)
  ate <- rep(0, TM + 2)

  for (s in 0:(TM + 1)) {
    prob <- PROB[(s + 1):length(PROB)]
    if (s < (TM + 1)) {
      Dact <- D + action(name = "do", nodes = prob)
      dat_do <- sim(Dact, n = accuracy, actions = "do", wide = TRUE)$do
      if (state == 1) {
        for (r in s:TM) {
          dat_do <- dat_do %>% filter(!!sym(paste("R", r, sep = "_")) %in% c("co", "de"))
        }
      }
    } else {
      dat_do <- sim(D, n = accuracy, wide = TRUE)
    }
    ate[s + 1] <- mean(dat_do[[paste("Y", TM + 1, sep = "_")]])
  }

  return(ate)
}

# Public plot wrappers ---------------------------------------------------------

plot_LIV_overall <- function(res, D, PROB) {
  TM <- nrow(res) - 2
  res$ate <- ATE_overall(D, TM, PROB)
  .liv_plot(res)
}

plot_LIV_local <- function(res, D, PROB, state = 1, TM = NULL) {
  if (is.null(TM)) {
    TM <- max(as.numeric(res$time)) - 1
  }

  res$ate <- ATE_local(D, TM, PROB, state = state)
  .liv_plot(res)
}
