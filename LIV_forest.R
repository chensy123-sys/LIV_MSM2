library(simcausal)
library(tibble)
library(dplyr)
library(randomForest)
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

.liv_forest_variable_types <- function(X, cts_num = 4) {
  continuous <- vapply(X, function(x) {
    is.numeric(x) && length(unique(x[!is.na(x)])) > cts_num
  }, logical(1))
  list(
    continuous = names(X)[continuous],
    discrete = names(X)[!continuous]
  )
}

.liv_forest_family_name <- function(family) {
  if (inherits(family, "family")) {
    return(family$family)
  }
  as.character(family)[1]
}

.liv_forest_numeric_response <- function(Y) {
  if (is.factor(Y)) {
    return(suppressWarnings(as.numeric(as.character(Y))))
  }
  as.numeric(Y)
}

.liv_forest_constant_fit <- function(Y, family) {
  value <- mean(.liv_forest_numeric_response(Y), na.rm = TRUE)
  if (!is.finite(value)) {
    stop("Cannot fit a constant learner because the response has no finite values.", call. = FALSE)
  }
  list(
    type = "constant",
    value = value,
    family = .liv_forest_family_name(family)
  )
}

.liv_forest_prepare_frame <- function(data,
                                      continuous_vars,
                                      discrete_vars,
                                      levels = NULL,
                                      medians = NULL) {
  data <- as.data.frame(data)[, c(continuous_vars, discrete_vars), drop = FALSE]

  if (is.null(medians)) {
    medians <- vapply(continuous_vars, function(var) {
      values <- as.numeric(data[[var]])
      value <- stats::median(values[is.finite(values)], na.rm = TRUE)
      if (is.finite(value)) value else 0
    }, numeric(1))
  }

  for (var in continuous_vars) {
    values <- as.numeric(data[[var]])
    values[!is.finite(values)] <- medians[[var]]
    data[[var]] <- values
  }

  if (is.null(levels)) {
    levels <- lapply(discrete_vars, function(var) {
      values <- unique(as.character(data[[var]]))
      values <- values[!is.na(values)]
      if (is.numeric(data[[var]]) || is.integer(data[[var]]) || is.logical(data[[var]])) {
        values <- unique(c(values, "0", "1"))
      }
      sort(unique(c(values, "__MISSING__", "__OTHER__")))
    })
    names(levels) <- discrete_vars
  }

  for (var in discrete_vars) {
    values <- as.character(data[[var]])
    values[is.na(values)] <- "__MISSING__"
    values[!values %in% levels[[var]]] <- "__OTHER__"
    data[[var]] <- factor(values, levels = levels[[var]])
  }

  list(data = data, levels = levels, medians = medians)
}

.liv_forest_matrix <- function(data,
                               continuous_vars,
                               discrete_vars,
                               levels = NULL,
                               medians = NULL,
                               feature_names = NULL) {
  prepared <- .liv_forest_prepare_frame(
    data,
    continuous_vars = continuous_vars,
    discrete_vars = discrete_vars,
    levels = levels,
    medians = medians
  )

  if (ncol(prepared$data) == 0) {
    matrix_data <- matrix(1, nrow = nrow(prepared$data), ncol = 1)
    colnames(matrix_data) <- ".intercept"
  } else {
    matrix_data <- stats::model.matrix(
      ~ . - 1,
      data = prepared$data,
      na.action = stats::na.pass
    )
    storage.mode(matrix_data) <- "double"
  }

  if (!is.null(feature_names)) {
    aligned <- matrix(0, nrow = nrow(matrix_data), ncol = length(feature_names))
    colnames(aligned) <- feature_names
    common <- intersect(colnames(matrix_data), feature_names)
    aligned[, common] <- matrix_data[, common, drop = FALSE]
    matrix_data <- aligned
  }

  list(
    matrix = matrix_data,
    levels = prepared$levels,
    medians = prepared$medians,
    feature_names = colnames(matrix_data)
  )
}

.liv_forest_fit <- function(Y, X, family, cts_num = 4) {
  X <- as.data.frame(X)
  Y <- .liv_forest_numeric_response(Y)
  keep <- is.finite(Y)
  X <- X[keep, , drop = FALSE]
  Y <- Y[keep]
  family_name <- .liv_forest_family_name(family)

  if (length(Y) == 0 || length(unique(Y)) == 1) {
    return(.liv_forest_constant_fit(Y, family_name))
  }
  if (family_name == "binomial" && !all(Y %in% c(0, 1))) {
    stop("Binomial random-forest responses must be coded as 0 and 1.", call. = FALSE)
  }

  types <- .liv_forest_variable_types(X, cts_num = cts_num)
  encoded <- .liv_forest_matrix(
    X,
    continuous_vars = types$continuous,
    discrete_vars = types$discrete
  )
  response <- if (family_name == "binomial") {
    factor(Y, levels = c(0, 1))
  } else {
    Y
  }

  params <- list(
    ntree = as.integer(getOption("liv.forest.ntree", 500L)),
    nodesize = as.integer(getOption(
      "liv.forest.nodesize",
      if (family_name == "binomial") 1L else 5L
    ))
  )
  mtry <- getOption("liv.forest.mtry", NULL)
  if (!is.null(mtry)) {
    params$mtry <- as.integer(mtry)
  }
  custom_params <- getOption("liv.forest.params", list())
  if (!is.list(custom_params)) {
    stop("Option 'liv.forest.params' must be a list.", call. = FALSE)
  }
  params <- utils::modifyList(params, custom_params)
  params$x <- NULL
  params$y <- NULL
  params$formula <- NULL
  params$data <- NULL

  object <- tryCatch(
    do.call(
      randomForest::randomForest,
      c(list(x = encoded$matrix, y = response), params)
    ),
    error = function(e) NULL
  )
  if (is.null(object)) {
    return(.liv_forest_constant_fit(Y, family_name))
  }

  list(
    type = "forest",
    object = object,
    family = family_name,
    fallback = mean(Y),
    continuous_vars = types$continuous,
    discrete_vars = types$discrete,
    discrete_levels = encoded$levels,
    continuous_medians = encoded$medians,
    feature_names = encoded$feature_names
  )
}

.liv_forest_predict <- function(fit, newdata) {
  newdata <- as.data.frame(newdata)

  if (fit$type == "constant") {
    pred <- rep(fit$value, nrow(newdata))
  } else {
    encoded <- .liv_forest_matrix(
      newdata,
      continuous_vars = fit$continuous_vars,
      discrete_vars = fit$discrete_vars,
      levels = fit$discrete_levels,
      medians = fit$continuous_medians,
      feature_names = fit$feature_names
    )
    pred <- tryCatch({
      if (fit$family == "binomial") {
        probabilities <- predict(fit$object, encoded$matrix, type = "prob")
        as.numeric(probabilities[, "1"])
      } else {
        as.numeric(predict(fit$object, encoded$matrix, type = "response"))
      }
    }, error = function(e) rep(fit$fallback, nrow(newdata)))
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
    SD <- sd(ETA_merge - ATE * ETA_den_merge) / abs(mean(ETA_den_merge)) / sqrt(n)
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
                             superLearners = "randomForest",
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

      fit_Y <- .liv_forest_fit(oY[train], fY[train, ], family = gaussian())
      fit_A <- .liv_forest_fit(oA[train], fA[train, ], family = binomial())
      fit_Z <- .liv_forest_fit(oZ[train], fZ[train, ], family = binomial())

      Y11 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 1, A = 1))
      Y10 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 1, A = 0))
      Y01 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 0, A = 1))
      Y00 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 0, A = 0))
      pi1 <- .liv_forest_predict(fit_A, fA %>% mutate(Z = 1))
      pi0 <- .liv_forest_predict(fit_A, fA %>% mutate(Z = 0))
      omega <- pi0
      kappa <- pi1 - pi0
      f <- .liv_forest_predict(fit_Z, fZ)
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
                           superLearners = "randomForest",
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

      fit_Y <- .liv_forest_fit(oY[train], fY[train, ], family = gaussian())
      if (length(unique(oden[train])) == 1) {
        fit_den <- NULL
      } else {
        fit_den <- .liv_forest_fit(oden[train], fY[train, ], family = gaussian())
      }
      fit_A <- .liv_forest_fit(oA[train], fA[train, ], family = binomial())
      fit_Z <- .liv_forest_fit(oZ[train], fZ[train, ], family = binomial())

      Y11 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 1, A = 1))
      Y10 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 1, A = 0))
      Y01 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 0, A = 1))
      Y00 <- .liv_forest_predict(fit_Y, fY %>% mutate(Z = 0, A = 0))

      if (is.null(fit_den)) {
        den11 <- den10 <- den01 <- den00 <- rep(oden[train][1], n)
      } else {
        den11 <- .liv_forest_predict(fit_den, fY %>% mutate(Z = 1, A = 1))
        den10 <- .liv_forest_predict(fit_den, fY %>% mutate(Z = 1, A = 0))
        den01 <- .liv_forest_predict(fit_den, fY %>% mutate(Z = 0, A = 1))
        den00 <- .liv_forest_predict(fit_den, fY %>% mutate(Z = 0, A = 0))
      }

      pA1 <- .liv_forest_predict(fit_A, fA %>% mutate(Z = 1))
      pA0 <- .liv_forest_predict(fit_A, fA %>% mutate(Z = 0))
      f <- .liv_forest_predict(fit_Z, fZ)
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
                    superLearners = "randomForest",
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
