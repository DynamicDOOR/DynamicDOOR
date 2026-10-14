#' fit_polr_probit: 
#' Fits a proportional-probit outcome model using MASS::polr.
#' If initialization fails, retries with zero slopes and evenly spaced
#' cumulative probabilities transformed to probit cutpoints.
#' Returns the fitted model with its Hessian.
fit_polr_probit <- function(formula, data) {
  # Avoiding a bad start due to outcome separation or extreme categories:
  # zero slopes + evenly spaced probit quantiles for cutpoints are always finite.
  tryCatch(
    MASS::polr(formula, data = data, method = "probit", Hess = TRUE),
    error = function(first_error) {
      init_error_msgs <- c(
        "initial value in 'vmmin' is not finite",
        "attempt to find suitable starting values failed"
      )
      if (!grepl(paste(init_error_msgs, collapse = "|"), conditionMessage(first_error))) {
        stop(first_error)
      }
      
      mf <- stats::model.frame(formula, data = data)
      y  <- stats::model.response(mf)
      X  <- stats::model.matrix(attr(mf, "terms"), mf)
      n_slopes <- ncol(X) - ("(Intercept)" %in% colnames(X))
      n_cuts   <- nlevels(y) - 1L
      stable_start <- c(rep(0, n_slopes), stats::qnorm(seq_len(n_cuts) / (n_cuts + 1L)))
      
      tryCatch(
        MASS::polr(formula, data = data, method = "probit", start = stable_start, Hess = TRUE),
        error = function(retry_error) {
          stop(sprintf(
            "polr failed. Default start: %s; stable-start retry: %s",
            conditionMessage(first_error), conditionMessage(retry_error)
          ), call. = FALSE)
        }
      )
    }
  )
}


#' fit_outcome_probs: Fits a proportional-probit, multinomial-logistic, or
#' classification-tree outcome model and predicts probabilities for newdata.
#' Aligns probability columns with outcome_levels and allows tree settings
#' to be modified through tree_control.
#' Returns the fitted model and predicted category probabilities.
fit_outcome_probs <- function(formula, data, newdata, outcome_levels=c('1','2', '3', '4'),
                              outcome_model='tree', tree_control = list()) {
  switch(outcome_model,
         polr = {
           data$Y <- ordered(data$Y, levels = outcome_levels)
           fit <- fit_polr_probit(formula, data)
           design <- stats::model.matrix(stats::delete.response(fit$terms), data,
                                         contrasts.arg = fit$contrasts, xlev = fit$xlevels)
           if (ncol(design) - ("(Intercept)" %in% colnames(design)) != length(fit$coefficients))
             stop("Proportional-probit model dropped coefficients because the design is rank-deficient; remove redundant covariates before fitting",
                  call. = FALSE)
           prob <- predict(fit, newdata = newdata, type = "probs")[, outcome_levels, drop = FALSE]
         },
         multinom = {
           data$Y <- factor(data$Y, levels = outcome_levels)
           fit <- nnet::multinom(formula, data = data, trace = FALSE, maxit = 500, MaxNWts = 10000, model = TRUE)
           prob <- predict(fit, newdata = newdata, type = "probs")[, outcome_levels, drop = FALSE]
         },
         tree = {
           data$Y <- factor(data$Y, levels = outcome_levels)
           control <- do.call(rpart::rpart.control, utils::modifyList(
             list(minsplit = 5, minbucket = 1, cp = 0, maxdepth = 30,
                  xval = 0, maxcompete = 0, maxsurrogate = 0), tree_control)) # override the default with tree_control
           fit <- rpart::rpart( formula, data = data, method = "class",
                                control = control )
           prob <- predict(fit, newdata = newdata, type = "prob")[, outcome_levels, drop = FALSE]
         }
  )
  list(fit = fit, prob = prob)
}


#' Training scores and Jacobian for the logistic propensity model
#'
#' Computes individual likelihood scores s_beta = X_beta * (Treatment - pi)
#' and their average Jacobian A_beta, using raw fitted logistic probabilities.
#' Only active_rows contribute; scores for all other participants are zero.
logistic_training_components <- function(fit, data, active_rows = seq_len(nrow(data))) {
  X_beta <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data,
                                contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  pi <- as.numeric(stats::predict(fit, newdata = data, type = "response"))
  active_rows <- as.integer(active_rows)
  A_beta <- -crossprod(X_beta[active_rows, , drop = FALSE], 
                       X_beta[active_rows, , drop = FALSE] * (pi[active_rows] * (1 - pi[active_rows]))) / length(active_rows)
  s_beta <- matrix(0, nrow(data), ncol(X_beta))
  s_beta[active_rows, ] <- X_beta[active_rows, , drop = FALSE] * (data$Treatment[active_rows] - pi[active_rows])
  list(A_beta = A_beta, s_beta = s_beta, X_beta = X_beta, fit = fit)
}


#' Training scores and Jacobian for the proportional-probit outcome model
#'
#' Computes the observed-category likelihood score g / q and its Jacobian
#' R / q - g %*% t(g) / q^2, where g and R are the probability gradient
#' and Hessian. Parameters are ordered as slopes followed by cutpoints.
#' Only participants selected by arm_indicator contribute; all other
#' score rows are zero.
polr_training_components <- function(fit, obs_data, arm_indicator, outcome_numeric, probabilities,
                                     expansion_n = nrow(obs_data)) {
  dq <- polr_probability_derivatives(fit, obs_data, probabilities)
  X_all <- stats::model.matrix(stats::delete.response(fit$terms), obs_data, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  if ("(Intercept)" %in% colnames(X_all))
    X_all <- X_all[, colnames(X_all) != "(Intercept)", drop = FALSE]
  
  n <- nrow(X_all); K <- dim(dq)[2]; p <- ncol(X_all); d <- dim(dq)[3]
  eta_all <- as.vector(X_all %*% fit$coefficients)
  cut <- c(-Inf, as.numeric(fit$zeta), Inf)
  slope_cols <- seq_len(p)
  s_theta <- matrix(0, n, d)
  A_theta <- matrix(0, d, d)
  
  for (i in which(arm_indicator == 1)) {
    k <- as.integer(outcome_numeric[i])
    x <- X_all[i, , drop = TRUE]
    c_lo <- cut[k] - eta_all[i]
    c_hi <- cut[k + 1L] - eta_all[i]
    q_i <- max(stats::pnorm(c_hi) - stats::pnorm(c_lo), .Machine$double.eps)
    
    # Probability gradient and likelihood score for the observed category
    g <- dq[i, k, ]
    s_theta[i, ] <- g / q_i
    
    # Probability Hessian for the observed category
    c_lo_f_lo <- if (is.finite(c_lo)) c_lo * stats::dnorm(c_lo) else 0
    c_hi_f_hi <- if (is.finite(c_hi)) c_hi * stats::dnorm(c_hi) else 0
    R <- matrix(0, d, d)
    
    if (p > 0L)
      R[slope_cols, slope_cols] <- (c_lo_f_lo - c_hi_f_hi) * tcrossprod(x)
    
    if (k < K) {
      upper_col <- p + k
      R[slope_cols, upper_col] <- c_hi_f_hi * x
      R[upper_col, slope_cols] <- c_hi_f_hi * x
      R[upper_col, upper_col] <- -c_hi_f_hi
    }
    if (k > 1L) {
      lower_col <- p + k - 1L
      R[slope_cols, lower_col] <- -c_lo_f_lo * x
      R[lower_col, slope_cols] <- -c_lo_f_lo * x
      R[lower_col, lower_col] <- c_lo_f_lo
    }
    
    # Likelihood-score Jacobian
    A_theta <- A_theta + R / q_i - tcrossprod(g) / q_i^2
  }
  
  list(A_theta = A_theta / expansion_n, s_theta = s_theta)
}

#' Training scores and Jacobian for the multinomial-logit outcome model
#'
#' Computes the likelihood score X * (I(Y = k) - q_k) for each nonreference
#' category and obtains its Jacobian from multinom_probability_derivatives().
#' Parameter blocks follow the nonreference categories in fit$lev, with
#' design-matrix coefficients ordered within each block.
#' Only participants selected by arm_indicator contribute; all other
#' score rows are zero.
multinom_training_components <- function(fit, obs_data, arm_indicator, outcome_numeric, probabilities,
                                         expansion_n = nrow(obs_data)) {
  fitted_categories <- as.integer(as.character(fit$lev))
  nonref <- fitted_categories[-1]
  X_all <- stats::model.matrix(stats::delete.response(fit$terms), obs_data, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  n <- nrow(obs_data); p <- ncol(X_all); d <- p * length(nonref)
  
  arm_rows <- which(arm_indicator == 1)
  X_train <- X_all[arm_rows, , drop = FALSE]
  y_train <- as.integer(outcome_numeric[arm_rows])
  n_arm <- length(arm_rows)
  
  # probabilities: raw fitted probabilities for all participants, columns 1,...,K
  q_train <- probabilities[arm_rows, , drop = FALSE]
  dq_train <- multinom_probability_derivatives(fit, obs_data[arm_rows, , drop = FALSE], q_train)
  
  s_theta <- matrix(0, n, d)
  A_theta <- matrix(0, d, d)
  
  for (j in seq_along(nonref)) {
    category <- nonref[j]
    cols <- (j - 1L) * p + seq_len(p)
    
    # Likelihood scores, with zeros outside the selected training arm
    s_theta[arm_rows, cols] <- X_train * (as.numeric(y_train == category) - q_train[, category])
    
    # Likelihood-score Jacobian
    dq_category <- matrix(dq_train[, category, ], nrow = n_arm, ncol = d)
    A_theta[cols, ] <- -crossprod(X_train, dq_category)
  }
  
  list(A_theta = A_theta / expansion_n, s_theta = s_theta)
}


#' fit_nuisance_ind: Fits a logistic propensity model on the full sample and
#' separate outcome models for each treatment arm and time point.
#' Returns fitted models, propensity components, and category probabilities
#' arranged as time point by subject by category.
#' With varc = TRUE, prepare nuisance components/Jacobian blocks. 
#' For all methods, propensity components are prepared
#' for polr and multinom, outcome components are prepared
#' fixed_propensity optionally fixes P(A = 1); no propensity model or scores are fitted.
fit_nuisance_ind <- function(obs_data, trt_predictors, outcome_predictors,
                             K = 4, n_timepoints = length(outcome_columns),
                             outcome_model = c("polr", "tree", "multinom"),
                             varc = TRUE, outcome_columns, tree_control = list(), fixed_propensity = NULL) {
  
  outcome_model <- match.arg(outcome_model)
  stopifnot(is.logical(varc), length(varc) == 1L, !is.na(varc), n_timepoints == length(outcome_columns))

  n <- nrow(obs_data)
  q1 <- q0 <- array(0, dim = c(n_timepoints, n, K))
  outcome_levels <- as.character(seq_len(K))
  
  # Estimated or fixed working propensity
  fit_pi <- NULL
  if (is.null(fixed_propensity)) {
    formula_pi <- reformulate(trt_predictors, response = "Treatment")
    fit_pi <- glm(formula_pi, family = binomial(), data = obs_data)
    pi <- predict(fit_pi, newdata = obs_data, type = "response")
  } else {
    stopifnot(is.numeric(fixed_propensity), length(fixed_propensity) == 1L,
              is.finite(fixed_propensity), fixed_propensity > 0, fixed_propensity < 1)
    pi <- rep(fixed_propensity, n)
  }
  
  # Outcome data: one row per participant and visit
  data_sim_long <- obs_data[rep(seq_len(n), times = n_timepoints), , drop = FALSE]
  data_sim_long$id <- rep(seq_len(n), times = n_timepoints)
  data_sim_long$time <- rep(seq_len(n_timepoints), each = n)
  data_sim_long$Y <- factor(unlist(lapply(obs_data[outcome_columns], as.character)), levels = outcome_levels)
  
  formula_outcome <- reformulate(outcome_predictors, response = "Y")
  fit_trt <- fit_control <- vector("list", n_timepoints)
  
  for (t in seq_len(n_timepoints)) {
    dat_t <- data_sim_long[data_sim_long$time == t, , drop = FALSE]
    dat_t_trt <- dat_t[dat_t$Treatment == 1, , drop = FALSE]
    dat_t_control <- dat_t[dat_t$Treatment == 0, , drop = FALSE]
    
    fitted_trt <- fit_outcome_probs(formula_outcome, dat_t_trt, dat_t, outcome_levels, outcome_model, tree_control = tree_control)
    fitted_control <- fit_outcome_probs(formula_outcome, dat_t_control, dat_t, outcome_levels, outcome_model, tree_control = tree_control)
    
    fit_trt[[t]] <- fitted_trt$fit
    fit_control[[t]] <- fitted_control$fit
    q1[t, , ] <- fitted_trt$prob
    q0[t, , ] <- fitted_control$prob
  }
  
  
  nuisance_components <- NULL
  if (varc) {
    nuisance_components <- tryCatch({
      components <- list(beta = if (is.null(fixed_propensity)) logistic_training_components(fit_pi, obs_data) else NULL,
                         theta1 = NULL, theta0 = NULL)
      
      if (outcome_model %in% c("polr", "multinom")) {
        outcome_training_components <- switch(outcome_model, polr = polr_training_components, multinom = multinom_training_components)
        components$theta1 <- vector("list", n_timepoints)
        components$theta0 <- vector("list", n_timepoints)
        Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
        
        for (t in seq_len(n_timepoints)) {
          components$theta1[[t]] <- outcome_training_components(fit_trt[[t]], obs_data, as.numeric(obs_data$Treatment == 1), Y_mat[, t], q1[t, , ])
          components$theta0[[t]] <- outcome_training_components(fit_control[[t]], obs_data, as.numeric(obs_data$Treatment == 0), Y_mat[, t], q0[t, , ])
        }
      }
      components
    }, analytic_correction_error = function(e) list(error = conditionMessage(e)))
  }
  
  # Record clipping rates before clipping
  clip <- probability_clip()
  clipping <- c(propensity = mean(pi < clip | pi > 1 - clip),
                outcome = mean(c(q1 < clip | q1 > 1 - clip, q0 < clip | q0 > 1 - clip)))
  
  pi <- clip01(pi, clip)
  q1 <- clip01(q1, clip)
  q0 <- clip01(q0, clip)
  
  for (t in seq_len(n_timepoints)) {
    q1[t, , ] <- q1[t, , ] / rowSums(q1[t, , ])
    q0[t, , ] <- q0[t, , ] / rowSums(q0[t, , ])
  }
  
  
  list(pi = pi, q1 = q1, q0 = q0, fit_pi = fit_pi, fit_trt = fit_trt, fit_control = fit_control,
       nuisance_components = nuisance_components, varc = varc, outcome_model = outcome_model,
       outcome_columns = outcome_columns, clipping = clipping, fixed_propensity = fixed_propensity)
}

#' fit_nuisance_ind_CF: Fits nuisance models using cross-fitting with supplied
#' folds or randomly generated folds. For each fold, fits a logistic propensity
#' model and separate outcome models by treatment arm and time point on the
#' training sample, then predicts probabilities for the validation sample.
#' Returns assembled validation predictions, fold assignments, and fitted models.
#' With varc = TRUE, prepare nuisance components/Jacobian blocks. 
#' For all methods, fold-specific propensity components are prepared
#' for polr and multinom, fold-specific outcome components are prepared
#' fixed_propensity optionally fixes P(A = 1) in every fold, omitting propensity score blocks.
fit_nuisance_ind_CF <- function(obs_data, V, trt_predictors, outcome_predictors,
                                K = 4, n_timepoints = length(outcome_columns), seed = 0,
                                outcome_model = c("polr", "tree", "multinom"),
                                folds = NULL, varc = TRUE, outcome_columns, tree_control = list(), fixed_propensity = NULL) {
  
  outcome_model <- match.arg(outcome_model)
  stopifnot(is.logical(varc), length(varc) == 1L, !is.na(varc), n_timepoints == length(outcome_columns))
  
  n <- nrow(obs_data)
  if (is.null(folds)) folds <- kfold_indices_cs(n, V = V, seed = seed)
  if (length(folds) != n || anyNA(folds) || V < 2L || length(unique(folds)) != V)
    stop("folds must assign every participant to one of V nonempty folds, with V >= 2")
  
  if (!is.null(fixed_propensity))
    stopifnot(is.numeric(fixed_propensity), length(fixed_propensity) == 1L,
              is.finite(fixed_propensity), fixed_propensity > 0, fixed_propensity < 1)
  pi <- if (is.null(fixed_propensity)) numeric(n) else rep(fixed_propensity, n)
  q1 <- q0 <- array(0, dim = c(n_timepoints, n, K))
  outcome_levels <- as.character(seq_len(K))
  fold_ids <- sort(unique(folds))
  
  fit_pi_folds <- setNames(vector("list", length(fold_ids)), fold_ids)
  fit_trt_folds <- setNames(vector("list", length(fold_ids)), fold_ids)
  fit_control_folds <- setNames(vector("list", length(fold_ids)), fold_ids)
  
  # Outcome data: one row per participant and visit
  data_sim_long <- obs_data[rep(seq_len(n), times = n_timepoints), , drop = FALSE]
  data_sim_long$id <- rep(seq_len(n), times = n_timepoints)
  data_sim_long$time <- rep(seq_len(n_timepoints), each = n)
  data_sim_long$Y <- factor(unlist(lapply(obs_data[outcome_columns], as.character)), levels = outcome_levels)
  if (is.null(fixed_propensity)) formula_pi <- reformulate(trt_predictors, response = "Treatment")
  formula_outcome <- reformulate(outcome_predictors, response = "Y")
  
  for (v in fold_ids) {
    key <- as.character(v)
    in_v <- which(folds == v)
    out_v <- which(folds != v)
    
    # Train outside I_v; predict on I_v
    if (is.null(fixed_propensity)) {
      fit_pi <- glm(formula_pi, family = binomial(), data = obs_data[out_v, , drop = FALSE])
      pi[in_v] <- predict(fit_pi, newdata = obs_data[in_v, , drop = FALSE], type = "response")
      fit_pi_folds[[key]] <- fit_pi
    }
    fit_trt_folds[[key]] <- vector("list", n_timepoints)
    fit_control_folds[[key]] <- vector("list", n_timepoints)
    
    for (t in seq_len(n_timepoints)) {
      train_t <- data_sim_long[data_sim_long$time == t & data_sim_long$id %in% out_v, , drop = FALSE]
      dat_t_trt <- train_t[train_t$Treatment == 1, , drop = FALSE]
      dat_t_control <- train_t[train_t$Treatment == 0, , drop = FALSE]
      eval_data <- obs_data[in_v, , drop = FALSE]
      
      fitted_trt <- fit_outcome_probs(formula_outcome, dat_t_trt, eval_data, outcome_levels, outcome_model, tree_control = tree_control)
      fitted_control <- fit_outcome_probs(formula_outcome, dat_t_control, eval_data, outcome_levels, outcome_model, tree_control = tree_control)
      
      fit_trt_folds[[key]][[t]] <- fitted_trt$fit
      fit_control_folds[[key]][[t]] <- fitted_control$fit
      q1[t, in_v, ] <- fitted_trt$prob
      q0[t, in_v, ] <- fitted_control$prob
    }
  }
  
  # Optional training components: propensity for every model, outcome for parametric models
  # components use raw model probabilities and are zero outside the training sample.
  # The correction calculation applies the n_v/(n-n_v) training weights.
  
  nuisance_components_folds <- NULL
  if (varc) {
    nuisance_components_folds <- tryCatch({
      components_folds <- setNames(vector("list", length(fold_ids)), fold_ids)
      
      if (outcome_model %in% c("polr", "multinom")) {
        outcome_training_components <- switch(outcome_model, polr = polr_training_components, multinom = multinom_training_components)
        Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
      }
      
      for (v in fold_ids) {
        key <- as.character(v)
        out_v <- which(folds != v)
        n_out <- length(out_v)
        
        components_folds[[key]] <- list(
          beta = if (is.null(fixed_propensity)) logistic_training_components(fit_pi_folds[[key]], obs_data, active_rows = out_v) else NULL,
          theta1 = NULL, theta0 = NULL)
        
        if (outcome_model %in% c("polr", "multinom")) {
          components_folds[[key]]$theta1 <- vector("list", n_timepoints)
          components_folds[[key]]$theta0 <- vector("list", n_timepoints)
          arm1 <- as.numeric(seq_len(n) %in% out_v & obs_data$Treatment == 1)
          arm0 <- as.numeric(seq_len(n) %in% out_v & obs_data$Treatment == 0)
          
          for (t in seq_len(n_timepoints)) {
            fit1 <- fit_trt_folds[[key]][[t]]
            fit0 <- fit_control_folds[[key]][[t]]
            q1_all <- predict(fit1, newdata = obs_data, type = "probs")[, outcome_levels, drop = FALSE]
            q0_all <- predict(fit0, newdata = obs_data, type = "probs")[, outcome_levels, drop = FALSE]
            
            components_folds[[key]]$theta1[[t]] <- outcome_training_components(fit1, obs_data, arm1, Y_mat[, t], q1_all, expansion_n = n_out)
            components_folds[[key]]$theta0[[t]] <- outcome_training_components(fit0, obs_data, arm0, Y_mat[, t], q0_all, expansion_n = n_out)
          }
        }
      }
      components_folds
    }, analytic_correction_error = function(e) list(error = conditionMessage(e)))
  }
  
  
  # Record clipping rates before clipping
  clip <- probability_clip()
  clipping <- c(propensity = mean(pi < clip | pi > 1 - clip),
                outcome = mean(c(q1 < clip | q1 > 1 - clip, q0 < clip | q0 > 1 - clip)))
  
  pi <- clip01(pi, clip)
  q1 <- clip01(q1, clip)
  q0 <- clip01(q0, clip)
  
  for (t in seq_len(n_timepoints)) {
    q1[t, , ] <- q1[t, , ] / rowSums(q1[t, , ])
    q0[t, , ] <- q0[t, , ] / rowSums(q0[t, , ])
  }
  
  
  
  list(pi = pi, q1 = q1, q0 = q0, folds = folds, fit_pi_folds = fit_pi_folds,
       fit_trt_folds = fit_trt_folds, fit_control_folds = fit_control_folds,
       nuisance_components_folds = nuisance_components_folds, varc = varc,
       outcome_model = outcome_model, outcome_columns = outcome_columns, clipping = clipping,
       fixed_propensity = fixed_propensity)
}
