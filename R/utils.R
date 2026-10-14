#' probability_clip: Returns the clipping threshold for fitted probabilities.
probability_clip <- function() getOption("dynamic_door_probability_clip", 1e-5)

#' polr_probability_derivatives: Computes first-order derivatives of initial
#' category probabilities with respect to the slopes and cutpoints of a fitted
#' proportional-probit model. Returns an array indexed by observation,
#' category, and parameter.
polr_probability_derivatives <- function(fit, newdata, probabilities) {
  X_eval <- stats::model.matrix(stats::delete.response(fit$terms), newdata, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  if ("(Intercept)" %in% colnames(X_eval))
    X_eval <- X_eval[, colnames(X_eval) != "(Intercept)", drop = FALSE]
  
  n_eval <- nrow(X_eval); K <- ncol(probabilities); p <- ncol(X_eval); d <- p + K - 1L
  zeta <- as.numeric(fit$zeta)
  eta_eval <- as.vector(X_eval %*% fit$coefficients)
  density_at_cut <- matrix(vapply(zeta, function(z) stats::dnorm(z - eta_eval), numeric(n_eval)), nrow = n_eval, ncol = K - 1L)
  density_ext <- cbind(0, density_at_cut, 0)
  dq_theta <- array(0, dim = c(n_eval, K, d))
  
  for (k in seq_len(K)) {
    if (p > 0L) dq_theta[, k, seq_len(p)] <- X_eval * (density_ext[, k] - density_ext[, k + 1L])
    if (k < K) dq_theta[, k, p + k] <- density_at_cut[, k]
    if (k > 1L) dq_theta[, k, p + k - 1L] <- -density_at_cut[, k - 1L]
  }
  
  dq_theta
}


#' multinom_probability_derivatives: Computes first-order derivatives of initial
#' category probabilities with respect to the coefficients for nonreference
#' categories of a fitted multinomial-logistic model. Returns an array indexed
#' by observation, category, and parameter.
multinom_probability_derivatives <- function(fit, newdata, probabilities) {
  K <- ncol(probabilities)
  fitted_categories <- as.integer(as.character(fit$lev))
  X_eval <- stats::model.matrix( stats::delete.response(fit$terms), newdata, contrasts.arg = fit$contrasts, xlev = fit$xlevels)
  n_eval <- nrow(X_eval)
  p <- ncol(X_eval)
  nonref <- fitted_categories[-1]
  d <- p * length(nonref)
  dq_theta <- array(0, dim = c(n_eval, K, d))
  
  for (k in seq_len(K)) {
    for (j in seq_along(nonref)) {
      category <- nonref[j]
      cols <- ((j - 1) * p + 1):(j * p)
      derivative_multiplier <- probabilities[, k] * (as.numeric(k == category) - probabilities[, category])
      dq_theta[, k, cols] <- X_eval * derivative_multiplier
    }
  }
  dq_theta
}

#' long_to_wide:
#' Convert one row per subject and visit to one row per subject. Supply the
#' column names for id, time, outcome, and baseline_columns (including treatment).
#' Baseline columns must be constant within subjects; duplicate subject/visit
#' rows are rejected. Subjects retain their order of first appearance.
#' times selects and orders visits; outcome_columns gives their output names,
#' defaulting to outcome_time. Outcome classes and factor levels are preserved.
#' Missing visits become NA, which fit_dynamic_door() does not accept.
long_to_wide <- function(data, id, time, outcome, baseline_columns,
                         times = sort(unique(data[[time]])),
                         outcome_columns = paste0(outcome, "_", times)) {
  data <- as.data.frame(data)
  data <- data[, unique(c(id, time, outcome, baseline_columns)), drop = FALSE]
  if (anyNA(data[c(id, time)]) || anyDuplicated(data[c(id, time)]))
    stop("Each subject/visit pair must be unique and have nonmissing id and time")
  if (!length(times) || anyNA(times) || anyDuplicated(as.character(times)) ||
      length(outcome_columns) != length(times) || anyNA(outcome_columns) ||
      anyDuplicated(c(id, baseline_columns, outcome_columns)))
    stop("Supply distinct times and one distinct outcome column name per visit, separate from id and baseline columns")
  
  wide <- unique(data[c(id, baseline_columns)])
  if (anyDuplicated(wide[[id]]))
    stop("baseline_columns must be constant within each subject")
  for (t in seq_along(times)) {
    rows <- which(as.character(data[[time]]) == as.character(times[t]))
    wide[[outcome_columns[t]]] <- data[[outcome]][rows[match(wide[[id]], data[[id]][rows])]]
  }
  rownames(wide) <- NULL
  wide
}

#' kfold_indices_cs:
#' Randomly assign subjects to approximately equal-sized cross-fitting folds.
#' Uses the specified seed for reproducibility.
kfold_indices_cs <- function(n, V = 5, seed = 1) {
  set.seed(seed); sample(rep(1:V, length.out = n))
}


#' complete_covariates:
#' Check that each predictor has at least two observed values and that
#' factor predictors retain all defined levels.
#' Returns TRUE/FALSE, or issue descriptions when details = TRUE.
complete_covariates <- function(data, predictors, details = FALSE) {
  issues <- character()
  for (variable in predictors) {
    observed <- unique(as.character(data[[variable]][!is.na(data[[variable]])]))
    if (length(observed) < 2L)
      issues <- c(issues, sprintf("covariate '%s' has fewer than two distinct values", variable))
    if (is.factor(data[[variable]])) {
      missing_levels <- setdiff(levels(data[[variable]]), observed)
      if (length(missing_levels))
        issues <- c(issues, sprintf("covariate '%s' is missing levels: %s", variable, paste(missing_levels, collapse = ", ")))
    }
    if (!details && length(issues)) return(FALSE)
  }
  if (details) issues else TRUE
}

#' complete_outcomes:
#' Check that both treatment arms are present and that every outcome column
#' contains all category codes 1,...,K within each arm, with no unexpected codes.
#' Returns TRUE/FALSE, or issue descriptions when details = TRUE.
complete_outcomes <- function(data, K = 4, outcome_columns,
                              details = FALSE, outcome_levels = as.character(seq_len(K))) {
  issues <- character()
  for (arm in c(0, 1)) {
    arm_data <- data[data$Treatment == arm, , drop = FALSE]
    if (!nrow(arm_data)) {
      if (!details) return(FALSE)
      issues <- c(issues, sprintf("treatment arm %s is absent", arm))
      next
    }
    for (outcome in outcome_columns) {
      observed <- unique(as.character(arm_data[[outcome]]))
      missing_levels <- setdiff(as.character(seq_len(K)), observed)
      unexpected_levels <- setdiff(observed, as.character(seq_len(K)))
      if (!details && (length(missing_levels) || length(unexpected_levels))) return(FALSE)
      if (length(missing_levels))
        issues <- c(issues, sprintf("outcome '%s', treatment arm %s, is missing categories: %s",
                                    outcome, arm, paste(outcome_levels[as.integer(missing_levels)], collapse = ", ")))
      if (length(unexpected_levels))
        issues <- c(issues, sprintf("outcome '%s', treatment arm %s, contains invalid category codes: %s",
                                    outcome, arm, paste(unexpected_levels, collapse = ", ")))
    }
  }
  if (details) issues else TRUE
}

#' nuisance_support_complete
#' Check covariate and outcome support in the full sample or each cross-fitting
#' training sample. Outcome predictors are checked separately within each arm.
#' Returns TRUE/FALSE, or issue descriptions when details = TRUE.
nuisance_support_complete <- function(dat, trt_predictors, outcome_predictors, K = 4, V = 5,
                                      seed = 1, folds = NULL, cf = TRUE, outcome_columns,
                                      details = FALSE, outcome_levels = as.character(seq_len(K))) {
  issues <- character()
  samples <- list("Full sample" = seq_len(nrow(dat)))
  if (cf) {
    if (is.null(folds)) folds <- kfold_indices_cs(nrow(dat), V, seed)
    samples <- list()
    for (f in sort(unique(folds)))
      samples[[paste0("Training sample excluding fold ", f)]] <- which(folds != f)
  }
  for (scope in names(samples)) {
    sample_data <- dat[samples[[scope]], , drop = FALSE]
    checks <- list(
      "Propensity model" = complete_covariates(sample_data, trt_predictors, details = details),
      "Outcomes" = complete_outcomes(sample_data, K, outcome_columns, details = details, outcome_levels = outcome_levels)
    )
    for (arm in c(0, 1)) {
      arm_data <- sample_data[sample_data$Treatment == arm, , drop = FALSE]
      if (nrow(arm_data))
        checks[[paste0("Outcome model, treatment arm ", arm)]] <- complete_covariates(arm_data, outcome_predictors, details = details)
    }
    if (details) {
      for (model in names(checks)) {
        if (length(checks[[model]]))
          issues <- c(issues, paste0(scope, ": ", model, ": ", checks[[model]]))
      }
    } else if (!all(unlist(checks))) {
      return(FALSE)
    }
  }
  if (details) issues else TRUE
}

#' clip01
#' Restrict probabilities to [eps, 1 - eps] to avoid values at zero or one.
clip01 <- function(x, eps = 1e-5) {
  pmin(pmax(x, eps), 1 - eps)
}


#' Solves X %*% z = b, using a generalized inverse when needed.
#' Uses solve() for well-conditioned matrices. For singular/ill-conditioned
#' score Jacobians or a failed solve(), warns and falls back to MASS::ginv().
#' Nonfinite inputs still raise an analytic_correction_error.
checked_solve <- function(x, b) {
  if (any(!is.finite(x)) || any(!is.finite(b)))
    stop(structure(list(message = "Analytic correction unavailable: nonfinite matrix or right-hand side", call = NULL),
                   class = c("analytic_correction_error", "error", "condition")))
  reciprocal_condition <- tryCatch(rcond(x), error = function(e) NA_real_)
  if (is.finite(reciprocal_condition) && reciprocal_condition >= 1e-12) {
    solution <- tryCatch(solve(x, b), error = function(e) NULL)
    if (!is.null(solution) && all(is.finite(solution))) return(solution)
  }
  warning("Singular/ill-conditioned score Jacobian or failed solve(); using MASS::ginv() for the analytic correction",
          call. = FALSE)
  solution <- MASS::ginv(x) %*% b
  if (is.null(dim(b))) solution <- drop(solution)
  solution
}
