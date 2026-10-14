#' AIPW outcome-parameter derivatives on the evaluation sample
#'
#' Computes G_at_theta as the average of h_weight times the category-probability
#' derivative. For AIPW, h_weight = 1 - I(A = a) / pi_a.
#' Uses the full sample without cross-fitting and the validation fold otherwise.
#'
aipw_val_components <- function(fit, newdata, probabilities, h_weight, outcome_model = c("polr", "multinom")) {
  outcome_model <- match.arg(outcome_model)
  dq <- switch(outcome_model,
               polr = polr_probability_derivatives(fit, newdata, probabilities),
               multinom = multinom_probability_derivatives(fit, newdata, probabilities)
  )
  G_at_theta <- apply(sweep(dq, 1, h_weight, FUN = "*"), c(2, 3), mean)
  list(G_at_theta = G_at_theta)
}


#' Full-sample AIPW estimation of the DOOR trajectory
#'
#' Estimates each arm's marginal category probabilities by averaging augmented
#' contributions, then calculates psi[t] = p1[t, ] %*% W %*% p0[t, ].
#' Returns ordinary participant influence values and, when nuis$varc is TRUE,
#' corrected influence values using the selected correction.
#' The proposition correction uses propensity and outcome scores, Jacobians,
#' and prediction derivatives; it requires polr or multinom outcome models.
#' With a fixed working propensity, only outcome parameters enter this correction.
#' The corollary correction uses only propensity components, but assumes 
#' correctly specified propensity score.
#' 
AIPW <- function(obs_data, nuis, K = dim(nuis$q1)[3], n_timepoints = dim(nuis$q1)[1],
                 outcome_columns = nuis$outcome_columns, correction = c("proposition", "corollary")) {
  start <- Sys.time()
  correction <- match.arg(correction)
  outcome_columns <- outcome_columns[seq_len(n_timepoints)]
  
  n <- nrow(obs_data); A <- obs_data$Treatment; pi <- nuis$pi
  Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
  W <- outer(seq_len(K), seq_len(K), function(k, l) (k > l) + 0.5 * (k == l))
  p1 <- p0 <- matrix(0, n_timepoints, K)
  psi <- numeric(n_timepoints)
  phi <- matrix(0, n, n_timepoints)
  
  # Augmented category contributions, pooled margins, and ordinary influence values
  for (t in seq_len(n_timepoints)) {
    q1 <- matrix(nuis$q1[t, , ], n, K)
    q0 <- matrix(nuis$q0[t, , ], n, K)
    indicators <- outer(Y_mat[, t], seq_len(K), "==")
    m1 <- q1 + (A / pi) * (indicators - q1)
    m0 <- q0 + ((1 - A) / (1 - pi)) * (indicators - q0)
    
    p1[t, ] <- colMeans(m1)
    p0[t, ] <- colMeans(m0)
    psi[t] <- as.numeric(p1[t, ] %*% W %*% p0[t, ])
    phi[, t] <- as.vector(sweep(m1, 2, p1[t, ], FUN = "-") %*% (W %*% p0[t, ])) +
      as.vector(sweep(m0, 2, p0[t, ], FUN = "-") %*% (t(W) %*% p1[t, ]))
  }
  
  base <- list(psi = psi, phi = phi, p1 = p1, p0 = p0, time = Sys.time() - start)
  if (!nuis$varc) return(base)
  
  beta_scores <- nuis$nuisance_components$beta
  
  # Corollary: propensity-score correction, using I_beta = -A_beta
  if (correction == "corollary") {
    if (is.null(beta_scores)) stop("The corollary requires fitted propensity components; use the proposition for a fixed working propensity.")
    s_beta <- beta_scores$s_beta
    s_beta <- sweep(s_beta, 2, colMeans(s_beta), FUN = "-")
    phi_centered <- sweep(phi, 2, colMeans(phi), FUN = "-")
    C <- crossprod(phi_centered, s_beta) / n
    I_beta <- -beta_scores$A_beta
    
    base$phi_corrected <- phi - t(C %*% checked_solve(I_beta, t(s_beta)))
    base$time <- Sys.time() - start
    return(base)
  }
  
  if (!nuis$outcome_model %in% c("polr", "multinom"))
    stop("The proposition correction requires a polr or multinom outcome model.")
  
  # Proposition: complete parametric stacked-score correction
  phi_corrected <- phi
  # Smooth derivative convention; active clipping is qualified in the appendix.
  if (!is.null(beta_scores)) dpi <- beta_scores$X_beta * (pi * (1 - pi))
  h1 <- 1 - A / pi
  h0 <- 1 - (1 - A) / (1 - pi)
  
  for (t in seq_len(n_timepoints)) {
    q1 <- matrix(nuis$q1[t, , ], n, K)
    q0 <- matrix(nuis$q0[t, , ], n, K)
    indicators <- outer(Y_mat[, t], seq_len(K), "==")
    
    theta1 <- aipw_val_components(nuis$fit_trt[[t]], obs_data, q1, h1, outcome_model = nuis$outcome_model)
    theta0 <- aipw_val_components(nuis$fit_control[[t]], obs_data, q0, h0, outcome_model = nuis$outcome_model)
    
    theta1_scores <- nuis$nuisance_components$theta1[[t]]
    theta0_scores <- nuis$nuisance_components$theta0[[t]]
    correction_p1 <- -t(theta1$G_at_theta %*% checked_solve(theta1_scores$A_theta, t(theta1_scores$s_theta)))
    correction_p0 <- -t(theta0$G_at_theta %*% checked_solve(theta0_scores$A_theta, t(theta0_scores$s_theta)))
    if (!is.null(beta_scores)) {
      G_1t_beta <- crossprod((indicators - q1) * (-A / pi^2), dpi) / n
      G_0t_beta <- crossprod((indicators - q0) * ((1 - A) / (1 - pi)^2), dpi) / n
      correction_p1 <- correction_p1 - t(G_1t_beta %*% checked_solve(beta_scores$A_beta, t(beta_scores$s_beta)))
      correction_p0 <- correction_p0 - t(G_0t_beta %*% checked_solve(beta_scores$A_beta, t(beta_scores$s_beta)))
    }
    
    phi_corrected[, t] <- phi[, t] + as.vector(correction_p1 %*% (W %*% p0[t, ])) +
      as.vector(correction_p0 %*% (t(W) %*% p1[t, ]))
  }
  
  base$phi_corrected <- phi_corrected
  base$time <- Sys.time() - start
  base
}


#' Cross-fitted AIPW estimation of the DOOR trajectory
#'
#' Uses each participant's held-out nuisance predictions to construct augmented
#' category contributions. Pools these contributions across participants before
#' calculating the DOOR trajectory and its ordinary influence values.
#'
#' With nuis_CF$varc = TRUE, the proposition correction combines validation-fold
#' G matrices with the corresponding training scores and Jacobians.
#' The corollary correction estimates C_v on the validation fold. 
#' I_beta_v = -A_beta_v and score contributions are from the training fold, 
#' weighted by n_v / (n - n_v).
AIPW_CF <- function(obs_data, nuis_CF, K = dim(nuis_CF$q1)[3], n_timepoints = dim(nuis_CF$q1)[1],
                    outcome_columns = nuis_CF$outcome_columns, correction = c("proposition", "corollary")) {
  start <- Sys.time()
  correction <- match.arg(correction)
  if (is.null(outcome_columns)) outcome_columns <- paste0("PrimaryDOOROutcome_collapsed_", seq_len(n_timepoints) * 5)
  outcome_columns <- outcome_columns[seq_len(n_timepoints)]
  
  n <- nrow(obs_data); A <- obs_data$Treatment; pi <- nuis_CF$pi
  Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
  folds <- nuis_CF$folds
  fold_ids <- sort(unique(folds))
  W <- outer(seq_len(K), seq_len(K), function(k, l) (k > l) + 0.5 * (k == l))
  p1 <- p0 <- matrix(0, n_timepoints, K)
  psi <- numeric(n_timepoints)
  phi <- matrix(0, n, n_timepoints)
  
  # Each participant uses their held-out nuisance predictions
  for (t in seq_len(n_timepoints)) {
    q1 <- matrix(nuis_CF$q1[t, , ], n, K)
    q0 <- matrix(nuis_CF$q0[t, , ], n, K)
    indicators <- outer(Y_mat[, t], seq_len(K), "==")
    m1 <- q1 + (A / pi) * (indicators - q1)
    m0 <- q0 + ((1 - A) / (1 - pi)) * (indicators - q0)
    
    p1[t, ] <- colMeans(m1)
    p0[t, ] <- colMeans(m0)
    psi[t] <- as.numeric(p1[t, ] %*% W %*% p0[t, ])
    phi[, t] <- as.vector(sweep(m1, 2, p1[t, ], FUN = "-") %*% (W %*% p0[t, ])) +
      as.vector(sweep(m0, 2, p0[t, ], FUN = "-") %*% (t(W) %*% p1[t, ]))
  }
  
  base <- list(psi = psi, phi = phi, p1 = p1, p0 = p0, time = Sys.time() - start)
  if (!nuis_CF$varc) return(base)
  
  # Corollary: pair each training score with its own information estimate
  # Corollary: estimate C_v on validation rows; apply correction to training rows
  if (correction == "corollary") {
    phi_corrected <- phi
    
    for (v in fold_ids) {
      key <- as.character(v)
      in_v <- which(folds == v)
      out_v <- which(folds != v)
      n_v <- length(in_v); n_out <- length(out_v)
      training_weight_v <- n_v / n_out
      beta_scores_v <- nuis_CF$nuisance_components_folds[[key]]$beta
      if (is.null(beta_scores_v)) stop("The corollary requires fitted propensity components; use the proposition for a fixed working propensity.")
      
      # Validation scores evaluated at the propensity model trained outside fold v
      X_beta_val <- beta_scores_v$X_beta[in_v, , drop = FALSE]
      raw_pi_val <- as.numeric(predict(beta_scores_v$fit, newdata = obs_data[in_v, , drop = FALSE], type = "response"))
      s_beta_val <- X_beta_val * (A[in_v] - raw_pi_val)
      s_beta_val <- sweep(s_beta_val, 2, colMeans(s_beta_val), FUN = "-")
      
      # Both phi_val and s_beta_val use nuisance fit v
      phi_val <- phi[in_v, , drop = FALSE]
      phi_val <- sweep(phi_val, 2, colMeans(phi_val), FUN = "-")
      C_v <- crossprod(phi_val, s_beta_val) / n_v
      
      # Information and score contributions from the corresponding training sample
      I_beta_v <- -beta_scores_v$A_beta
      s_beta_train <- beta_scores_v$s_beta[out_v, , drop = FALSE]
      s_beta_train <- sweep(s_beta_train, 2, colMeans(s_beta_train), FUN = "-")
      
      # Accumulate corrections only for participants used to fit beta_v
      phi_corrected[out_v, ] <- phi_corrected[out_v, , drop = FALSE] -
        training_weight_v * t(C_v %*% checked_solve(I_beta_v, t(s_beta_train)))
    }
    
    base$phi_corrected <- phi_corrected
    base$time <- Sys.time() - start
    return(base)
  }
  
  if (!nuis_CF$outcome_model %in% c("polr", "multinom"))
    stop("The proposition correction requires a polr or multinom outcome model.")
  
  # Proposition: validation-fold G blocks with training scores and Jacobians
  phi_corrected <- phi
  
  for (v in fold_ids) {
    key <- as.character(v)
    in_v <- which(folds == v)
    n_v <- length(in_v)
    training_weight_v <- n_v / (n - n_v)
    val_data <- obs_data[in_v, , drop = FALSE]
    A_v <- A[in_v]; pi_v <- pi[in_v]
    scores_v <- nuis_CF$nuisance_components_folds[[key]]
    beta_scores_v <- scores_v$beta
    
    if (!is.null(beta_scores_v)) {
      X_beta_v <- beta_scores_v$X_beta[in_v, , drop = FALSE]
      # Smooth derivative convention; active clipping is qualified in the appendix.
      dpi_v <- X_beta_v * (pi_v * (1 - pi_v))
    }
    h1_v <- 1 - A_v / pi_v
    h0_v <- 1 - (1 - A_v) / (1 - pi_v)
    
    for (t in seq_len(n_timepoints)) {
      q1_v <- matrix(nuis_CF$q1[t, in_v, ], n_v, K)
      q0_v <- matrix(nuis_CF$q0[t, in_v, ], n_v, K)
      Y_v <- Y_mat[in_v, t]
      indicators_v <- outer(Y_v, seq_len(K), "==")
      
      fit1 <- nuis_CF$fit_trt_folds[[key]][[t]]
      fit0 <- nuis_CF$fit_control_folds[[key]][[t]]
      theta1_v <- aipw_val_components(fit1, val_data, q1_v, h1_v, outcome_model = nuis_CF$outcome_model)
      theta0_v <- aipw_val_components(fit0, val_data, q0_v, h0_v, outcome_model = nuis_CF$outcome_model)
      
      theta1_scores_v <- scores_v$theta1[[t]]
      theta0_scores_v <- scores_v$theta0[[t]]
      correction_p1 <- -t(theta1_v$G_at_theta %*% checked_solve(theta1_scores_v$A_theta, t(theta1_scores_v$s_theta)))
      correction_p0 <- -t(theta0_v$G_at_theta %*% checked_solve(theta0_scores_v$A_theta, t(theta0_scores_v$s_theta)))
      if (!is.null(beta_scores_v)) {
        G_1t_beta_v <- crossprod((indicators_v - q1_v) * (-A_v / pi_v^2), dpi_v) / n_v
        G_0t_beta_v <- crossprod((indicators_v - q0_v) * ((1 - A_v) / (1 - pi_v)^2), dpi_v) / n_v
        correction_p1 <- correction_p1 - t(G_1t_beta_v %*% checked_solve(beta_scores_v$A_beta, t(beta_scores_v$s_beta)))
        correction_p0 <- correction_p0 - t(G_0t_beta_v %*% checked_solve(beta_scores_v$A_beta, t(beta_scores_v$s_beta)))
      }
      
      phi_corrected[, t] <- phi_corrected[, t] + training_weight_v *
        (as.vector(correction_p1 %*% (W %*% p0[t, ])) + as.vector(correction_p0 %*% (t(W) %*% p1[t, ])))
    }
  }
  
  base$phi_corrected <- phi_corrected
  base$time <- Sys.time() - start
  base
}
