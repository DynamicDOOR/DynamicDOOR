# TMLE likelihood evaluation: "original" computes probabilities before logs;
# "log" uses normalized log probabilities. Both use the same targeting model.
#' tilt_q1: Applies an exponential tilt to treatment-arm probabilities, with
#' category K as the reference and eps scaled by 1 / pi_hat.
#' Returns updated probabilities for all rows, each summing to one.
tilt_q1 <- function(eps, pi_hat, q1){
  n <- dim(q1)[1]; K <- dim(q1)[2]
  log_weight <- log(q1)
  log_weight[,1:(K-1)] <- log_weight[,1:(K-1),drop=FALSE] + outer(1 / pi_hat, eps)
  weight <- exp(log_weight - apply(log_weight, 1, max))
  weight / rowSums(weight)
}

#' tilt_arm: Performs treatment-arm TMLE targeting using tilt_q1.
#' Estimates epsilon using the treated observations and BFGS optimization,
#' then updates treatment-arm category probabilities for all rows.
#' Returns epsilon, targeted probabilities q1_epsilon, the gradient of
#' the negative log-likelihood, and the optimizer convergence code.
tilt_arm <- function(q1, pi_hat, A, Y, tol = 1e-8, maxit = 100,
                     tmle_scale = "original"){
  tmle_scale <- match.arg(tmle_scale, choices = c("original", "log"))
  n <- dim(q1)[1]; K <- dim(q1)[2]
  idx <- cbind(seq_len(n), Y) 
  if (!setequal(Y[A == 1], seq_len(K)))
    warning("TMLE treatment targeting has an absent category; continuing optimization, but a finite unique unpenalized solution is not guaranteed",
            call. = FALSE)
  
  negll <- function(eps){
    if (tmle_scale == "original") {
      q1_epsilon <- tilt_q1(eps, pi_hat, q1)
      return(-sum(log(q1_epsilon[idx][A == 1])) / n)
    }
    # Normalize on the log scale to avoid probability underflow.
    log_weight <- log(q1)
    log_weight[,1:(K-1)] <- log_weight[,1:(K-1),drop=FALSE] + outer(1 / pi_hat, eps)
    shifted <- log_weight - apply(log_weight, 1, max)
    log_q1_epsilon <- shifted - log(rowSums(exp(shifted)))
    -sum(log_q1_epsilon[idx][A == 1]) / n
  }
  grad <- function(eps){
    q1_epsilon <- tilt_q1(eps, pi_hat, q1)
    g <- numeric(K-1)
    for(k in 1:(K-1)){
      g[k] <- -mean( A / pi_hat * ( (Y == k) - q1_epsilon[,k] ) )
    }
    return(g)
  }
  
  fit <- optim(par = rep(0, K-1), fn = negll, gr = grad, method = "BFGS", control = list(reltol = tol, maxit = maxit))
  q1_epsilon <- tilt_q1(fit$par, pi_hat, q1)
  score <- grad(fit$par) 
  
  list( epsilon = fit$par, q1_epsilon = q1_epsilon, score = score, convergence = fit$convergence )
}

#' tilt_q0: Applies an exponential tilt to control-arm probabilities, with
#' category K as the reference and eps scaled by 1 / (1 - pi_hat).
#' Returns updated probabilities for all rows, each summing to one.
tilt_q0 <- function(eps, pi_hat, q0){
  n <- dim(q0)[1]; K <- dim(q0)[2]
  log_weight <- log(q0)
  log_weight[,1:(K-1)] <- log_weight[,1:(K-1),drop=FALSE] + outer(1 / (1 - pi_hat), eps)
  row_max <- apply(log_weight, 1, max)
  weight <- exp(log_weight - row_max)
  weight / rowSums(weight)
}

#' tilt_control: Performs control-arm TMLE targeting using tilt_q0.
#' Estimates epsilon using the control observations and BFGS optimization,
#' then updates control-arm category probabilities for all rows.
#' Returns epsilon, targeted probabilities q0_epsilon, the gradient of
#' the negative log-likelihood, and the optimizer convergence code.
tilt_control <- function(q0, pi_hat, A, Y, tol = 1e-8, maxit = 100,
                         tmle_scale = "original"){
  tmle_scale <- match.arg(tmle_scale, choices = c("original", "log"))
  n <- dim(q0)[1]; K <- dim(q0)[2]
  idx <- cbind(seq_len(n), Y) 
  if (!setequal(Y[A == 0], seq_len(K)))
    warning("TMLE control targeting has an absent category; continuing optimization, but a finite unique unpenalized solution is not guaranteed",
            call. = FALSE)
  negll <- function(eps){
    if (tmle_scale == "original") {
      q0_epsilon <- tilt_q0(eps, pi_hat, q0)
      return(-sum(log(q0_epsilon[idx][A == 0])) / n)
    }
    # Normalize on the log scale to avoid probability underflow.
    log_weight <- log(q0)
    log_weight[,1:(K-1)] <- log_weight[,1:(K-1),drop=FALSE] + outer(1 / (1 - pi_hat), eps)
    shifted <- log_weight - apply(log_weight, 1, max)
    log_q0_epsilon <- shifted - log(rowSums(exp(shifted)))
    -sum(log_q0_epsilon[idx][A == 0]) / n
  }
  grad <- function(eps){
    q0_epsilon <- tilt_q0(eps, pi_hat, q0)
    g <- numeric(K-1)
    for(k in 1:(K-1)){
      g[k] <- -mean( (1-A) / (1-pi_hat) * ( (Y == k) - q0_epsilon[,k] ) )
    }
    return(g)
  }
  
  fit <- optim(par = rep(0, K-1), fn = negll, gr = grad, method = "BFGS", control = list(reltol = tol, maxit = maxit))
  q0_epsilon <- tilt_q0(fit$par, pi_hat, q0)
  
  score <- grad(fit$par) 
  list( epsilon = fit$par, q0_epsilon = q0_epsilon, score = score, convergence = fit$convergence )
}



#' TMLE derivative blocks for one arm, visit, and nuisance fit
#'
#' Evaluates the targeting scores and derivatives on the supplied evaluation rows.
#' With blocks_only = FALSE, also constructs the evaluation and training
#' contributions for the parametric proposition. Training scores retain their
#' full participant indexing; the caller supplies cross-fitting weights.
#' With blocks_only = TRUE, returns the blocks needed for pooled targeting.
#' A fixed propensity has no beta parameters: its derivative blocks have zero columns.
tmle_arm_score_components <- function(obs_data_val, Y_val, arm, pi_val, fit_q,
                                      q_at_initial, q_at_epsilon, epsilon, X_beta_val,
                                      outcome_model = c("polr", "multinom"),
                                      beta_scores, theta_scores, blocks_only = FALSE,
                                      p_at = colMeans(q_at_epsilon)) {
  outcome_model <- match.arg(outcome_model)
  n_val <- nrow(obs_data_val); K <- ncol(q_at_epsilon); r <- K - 1L
  arm_indicator_val <- as.numeric(obs_data_val$Treatment == arm)
  if (is.null(X_beta_val)) X_beta_val <- matrix(0, n_val, 0L)
  d_beta <- ncol(X_beta_val)
  pi_a_val <- if (arm == 1) pi_val else 1 - pi_val
  h_val <- 1 / pi_a_val

  # Smooth derivative evaluated at clipped pi, as documented in the appendix.
  # This is not the exact derivative of clipping on rows where clipping is active.
  dpi <- X_beta_val * (pi_val * (1 - pi_val))
  dh_val <- if (arm == 1) -dpi / pi_a_val^2 else dpi / pi_a_val^2
  dq_initial <- switch(outcome_model,
    polr = polr_probability_derivatives(fit_q, obs_data_val, q_at_initial),
    multinom = multinom_probability_derivatives(fit_q, obs_data_val, q_at_initial)
  )
  d_theta <- dim(dq_initial)[3]
  dq_theta <- array(0, c(n_val, K, d_theta))
  dq_beta <- array(0, c(n_val, K, d_beta))
  dq_eps <- array(0, c(n_val, K, r))

  for (i in seq_len(n_val)) {
    J_q <- diag(q_at_epsilon[i, ], nrow = K) - tcrossprod(q_at_epsilon[i, ])
    dlogq_theta <- matrix(dq_initial[i, , ], K, d_theta) / q_at_initial[i, ]
    dq_theta[i, , ] <- J_q %*% dlogq_theta
    dq_eps[i, , ] <- h_val[i] * J_q[, seq_len(r), drop = FALSE]
    dq_beta[i, , ] <- tcrossprod(as.vector(J_q[, seq_len(r), drop = FALSE] %*% epsilon), dh_val[i, ])
  }

  indicators <- outer(Y_val, seq_len(r), "==")
  resid <- indicators - q_at_epsilon[, seq_len(r), drop = FALSE]
  s_epsilon_at <- resid * (arm_indicator_val * h_val)
  A_at_epsilon <- matrix(0, r, r)
  A_at_epsilon_beta <- matrix(0, r, d_beta)
  A_at_epsilon_theta <- matrix(0, r, d_theta)

  for (i in seq_len(n_val)) {
    in_arm <- arm_indicator_val[i]
    A_at_epsilon <- A_at_epsilon - in_arm * h_val[i] * matrix(dq_eps[i, seq_len(r), ], r, r)
    A_at_epsilon_theta <- A_at_epsilon_theta - in_arm * h_val[i] * matrix(dq_theta[i, seq_len(r), ], r, d_theta)
    A_at_epsilon_beta <- A_at_epsilon_beta + in_arm * tcrossprod(resid[i, ], dh_val[i, ]) -
      in_arm * h_val[i] * matrix(dq_beta[i, seq_len(r), ], r, d_beta)
  }

  A_at_epsilon <- A_at_epsilon / n_val
  A_at_epsilon_beta <- A_at_epsilon_beta / n_val
  A_at_epsilon_theta <- A_at_epsilon_theta / n_val
  G_at_beta <- apply(dq_beta, c(2, 3), mean)
  G_at_theta <- apply(dq_theta, c(2, 3), mean)
  G_at_epsilon <- apply(dq_eps, c(2, 3), mean)

  blocks <- list(G_at_beta = G_at_beta, G_at_theta = G_at_theta, G_at_epsilon = G_at_epsilon,
                 A_at_epsilon = A_at_epsilon, A_at_epsilon_beta = A_at_epsilon_beta,
                 A_at_epsilon_theta = A_at_epsilon_theta, s_epsilon_at = s_epsilon_at)
  if (blocks_only) return(blocks)

  if (is.null(theta_scores) || (d_beta > 0L && is.null(beta_scores)))
    stop("The TMLE proposition correction requires training components for each estimated nuisance parameter.")

  B_at <- t(checked_solve(t(A_at_epsilon), t(G_at_epsilon)))
  phi_at_evaluation <- sweep(q_at_epsilon, 2, p_at, FUN = "-") - s_epsilon_at %*% t(B_at)
  L_at_beta <- G_at_beta - B_at %*% A_at_epsilon_beta
  L_at_theta <- G_at_theta - B_at %*% A_at_epsilon_theta
  phi_at_training <- -t(L_at_theta %*% checked_solve(theta_scores$A_theta, t(theta_scores$s_theta)))
  if (!is.null(beta_scores))
    phi_at_training <- phi_at_training - t(L_at_beta %*% checked_solve(beta_scores$A_beta, t(beta_scores$s_beta)))

  c(list(phi_at_evaluation = phi_at_evaluation, phi_at_training = phi_at_training,
         B_at = B_at, L_at_beta = L_at_beta, L_at_theta = L_at_theta), blocks)
}

#' Full-sample TMLE estimation of the DOOR trajectory
#'
#' Fits one fluctuation per arm and visit, averages targeted category
#' probabilities, and applies the DOOR matrix product. Ordinary influence
#' values use augmented contributions evaluated at the targeted predictions.
#' The proposition uses complete parametric score blocks; the corollary uses
#' only propensity scores and I_beta = -A_beta under its stated assumptions.
#' Corrections require nuis$varc = TRUE. compute_correction = FALSE permits
#' reuse of the targeting step by TMLE_CF_pooled().
#' Returns psi, phi, p1, p0, targeted predictions, fluctuations, targeting
#' diagnostics, elapsed time, and phi_corrected when correction is enabled.
TMLE <- function(obs_data, nuis, K = dim(nuis$q1)[3], n_timepoints = dim(nuis$q1)[1],
                 tol = 1e-8, maxit = 1000, outcome_columns = nuis$outcome_columns,
                 compute_correction = TRUE, correction = c("proposition", "corollary"), tmle_scale = "original") {
  start <- Sys.time()
  correction <- match.arg(correction)
  if (is.null(outcome_columns)) outcome_columns <- paste0("PrimaryDOOROutcome_collapsed_", seq_len(n_timepoints) * 5)
  outcome_columns <- outcome_columns[seq_len(n_timepoints)]
  n <- nrow(obs_data); A <- obs_data$Treatment; pi <- nuis$pi
  Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
  W <- outer(seq_len(K), seq_len(K), function(k, l) (k > l) + 0.5 * (k == l))
  q1_epsilon <- q0_epsilon <- array(0, c(n_timepoints, n, K))
  epsilon1 <- epsilon0 <- matrix(0, n_timepoints, K - 1L)
  target_score_max <- 0; target_convergence <- 0L
  p1 <- p0 <- matrix(0, n_timepoints, K)
  psi <- numeric(n_timepoints)
  phi <- matrix(0, n, n_timepoints)

  for (t in seq_len(n_timepoints)) {
    q1_initial <- matrix(nuis$q1[t, , ], n, K)
    q0_initial <- matrix(nuis$q0[t, , ], n, K)
    res1 <- tilt_arm(q1_initial, pi, A, Y_mat[, t], tol = tol, maxit = maxit, tmle_scale = tmle_scale)
    res0 <- tilt_control(q0_initial, pi, A, Y_mat[, t], tol = tol, maxit = maxit, tmle_scale = tmle_scale)
    q1 <- res1$q1_epsilon; q0 <- res0$q0_epsilon
    q1_epsilon[t, , ] <- q1
    q0_epsilon[t, , ] <- q0
    epsilon1[t, ] <- res1$epsilon
    epsilon0[t, ] <- res0$epsilon
    target_score_max <- max(target_score_max, abs(res1$score), abs(res0$score))
    target_convergence <- max(target_convergence, res1$convergence, res0$convergence)

    p1[t, ] <- colMeans(q1)
    p0[t, ] <- colMeans(q0)
    psi[t] <- as.numeric(p1[t, ] %*% W %*% p0[t, ])
    indicators <- outer(Y_mat[, t], seq_len(K), "==")
    m1 <- q1 + (A / pi) * (indicators - q1)
    m0 <- q0 + ((1 - A) / (1 - pi)) * (indicators - q0)
    phi[, t] <- as.vector(sweep(m1, 2, p1[t, ], FUN = "-") %*% (W %*% p0[t, ])) +
      as.vector(sweep(m0, 2, p0[t, ], FUN = "-") %*% (t(W) %*% p1[t, ]))
  }

  base <- list(psi = psi, phi = phi, p1 = p1, p0 = p0, q1_epsilon = q1_epsilon, q0_epsilon = q0_epsilon,
               epsilon1 = epsilon1, epsilon0 = epsilon0, target_score_max = target_score_max,
               target_convergence = target_convergence, time = Sys.time() - start)
  if (!compute_correction || !nuis$varc) return(base)
  beta_scores <- nuis$nuisance_components$beta

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

  phi_corrected <- matrix(0, n, n_timepoints)
  for (t in seq_len(n_timepoints)) {
    arm1 <- tmle_arm_score_components(
      obs_data, Y_mat[, t], 1, pi, nuis$fit_trt[[t]],
      matrix(nuis$q1[t, , ], n, K), matrix(q1_epsilon[t, , ], n, K), epsilon1[t, ],
      X_beta_val = beta_scores$X_beta, outcome_model = nuis$outcome_model,
      beta_scores = beta_scores, theta_scores = nuis$nuisance_components$theta1[[t]], p_at = p1[t, ])
    arm0 <- tmle_arm_score_components(
      obs_data, Y_mat[, t], 0, pi, nuis$fit_control[[t]],
      matrix(nuis$q0[t, , ], n, K), matrix(q0_epsilon[t, , ], n, K), epsilon0[t, ],
      X_beta_val = beta_scores$X_beta, outcome_model = nuis$outcome_model,
      beta_scores = beta_scores, theta_scores = nuis$nuisance_components$theta0[[t]], p_at = p0[t, ])

    phi_corrected[, t] <- as.vector((arm1$phi_at_evaluation + arm1$phi_at_training) %*% (W %*% p0[t, ])) +
      as.vector((arm0$phi_at_evaluation + arm0$phi_at_training) %*% (t(W) %*% p1[t, ]))
  }

  base$phi_corrected <- phi_corrected
  base$time <- Sys.time() - start
  base
}

#' Cross-fitted TMLE with separate targeting in each validation fold
#'
#' Pools targeted category probabilities before calculating the DOOR trajectory.
#' psi_fold is retained as a diagnostic. With nuis_CF$varc = TRUE, the proposition
#' combines validation targeting blocks with corresponding training scores.
#' The corollary estimates C_v on validation rows at nuisance fit v, then applies
#' the corresponding training scores and I_beta_v = -A_beta_v with weight
#' n_v / (n - n_v). Its assumptions are the same as for the AIPW corollary.
TMLE_CF <- function(obs_data, nuis_CF, K = dim(nuis_CF$q1)[3], n_timepoints = dim(nuis_CF$q1)[1],
                    tol = 1e-8, maxit = 1000, outcome_columns = nuis_CF$outcome_columns,
                    correction = c("proposition", "corollary"), tmle_scale = "original") {
  start <- Sys.time()
  correction <- match.arg(correction)
  if (is.null(outcome_columns)) outcome_columns <- paste0("PrimaryDOOROutcome_collapsed_", seq_len(n_timepoints) * 5)
  outcome_columns <- outcome_columns[seq_len(n_timepoints)]
  n <- nrow(obs_data); A <- obs_data$Treatment; pi <- nuis_CF$pi
  Y_mat <- vapply(obs_data[outcome_columns], function(x) as.integer(as.character(x)), integer(n))
  folds <- nuis_CF$folds; fold_ids <- sort(unique(folds)); V <- length(fold_ids)
  W <- outer(seq_len(K), seq_len(K), function(k, l) (k > l) + 0.5 * (k == l))
  q1_epsilon <- q0_epsilon <- array(0, c(n_timepoints, n, K))
  epsilon1 <- epsilon0 <- array(0, c(V, n_timepoints, K - 1L))
  psi_fold <- matrix(0, V, n_timepoints)
  target_score_max <- 0; target_convergence <- 0L

  # Separate targeting in each validation fold
  for (v_index in seq_along(fold_ids)) {
    in_v <- which(folds == fold_ids[v_index])
    n_v <- length(in_v)
    for (t in seq_len(n_timepoints)) {
      q1_initial <- matrix(nuis_CF$q1[t, in_v, ], n_v, K)
      q0_initial <- matrix(nuis_CF$q0[t, in_v, ], n_v, K)
      res1 <- tilt_arm(q1_initial, pi[in_v], A[in_v], Y_mat[in_v, t], tol = tol, maxit = maxit, tmle_scale = tmle_scale)
      res0 <- tilt_control(q0_initial, pi[in_v], A[in_v], Y_mat[in_v, t], tol = tol, maxit = maxit, tmle_scale = tmle_scale)
      q1_epsilon[t, in_v, ] <- res1$q1_epsilon
      q0_epsilon[t, in_v, ] <- res0$q0_epsilon
      epsilon1[v_index, t, ] <- res1$epsilon
      epsilon0[v_index, t, ] <- res0$epsilon
      target_score_max <- max(target_score_max, abs(res1$score), abs(res0$score))
      target_convergence <- max(target_convergence, res1$convergence, res0$convergence)
      psi_fold[v_index, t] <- as.numeric(colMeans(res1$q1_epsilon) %*% W %*% colMeans(res0$q0_epsilon))
    }
  }

  # Pool targeted margins; calculate ordinary influence values explicitly
  p1 <- p0 <- matrix(0, n_timepoints, K)
  psi <- numeric(n_timepoints)
  phi <- matrix(0, n, n_timepoints)
  for (t in seq_len(n_timepoints)) {
    q1 <- matrix(q1_epsilon[t, , ], n, K)
    q0 <- matrix(q0_epsilon[t, , ], n, K)
    p1[t, ] <- colMeans(q1)
    p0[t, ] <- colMeans(q0)
    psi[t] <- as.numeric(p1[t, ] %*% W %*% p0[t, ])
    indicators <- outer(Y_mat[, t], seq_len(K), "==")
    m1 <- q1 + (A / pi) * (indicators - q1)
    m0 <- q0 + ((1 - A) / (1 - pi)) * (indicators - q0)
    phi[, t] <- as.vector(sweep(m1, 2, p1[t, ], FUN = "-") %*% (W %*% p0[t, ])) +
      as.vector(sweep(m0, 2, p0[t, ], FUN = "-") %*% (t(W) %*% p1[t, ]))
  }

  base <- list(psi = psi, phi = phi, p1 = p1, p0 = p0, psi_fold = psi_fold,
               q1_epsilon = q1_epsilon, q0_epsilon = q0_epsilon, epsilon1 = epsilon1, epsilon0 = epsilon0,
               fold_ids = fold_ids, target_score_max = target_score_max,
               target_convergence = target_convergence, time = Sys.time() - start)
  if (!nuis_CF$varc) return(base)

  # Corollary: validation covariance, paired training information and scores
  if (correction == "corollary") {
    phi_corrected <- phi
    for (v in fold_ids) {
      key <- as.character(v)
      in_v <- which(folds == v); out_v <- which(folds != v)
      n_v <- length(in_v); n_out <- length(out_v)
      beta_scores_v <- nuis_CF$nuisance_components_folds[[key]]$beta
      if (is.null(beta_scores_v)) stop("The corollary requires fitted propensity components; use the proposition for a fixed working propensity.")

      X_beta_val <- beta_scores_v$X_beta[in_v, , drop = FALSE]
      raw_pi_val <- as.numeric(predict(beta_scores_v$fit, newdata = obs_data[in_v, , drop = FALSE], type = "response"))
      s_beta_val <- X_beta_val * (A[in_v] - raw_pi_val)
      s_beta_val <- sweep(s_beta_val, 2, colMeans(s_beta_val), FUN = "-")
      phi_val <- phi[in_v, , drop = FALSE]
      phi_val <- sweep(phi_val, 2, colMeans(phi_val), FUN = "-")
      C_v <- crossprod(phi_val, s_beta_val) / n_v

      I_beta_v <- -beta_scores_v$A_beta
      s_beta_train <- beta_scores_v$s_beta[out_v, , drop = FALSE]
      s_beta_train <- sweep(s_beta_train, 2, colMeans(s_beta_train), FUN = "-")
      phi_corrected[out_v, ] <- phi_corrected[out_v, , drop = FALSE] -
        (n_v / n_out) * t(C_v %*% checked_solve(I_beta_v, t(s_beta_train)))
    }
    base$phi_corrected <- phi_corrected
    base$time <- Sys.time() - start
    return(base)
  }

  if (!nuis_CF$outcome_model %in% c("polr", "multinom"))
    stop("The proposition correction requires a polr or multinom outcome model.")

  # Proposition: separate B_at and L_at blocks for each targeting fit
  phi_corrected <- matrix(0, n, n_timepoints)
  for (v_index in seq_along(fold_ids)) {
    v <- fold_ids[v_index]; key <- as.character(v)
    in_v <- which(folds == v)
    n_v <- length(in_v); training_weight_v <- n_v / (n - n_v)
    scores_v <- nuis_CF$nuisance_components_folds[[key]]
    beta_scores_v <- scores_v$beta

    for (t in seq_len(n_timepoints)) {
      arm1 <- tmle_arm_score_components(
        obs_data[in_v, , drop = FALSE], Y_mat[in_v, t], 1, pi[in_v], nuis_CF$fit_trt_folds[[key]][[t]],
        matrix(nuis_CF$q1[t, in_v, ], n_v, K), matrix(q1_epsilon[t, in_v, ], n_v, K), epsilon1[v_index, t, ],
        X_beta_val = beta_scores_v$X_beta[in_v, , drop = FALSE], outcome_model = nuis_CF$outcome_model,
        beta_scores = beta_scores_v, theta_scores = scores_v$theta1[[t]], p_at = p1[t, ])
      arm0 <- tmle_arm_score_components(
        obs_data[in_v, , drop = FALSE], Y_mat[in_v, t], 0, pi[in_v], nuis_CF$fit_control_folds[[key]][[t]],
        matrix(nuis_CF$q0[t, in_v, ], n_v, K), matrix(q0_epsilon[t, in_v, ], n_v, K), epsilon0[v_index, t, ],
        X_beta_val = beta_scores_v$X_beta[in_v, , drop = FALSE], outcome_model = nuis_CF$outcome_model,
        beta_scores = beta_scores_v, theta_scores = scores_v$theta0[[t]], p_at = p0[t, ])

      phi_corrected[in_v, t] <- phi_corrected[in_v, t] +
        as.vector(arm1$phi_at_evaluation %*% (W %*% p0[t, ])) +
        as.vector(arm0$phi_at_evaluation %*% (t(W) %*% p1[t, ]))
      phi_corrected[, t] <- phi_corrected[, t] + training_weight_v *
        (as.vector(arm1$phi_at_training %*% (W %*% p0[t, ])) +
           as.vector(arm0$phi_at_training %*% (t(W) %*% p1[t, ])))
    }
  }

  base$phi_corrected <- phi_corrected
  base$time <- Sys.time() - start
  base
}


#' Cross-fitted TMLE-P with pooled targeting
#'
#' Fits one common fluctuation per arm and visit over all held-out predictions.
#' The proposition pools G_at_epsilon and A_at_epsilon before forming their
#' matrix ratio, while retaining fold-specific nuisance training corrections.
#' The corollary uses validation covariances and the corresponding training
#' propensity scores and information, as in TMLE_CF().
TMLE_CF_pooled <- function(obs_data, nuis_CF, K = dim(nuis_CF$q1)[3], n_timepoints = dim(nuis_CF$q1)[1],
                           tol = 1e-8, maxit = 1000, outcome_columns = nuis_CF$outcome_columns,
                           correction = c("proposition", "corollary"), tmle_scale = "original") {
  start <- Sys.time()
  correction <- match.arg(correction)
  if (is.null(outcome_columns)) outcome_columns <- paste0("PrimaryDOOROutcome_collapsed_", seq_len(n_timepoints) * 5)
  outcome_columns <- outcome_columns[seq_len(n_timepoints)]

  # Common targeting over all held-out predictions
  fit <- TMLE(obs_data, nuis_CF, K, n_timepoints, tol, maxit, outcome_columns,
              compute_correction = FALSE, tmle_scale = tmle_scale)
  fit$targeting <- "pooled_out_of_fold"
  if (!nuis_CF$varc) {
    fit$time <- Sys.time() - start
    return(fit)
  }

  n <- nrow(obs_data); A <- obs_data$Treatment
  folds <- nuis_CF$folds; fold_ids <- sort(unique(folds))
  fold_proportions <- vapply(fold_ids, function(v) mean(folds == v), numeric(1))

  # Corollary: validation covariance, paired training information and scores
  if (correction == "corollary") {
    phi_corrected <- fit$phi
    for (v in fold_ids) {
      key <- as.character(v)
      in_v <- which(folds == v); out_v <- which(folds != v)
      n_v <- length(in_v); n_out <- length(out_v)
      beta_scores_v <- nuis_CF$nuisance_components_folds[[key]]$beta
      if (is.null(beta_scores_v)) stop("The corollary requires fitted propensity components; use the proposition for a fixed working propensity.")

      X_beta_val <- beta_scores_v$X_beta[in_v, , drop = FALSE]
      raw_pi_val <- as.numeric(predict(beta_scores_v$fit, newdata = obs_data[in_v, , drop = FALSE], type = "response"))
      s_beta_val <- X_beta_val * (A[in_v] - raw_pi_val)
      s_beta_val <- sweep(s_beta_val, 2, colMeans(s_beta_val), FUN = "-")
      phi_val <- fit$phi[in_v, , drop = FALSE]
      phi_val <- sweep(phi_val, 2, colMeans(phi_val), FUN = "-")
      C_v <- crossprod(phi_val, s_beta_val) / n_v

      I_beta_v <- -beta_scores_v$A_beta
      s_beta_train <- beta_scores_v$s_beta[out_v, , drop = FALSE]
      s_beta_train <- sweep(s_beta_train, 2, colMeans(s_beta_train), FUN = "-")
      phi_corrected[out_v, ] <- phi_corrected[out_v, , drop = FALSE] -
        (n_v / n_out) * t(C_v %*% checked_solve(I_beta_v, t(s_beta_train)))
    }
    fit$phi_corrected <- phi_corrected
    fit$time <- Sys.time() - start
    return(fit)
  }

  if (!nuis_CF$outcome_model %in% c("polr", "multinom"))
    stop("The proposition correction requires a polr or multinom outcome model.")

  W <- outer(seq_len(K), seq_len(K), function(k, l) (k > l) + 0.5 * (k == l))
  fit$phi_corrected <- matrix(0, n, n_timepoints)

  for (t in seq_len(n_timepoints)) {
    Y <- as.integer(as.character(obs_data[[outcome_columns[t]]]))
    phi_p_corrected <- vector("list", 2L)

    for (a in 0:1) {
      q_at_epsilon <- matrix(if (a == 1) fit$q1_epsilon[t, , ] else fit$q0_epsilon[t, , ], n, K)
      q_at_initial <- matrix(if (a == 1) nuis_CF$q1[t, , ] else nuis_CF$q0[t, , ], n, K)
      epsilon_at <- if (a == 1) fit$epsilon1[t, ] else fit$epsilon0[t, ]
      p_at <- if (a == 1) fit$p1[t, ] else fit$p0[t, ]

      # Fold-specific derivatives at the shared fitted fluctuation
      blocks_at_v <- lapply(fold_ids, function(v) {
        key <- as.character(v)
        in_v <- which(folds == v)
        scores_v <- nuis_CF$nuisance_components_folds[[key]]
        fit_at_v <- if (a == 1) nuis_CF$fit_trt_folds[[key]][[t]] else nuis_CF$fit_control_folds[[key]][[t]]
        theta_scores_at_v <- if (a == 1) scores_v$theta1[[t]] else scores_v$theta0[[t]]
        tmle_arm_score_components(
          obs_data[in_v, , drop = FALSE], Y[in_v], a, nuis_CF$pi[in_v], fit_at_v,
          q_at_initial[in_v, , drop = FALSE], q_at_epsilon[in_v, , drop = FALSE], epsilon_at,
          X_beta_val = scores_v$beta$X_beta[in_v, , drop = FALSE], outcome_model = nuis_CF$outcome_model,
          beta_scores = scores_v$beta, theta_scores = theta_scores_at_v, blocks_only = TRUE)
      })

      # Pool the two matrices separately before solving
      G_at_epsilon_pooled <- Reduce(`+`, Map(function(blocks, weight) weight * blocks$G_at_epsilon,
                                           blocks_at_v, fold_proportions))
      A_at_epsilon_pooled <- Reduce(`+`, Map(function(blocks, weight) weight * blocks$A_at_epsilon,
                                           blocks_at_v, fold_proportions))
      B_at_TMLE_P <- t(checked_solve(t(A_at_epsilon_pooled), t(G_at_epsilon_pooled)))

      phi_at_corrected <- sweep(q_at_epsilon, 2, p_at, FUN = "-")
      for (v_index in seq_along(fold_ids)) {
        v <- fold_ids[v_index]; key <- as.character(v)
        in_v <- which(folds == v)
        n_v <- length(in_v); training_weight_v <- n_v / (n - n_v)
        scores_v <- nuis_CF$nuisance_components_folds[[key]]
        beta_scores_v <- scores_v$beta
        theta_scores_at_v <- if (a == 1) scores_v$theta1[[t]] else scores_v$theta0[[t]]
        blocks_at <- blocks_at_v[[v_index]]

        phi_at_corrected[in_v, ] <- phi_at_corrected[in_v, , drop = FALSE] -
          blocks_at$s_epsilon_at %*% t(B_at_TMLE_P)
        L_at_beta_v <- blocks_at$G_at_beta - B_at_TMLE_P %*% blocks_at$A_at_epsilon_beta
        L_at_theta_v <- blocks_at$G_at_theta - B_at_TMLE_P %*% blocks_at$A_at_epsilon_theta
        phi_at_corrected <- phi_at_corrected - training_weight_v *
          t(L_at_theta_v %*% checked_solve(theta_scores_at_v$A_theta, t(theta_scores_at_v$s_theta)))
        if (!is.null(beta_scores_v))
          phi_at_corrected <- phi_at_corrected - training_weight_v *
            t(L_at_beta_v %*% checked_solve(beta_scores_v$A_beta, t(beta_scores_v$s_beta)))
      }
      phi_p_corrected[[a + 1L]] <- phi_at_corrected
    }

    fit$phi_corrected[, t] <- as.vector(phi_p_corrected[[2]] %*% (W %*% fit$p0[t, ])) +
      as.vector(phi_p_corrected[[1]] %*% (t(W) %*% fit$p1[t, ]))
  }

  fit$time <- Sys.time() - start
  fit
}
