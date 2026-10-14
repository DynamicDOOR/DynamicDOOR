# Simple longitudinal toy data generation.
cum_probs <- function(alpha, eta) {
  sapply(alpha, function(a) pnorm(a + eta))
}

to_cat_probs <- function(Fcum, K) {
  Fcum <- as.matrix(Fcum)
  cbind(
    Fcum[, 1, drop = FALSE],
    if (K > 2) {
      sapply(2:(K - 1), function(k) Fcum[, k] - Fcum[, k - 1])
    } else {
      NULL
    },
    1 - Fcum[, K - 1]
  )
}


gen_toy_data <- function(n, K = 4, eif = FALSE, seed = 0, true_values = NULL, Sigma_time=NULL, propensity_scale = 1) {
  stopifnot(length(propensity_scale) == 1L, is.finite(propensity_scale), propensity_scale > 0, K >= 2L)
  set.seed(seed)
  
  # Generate the covariates: two binary, two continuous and one nominal
  X1 <- rbinom(n, size = 1, prob = 0.20)
  X2 <- rbinom(n, size = 1, prob = 0.60)
  X3 <- factor( sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.40, 0.35, 0.25)), levels = c("A", "B", "C") )
  X3B <- as.numeric(X3 == "B")
  X3C <- as.numeric(X3 == "C")
  X4 <- rnorm(n, mean = -1, sd = 1)
  X5 <- rnorm(n, mean = 1, sd = 1)
  
  times <- c(5, 10, 15, 20, 25, 30)
  n_timepoints <- length(times)
  alpha0 <- seq(-1, 1, length.out = K - 1)
  alpha1 <- seq(-1, 1, length.out = K - 1)
  long_dat <- data.frame(
    id = rep(seq_len(n), times = n_timepoints),
    X1 = rep(X1, times = n_timepoints),
    X2 = rep(X2, times = n_timepoints),
    X3 = rep(X3, times = n_timepoints),
    X4 = rep(X4, times = n_timepoints),
    X5 = rep(X5, times = n_timepoints),
    time = factor(rep(times, each = n), levels = times)
  )
  
  # Day 5 is the reference. These are coefficients for the factor levels.
  time_coef0 <- c( `5` = 0.00, `10` = 0.12, `15` = 0.24, `20` = 0.36, `25` = 0.48, `30` = 0.60 )
  time_coef1 <- c( `5` = 0.00, `10` = 0.08, `15` = 0.16, `20` = 0.24, `25` = 0.32, `30` = 0.40 )
  
  # control
  eta0_base <- 0.20 + 0.20 * X1 + 0.10 * X2 + 0.30 * X3B - 0.30 * X3C + 1.20 * X4 + 0.70 * X5
  eta0_long <- rep(eta0_base, times = n_timepoints) + unname(time_coef0[as.character(long_dat$time)])
  
  # treatment
  eta1_base <- 0.50 + 0.50 * X1 + 0.30 * X2 + 0.25 * X3B - 0.25 * X3C + 0.90 * X4 - 0.50 * X5
  eta1_long <- rep(eta1_base, times = n_timepoints) + unname(time_coef1[as.character(long_dat$time)])
  
  pY0 <- array(0, dim = c(n_timepoints, n, K))
  pY1 <- array(0, dim = c(n_timepoints, n, K))
  Y0 <- Y1 <- matrix(0, nrow = n, ncol = n_timepoints)
  if (is.null(Sigma_time)){
    Sigma_time <- outer( seq_len(n_timepoints),  seq_len(n_timepoints), function(j, k) 0.60^abs(j - k) )
  }
  
  set.seed(seed + 1000L)
  error0 <- MASS::mvrnorm( n = n, mu = rep(0, n_timepoints), Sigma = Sigma_time )
  set.seed(seed + 2000L)
  error1 <- MASS::mvrnorm( n = n, mu = rep(0, n_timepoints), Sigma = Sigma_time )
  
  for (t in seq_len(n_timepoints)) {
    rows_t <- long_dat$time == times[t]
    eta0 <- eta0_long[rows_t]; eta1 <- eta1_long[rows_t]
    Fcum0 <- cum_probs(alpha0, eta0)
    Fcum1 <- cum_probs(alpha1, eta1)
    pY0[t, , ] <- to_cat_probs(Fcum0, K)
    pY1[t, , ] <- to_cat_probs(Fcum1, K)
    
    pY0[t, , ] <- pmax(pY0[t, , ], 0)
    pY1[t, , ] <- pmax(pY1[t, , ], 0)
    pY0[t, , ] <- pY0[t, , ] / rowSums(pY0[t, , ])
    pY1[t, , ] <- pY1[t, , ] / rowSums(pY1[t, , ])
    
    Y0[, t] <- cut( -eta0 + error0[, t], breaks = c(-Inf, alpha0, Inf), labels = FALSE, right = TRUE )
    Y1[, t] <- cut( -eta1 + error1[, t], breaks = c(-Inf, alpha1, Inf), labels = FALSE, right = TRUE )
  }
  
  pA <- plogis(propensity_scale * ( 0.80 + 0.80 * X1 - 0.50 * X2 +  0.40 * X3B - 0.40 * X3C + 1.00 * X4 + 0.30 * X5))
  set.seed(seed + 101L)
  A <- rbinom(n, 1, pA)
  
  # Potential outcome framework
  observed_dat <- data.frame(X1, X2, X3, X4, X5, Treatment = A)
  for (t in seq_len(n_timepoints)) {
    Y <- ifelse(A == 1, Y1[, t], Y0[, t])
    observed_dat[[paste0("PrimaryDOOROutcome_collapsed_", times[t])]] <-
      factor(Y, levels = seq_len(K))
  }
  long_dat$Treatment <- rep(A, times = n_timepoints)
  long_dat$Y0 <- factor(as.vector(Y0), levels = seq_len(K))
  long_dat$Y1 <- factor(as.vector(Y1), levels = seq_len(K))
  long_dat$Y <- ifelse( long_dat$Treatment == 1, as.integer(long_dat$Y1), as.integer(long_dat$Y0) )
  long_dat$Y <- factor(long_dat$Y, levels = seq_len(K))
  
  result <- list( observed_data = observed_dat, long_data = long_dat, propensity = pA, 
                  qY1 = pY1, qY0 = pY0,  K = K, alpha0 = alpha0,  alpha1 = alpha1 )
  
  if (eif) {
    if (
      is.null(true_values) || !all(c("true_psi", "true_p0", "true_p1") %in% names(true_values))
    ) {
      stop("true_values must contain true_psi, true_p0, and true_p1.")
    }
    true_psi <- true_values$true_psi
    true_p0 <- true_values$true_p0
    true_p1 <- true_values$true_p1
    phi <- matrix(0, nrow = n, ncol = n_timepoints)
    for (t in seq_len(n_timepoints)) {
      p0_t <- true_p0[t, ]
      p1_t <- true_p1[t, ]
      F0 <- cumsum(p0_t)
      F1 <- cumsum(p1_t)
      M0 <- M1 <- numeric(K)
      
      for (k in seq_len(K)) {
        M0[k] <- if (k == 1) { 0.5 * p0_t[k] } else {
          F0[k - 1] + 0.5 * p0_t[k]
        }
        M1[k] <- (1 - F1[k]) + 0.5 * p1_t[k]
      }
      
      eta <- as.vector(pY1[t, , ] %*% M0)
      H0 <- as.vector(pY0[t, , ] %*% M1)
      Y <- as.integer(as.character(
        observed_dat[[paste0("PrimaryDOOROutcome_collapsed_", times[t])]]
      ))
      phi[, t] <- A / pA * (M0[Y] - eta) + (1 - A) / (1 - pA) * (M1[Y] - H0) + eta + H0 - 2 * true_psi[t]
    }
    
    result$phi_mat <- phi
    result$sd_eif <- apply(phi, 2, sd)
    result$sd_eif_avg <- sd(rowMeans(phi))
  }
  
  result
}

true_psi_p_gen <- function(M, K = 4, seed = 0, propensity_scale = 1) {
  dat <- gen_toy_data(M, K = K, seed = seed, propensity_scale = propensity_scale)
  n_timepoints <- dim(dat$qY0)[1]
  true_p0 <- true_p1 <- matrix(0, nrow = n_timepoints, ncol = K)
  true_psi <- numeric(n_timepoints)
  
  for (t in seq_len(n_timepoints)) {
    true_p0[t, ] <- colMeans(dat$qY0[t, , ])
    true_p1[t, ] <- colMeans(dat$qY1[t, , ])
    mat <- outer(true_p1[t, ], true_p0[t, ])
    true_psi[t] <- sum(mat[lower.tri(mat)]) + 0.5 * sum(diag(mat))
  }
  
  list( true_psi = true_psi, true_psi_avg = mean(true_psi), true_p0 = true_p0, true_p1 = true_p1, p0 = true_p0, p1 = true_p1 )
}
