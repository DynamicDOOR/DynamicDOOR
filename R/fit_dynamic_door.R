#' Fit a dynamic DOOR trajectory and its weighted temporal average
#'
#' Larger categories are more desirable. The target at each visit is
#' P(Y(1) > Y(0)) + 0.5 P(Y(1) = Y(0)), using independent marginal draws.
#' This function fits the nuisance models, estimates the trajectory, computes
#' ordinary and corrected inference, and optionally displays the trajectory.
#'
#' @param obs_data One row per independent participant, with complete analysis columns.
#' @param treatment Name of the treatment column, coded 0/1.
#' @param outcome_columns Outcome columns in visit order.
#' @param trt_predictors,outcome_predictors Character vectors of baseline covariate names.
#' @param K,n_timepoints Category and visit counts; inferred when omitted. K must be >= 3.
#' @param method "AIPW", "TMLE", or "TMLE-P". With cf = FALSE, AIPW and
#'   TMLE use full-sample fitting; cf = TRUE gives AIPW-CF and TMLE-CF.
#'   TMLE-P requires cf = TRUE and pools targeting across evaluation folds.
#' @param outcome Required outcome learner: "polr", "multinom", or "tree".
#' @param varc Required TRUE/FALSE: whether to compute corrected inference.
#' @param correction Required: NULL with varc = FALSE; otherwise "proposition",
#'   "corollary", or "bootstrap". "proposition" uses the parametric correction
#'   for polr/multinom; "corollary" uses the projection under a correctly specified
#'   fitted propensity; "bootstrap" uses bootstrap covariance. Tree outcomes
#'   cannot use "proposition". The correction is never selected automatically.
#' @param cf,V,folds Use cross-fitting, with V folds and optional supplied assignments.
#' @param seed Seed for fold assignment and bootstrap resampling.
#' @param times Numeric visit times; defaults to 1,...,n_timepoints.
#' @param outcome_levels Category labels from least to most desirable. Defaults to
#'   1:K, ordered-factor levels, or sorted observed numeric values.
#' @param weights Nonnegative temporal weights; normalized to sum to one.
#' @param ci_trajectory "pointwise", "simul", or "both". "simul" returns the
#'   chi-square joint Wald region and its coordinate projections for plotting.
#'   Joint membership uses the full inverse covariance, not the projected bounds.
#'   A singular covariance makes the joint region unavailable; no fallback is used.
#'   The weighted average always has a pointwise normal Wald interval.
#' @param conf_level Confidence level (default 0.95).
#' @param null_trajectory Scalar or vector for the joint Wald test; default 0.5 at every visit.
#' @param plot Display and retain a trajectory plot. Corrected intervals are used when varc = TRUE.
#' @param B Bootstrap replicates. Bootstrap estimates supply covariance only;
#'   intervals remain centered at the original estimate, with normal/chi-square critical values.
#' @param fixed_propensity Optional fixed working P(A = 1); omits propensity fitting
#'   and all beta-parameter correction terms. Incompatible with "corollary".
#' @param scale Standardize numeric predictors before fitting (default TRUE).
#' @param tree_control Named rpart.control overrides.
#' @param tol,maxit,tmle_scale TMLE controls; tmle_scale is "original" or "log".
#' @return A dynamic_door object: results, covariance, joint regions, joint_tests,
#'   fit, nuisance, settings, diagnostics, bootstrap, and plot. Results contain
#'   ordinary and corrected SEs, ci_lower/upper for pointwise inference, and
#'   simul_lower/upper for ellipsoid projections. Unrequested bounds are NA.
fit_dynamic_door <- function(obs_data, treatment, outcome_columns,
                             trt_predictors, outcome_predictors,
                             K = NULL, n_timepoints = NULL,
                             method = "AIPW", outcome, varc,
                             cf = TRUE, seed = 0, times = NULL,
                             outcome_levels = NULL, V = 5,
                             B = 100, conf_level = 0.95,
                             weights = NULL, tol = 1e-8, maxit = 1000,
                             scale = TRUE, tree_control = list(),
                             tmle_scale = "original", ci_trajectory = "pointwise",
                             plot = TRUE, correction = NULL,
                             fixed_propensity = NULL, folds = NULL, null_trajectory = 0.5) {
  call <- match.call()

  # Check the specification.
  if (missing(outcome) || missing(varc) || missing(correction))
    stop("Specify outcome, varc, and correction explicitly")
  if (!is.character(outcome) || length(outcome) != 1L || is.na(outcome))
    stop("outcome must be a single character string")
  if (identical(ci_trajectory, "simultaneous")) ci_trajectory <- "simul"
  ci_trajectory <- match.arg(ci_trajectory, c("pointwise", "simul", "both"))
  tmle_scale <- match.arg(tmle_scale, c("original", "log"))
  method <- match.arg(method, c("AIPW", "TMLE", "TMLE-P"))
  outcome <- match.arg(outcome, c("polr", "multinom", "tree"))
  if (!is.null(correction)) {
    if (!is.character(correction) || length(correction) != 1L || is.na(correction))
      stop("correction must be NULL or a single character string")
    correction <- match.arg(correction, c("proposition", "corollary", "bootstrap"))
  }
  for (flag in list(varc = varc, cf = cf, scale = scale, plot = plot))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) stop("varc, cf, scale, and plot must be TRUE/FALSE")
  if (method == "TMLE-P" && !cf) stop("TMLE-P requires cf = TRUE; use TMLE for full-sample targeting")
  if (identical(correction, "proposition") && outcome == "tree") stop("The proposition requires polr or multinom; it cannot be used with tree")
  if (varc == is.null(correction))
    stop("Use correction = NULL with varc = FALSE; otherwise set varc = TRUE and choose a correction")
  if (!is.null(fixed_propensity)) {
    if (!is.numeric(fixed_propensity) || length(fixed_propensity) != 1L || !is.finite(fixed_propensity) ||
        fixed_propensity <= 0 || fixed_propensity >= 1) stop("fixed_propensity must be between zero and one")
    if (identical(correction, "corollary")) stop("The corollary requires an estimated, correctly specified propensity")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be between zero and one")
  if (!is.character(treatment) || length(treatment) != 1L || is.na(treatment) ||
      !is.character(trt_predictors) || !is.character(outcome_predictors) || anyNA(c(trt_predictors, outcome_predictors)))
    stop("Treatment and predictor arguments must contain column names")
  if (!is.character(outcome_columns) || !length(outcome_columns) || anyNA(outcome_columns) || anyDuplicated(outcome_columns))
    stop("Supply distinct outcome_columns in visit order")
  if (is.null(n_timepoints)) n_timepoints <- length(outcome_columns)
  if (!is.numeric(n_timepoints) || length(n_timepoints) != 1L || is.na(n_timepoints) || n_timepoints != length(outcome_columns))
    stop("n_timepoints must equal length(outcome_columns)")
  if (is.null(times)) times <- seq_len(n_timepoints)
  if (!is.numeric(times) || length(times) != n_timepoints || any(!is.finite(times)) || any(diff(times) <= 0))
    stop("times must be finite, increasing, and match outcome_columns")
  predictors <- union(trt_predictors, outcome_predictors)
  required <- unique(c(treatment, outcome_columns, predictors))
  if (is.null(weights)) weights <- rep(1 / n_timepoints, n_timepoints)
  if (!is.numeric(weights) || length(weights) != n_timepoints || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) stop("weights must be nonnegative with positive sum, one per visit")
  weights <- weights / sum(weights)
  if (cf && (length(V) != 1L || !is.finite(V) || V < 2 || V != floor(V))) stop("V must be an integer of at least two")
  if (identical(correction, "bootstrap") && (length(B) != 1L || !is.finite(B) || B < 2 || B != floor(B)))
    stop("B must be an integer of at least two")
  if (length(null_trajectory) == 1L) null_trajectory <- rep(null_trajectory, n_timepoints)
  if (length(null_trajectory) != n_timepoints || any(!is.finite(null_trajectory)))
    stop("null_trajectory must be scalar or one value per visit")

  # Check the data used for fitting.
  obs_data <- as.data.frame(obs_data)
  if (!is.data.frame(obs_data) || nrow(obs_data) < 2L || anyDuplicated(names(obs_data)))
    stop("obs_data must be a data frame with at least two subjects and unique column names")
  if (!all(required %in% names(obs_data)))
    stop("Missing analysis columns: ", paste(setdiff(required, names(obs_data)), collapse = ", "))
  if (treatment %in% outcome_columns || any(predictors %in% c(treatment, outcome_columns)))
    stop("Treatment, outcomes, and baseline predictors must use distinct columns")
  if ("Treatment" %in% c(outcome_columns, predictors))
    stop("Please avoid the name `Treatment` in outcome and predictors")
  if (any(outcome_predictors %in% c("Y", "time", "id")))
    stop("Please avoid the name `Y`, `time`, or `id` in the data")
  for (column in required) {
    x <- obs_data[[column]]
    invalid <- if (is.numeric(x)) any(!is.finite(x)) else anyNA(x)
    if (invalid) stop("Missing or nonfinite values in column: ", column)
  }
  A <- as.character(obs_data[[treatment]])
  if (!all(A %in% c("0", "1"))) stop("Treatment must be coded 0/1")
  if (is.null(outcome_levels)) {
    if (!is.null(K)) {
      if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 3 || K != floor(K))
        stop("K must be an integer of at least three")
      outcome_levels <- as.character(seq_len(K))
    } else if (is.ordered(obs_data[[outcome_columns[1L]]])) {
      outcome_levels <- levels(obs_data[[outcome_columns[1L]]])
    } else {
      observed <- unique(unlist(lapply(obs_data[outcome_columns], as.character)))
      numeric_levels <- suppressWarnings(as.numeric(observed))
      if (anyNA(numeric_levels)) stop("Supply outcome_levels from least to most desirable")
      outcome_levels <- observed[order(numeric_levels)]
    }
  }
  outcome_levels <- as.character(outcome_levels)
  if (anyNA(outcome_levels) || anyDuplicated(outcome_levels) || length(outcome_levels) < 3L)
    stop("outcome_levels must contain at least three distinct ordered labels")
  if (is.null(K)) K <- length(outcome_levels)
  if (!is.numeric(K) || length(K) != 1L || is.na(K) || K != length(outcome_levels))
    stop("K must be numeric and equal to length(outcome_levels)")
  
  # Prepare the analysis data.
  if (!is.null(seed)) set.seed(seed)
  
  dat <- obs_data
  dat$Treatment <- as.integer(A)
  for (column in outcome_columns) {
    y <- match(as.character(obs_data[[column]]), outcome_levels)
    if (anyNA(y)) stop("Outcome values outside outcome_levels in: ", column)
    dat[[column]] <- factor(y, levels = seq_len(K))
  }
  for (column in predictors) {
    if (is.character(dat[[column]])) dat[[column]] <- factor(dat[[column]])
    if (is.factor(dat[[column]])) dat[[column]] <- droplevels(dat[[column]])
    if (scale && is.numeric(dat[[column]])) {
      if (stats::sd(dat[[column]]) == 0) stop("Constant predictor: ", column)
      dat[[column]] <- as.numeric(base::scale(dat[[column]]))
    }
  }
  n <- nrow(dat)
  
  if (cf) {
    if (is.null(folds)) folds <- kfold_indices_cs(n, V, seed)
    if (length(folds) != n || anyNA(folds) || length(unique(folds)) != V)
      stop("folds must assign each participant to one of V nonempty folds")
  } else folds <- NULL
  support_issues <- nuisance_support_complete(dat, if (is.null(fixed_propensity)) trt_predictors else character(), outcome_predictors, K, V,
                                              seed, folds, cf = cf, outcome_columns = outcome_columns,
                                              details = TRUE, outcome_levels = outcome_levels)
  if (length(support_issues))
    stop(paste(c("Insufficient data support before nuisance fitting:",
                 paste0("- ", support_issues)), collapse = "\n"), call. = FALSE)

  # Fit the nuisance models.
  analytic <- varc && correction %in% c("proposition", "corollary")
  if (cf) {
    nuisance <- fit_nuisance_ind_CF(dat, V, trt_predictors, outcome_predictors, K, n_timepoints,
                                    seed = seed, outcome_model = outcome, folds = folds, varc = analytic,
                                    outcome_columns = outcome_columns, tree_control = tree_control,
                                    fixed_propensity = fixed_propensity)
  } else {
    nuisance <- fit_nuisance_ind(dat, trt_predictors, outcome_predictors, K, n_timepoints,
                                 outcome_model = outcome, varc = analytic, outcome_columns = outcome_columns,
                                 tree_control = tree_control, fixed_propensity = fixed_propensity)
  }
  # Use the original direct switch, with pooled targeting added for TMLE-P.
  # Current estimators read varc from nuisance; it is not a positional argument.
  fit_correction <- if (analytic) correction else "proposition" # placeholder when varc=FALSE
  fit <- switch(if (method == "TMLE-P") "TMLE-P" else paste0(method, if (cf) "_CF" else ""),
    AIPW = AIPW(dat, nuisance, K, n_timepoints, outcome_columns, correction = fit_correction),
    AIPW_CF = AIPW_CF(dat, nuisance, K, n_timepoints, outcome_columns, correction = fit_correction),
    TMLE = TMLE(dat, nuisance, K, n_timepoints, tol, maxit, outcome_columns,
                correction = fit_correction, tmle_scale = tmle_scale),
    TMLE_CF = TMLE_CF(dat, nuisance, K, n_timepoints, tol, maxit, outcome_columns,
                      correction = fit_correction, tmle_scale = tmle_scale),
    "TMLE-P" = TMLE_CF_pooled(dat, nuisance, K, n_timepoints, tol, maxit, outcome_columns,
                              correction = fit_correction, tmle_scale = tmle_scale))
  estimator <- if (method == "TMLE-P") "tmle_cf_pooled" else paste0(tolower(method), if (cf) "_cf" else "")

  # Keep the original estimate, influence, standard-error, and results calculations.
  if (!is.matrix(fit$phi) || !identical(dim(fit$phi), c(as.integer(n), as.integer(n_timepoints))) || any(!is.finite(fit$phi)))
    stop("Expected finite participant influence values with n rows and one column per visit")
  estimate <- c(fit$psi, sum(weights * fit$psi))
  influence <- cbind(fit$phi, as.vector(fit$phi %*% weights))
  se <- apply(influence, 2, stats::sd) / sqrt(n)
  z <- stats::qnorm((1 + conf_level) / 2)
  results <- data.frame(target = c(rep("time", n_timepoints), "average"),
                        time = c(times, NA_real_), estimate = estimate,
                        se = se, ci_lower = NA_real_, ci_upper = NA_real_,
                        simul_lower = NA_real_, simul_upper = NA_real_,
                        se_corrected = NA_real_, ci_lower_corrected = NA_real_, ci_upper_corrected = NA_real_,
                        simul_lower_corrected = NA_real_, simul_upper_corrected = NA_real_)
  covariance <- list(original = stats::cov(fit$phi) / n)
  bootstrap <- NULL
  if (varc && analytic) {
    if (!is.matrix(fit$phi_corrected) || !identical(dim(fit$phi_corrected), dim(fit$phi)) || any(!is.finite(fit$phi_corrected)))
      stop("Corrected influence values are unavailable; no substitution of ordinary inference is made")
    corrected <- cbind(fit$phi_corrected, as.vector(fit$phi_corrected %*% weights))
    results$se_corrected <- apply(corrected, 2, stats::sd) / sqrt(n)
    covariance$corrected <- stats::cov(fit$phi_corrected) / n
  } else if (varc) {
    bootstrap <- bootstrap_tree(dat, trt_predictors, outcome_predictors, B = B, V = V, K = K,
                                n_timepoints = n_timepoints, seed = seed, estimators = estimator,
                                tol = tol, maxit = maxit, folds = folds, outcome_columns = outcome_columns,
                                outcome_model = outcome, tree_control = tree_control, tmle_scale = tmle_scale,
                                fixed_propensity = fixed_propensity)
    boot <- bootstrap$estimates[[estimator]]
    if (!is.matrix(boot) || ncol(boot) != n_timepoints) stop("Bootstrap estimates must have one column per visit")
    valid <- apply(boot, 1, function(x) all(is.finite(x)))
    if (sum(valid) < 2L) stop("Fewer than two valid bootstrap replicates; corrected inference unavailable")
    if (!all(valid)) warning(sprintf("Only %d of %d bootstrap replicates were valid", sum(valid), length(valid)), call. = FALSE)
    boot <- cbind(boot[valid, , drop = FALSE], as.vector(boot[valid, , drop = FALSE] %*% weights))
    results$se_corrected <- apply(boot, 2, stats::sd)
    covariance$corrected <- stats::cov(boot[, seq_len(n_timepoints), drop = FALSE])
  }

  # The temporal summary always retains its pointwise normal Wald interval.
  point_rows <- if (ci_trajectory %in% c("pointwise", "both")) seq_len(n_timepoints + 1L) else n_timepoints + 1L
  results$ci_lower[point_rows] <- (estimate - z * se)[point_rows]
  results$ci_upper[point_rows] <- (estimate + z * se)[point_rows]
  if (varc) {
    results$ci_lower_corrected[point_rows] <- (estimate - z * results$se_corrected)[point_rows]
    results$ci_upper_corrected[point_rows] <- (estimate + z * results$se_corrected)[point_rows]
  }

  # Joint Wald region and its coordinate projections, using a chi-square critical value.
  joint <- list()
  joint_tests <- data.frame(variance = character(), statistic = numeric(), df = integer(),
                            critical = numeric(), p_value = numeric(), null_in_region = logical(), status = character())
  if (ci_trajectory %in% c("simul", "both")) {
    critical <- stats::qchisq(conf_level, df = n_timepoints)
    for (variance in names(covariance)) {
      cov_est <- covariance[[variance]]
      suffix <- if (variance == "corrected") "_corrected" else ""
      visit_se <- results[[paste0("se", suffix)]][seq_len(n_timepoints)]
      root <- if (all(is.finite(cov_est)) && rcond(cov_est) >= 1e-12)
        tryCatch(chol((cov_est + t(cov_est)) / 2), error = function(e) NULL) else NULL
      valid <- !is.null(root)
      precision <- if (valid) chol2inv(root) else NULL
      statistic <- if (valid) sum(forwardsolve(t(root), fit$psi - null_trajectory)^2) else NA_real_
      status <- if (valid) "ok" else "singular covariance; joint inference unavailable"
      joint[[variance]] <- list(center = fit$psi, covariance = cov_est, precision = precision, df = n_timepoints,
                                 critical = critical, conf_level = conf_level, null = null_trajectory, status = status)
      joint_tests <- rbind(joint_tests, data.frame(variance = variance, statistic = statistic, df = n_timepoints, critical = critical,
                          p_value = stats::pchisq(statistic, df = n_timepoints, lower.tail = FALSE),
                          null_in_region = statistic <= critical, status = status))
      if (valid) {
        results[[paste0("simul_lower", suffix)]][seq_len(n_timepoints)] <- fit$psi - sqrt(critical) * visit_se
        results[[paste0("simul_upper", suffix)]][seq_len(n_timepoints)] <- fit$psi + sqrt(critical) * visit_se
      } else warning(paste(variance, status), call. = FALSE)
    }
  }

  # Return estimation, inference, and diagnostics together.
  if (!is.null(fit$target_convergence) && fit$target_convergence != 0L)
    warning("Some TMLE targeting optimizations did not converge; inspect diagnostics", call. = FALSE)
  object <- structure(list(call = call, results = results, covariance = covariance, joint = joint,
    joint_tests = joint_tests, fit = fit, nuisance = nuisance,
    settings = list(method = method, estimator = estimator, outcome = outcome, outcome_model = outcome,
                    varc = varc, cf = cf, scale = scale, seed = seed, treatment = treatment,
                    tmle_scale = tmle_scale, tree_control = tree_control, outcome_columns = outcome_columns,
                    outcome_levels = outcome_levels, trt_predictors = trt_predictors, outcome_predictors = outcome_predictors,
                    K = K, n_timepoints = n_timepoints, n = n, times = times, weights = weights,
                    V = if (cf) V else NULL, folds = folds, B = if (!is.null(bootstrap)) B else NULL,
                    ci_trajectory = ci_trajectory, conf_level = conf_level, correction = correction,
                    fixed_propensity = fixed_propensity, probability_floor = probability_clip(), tol = tol, maxit = maxit),
    diagnostics = list(target_convergence = fit$target_convergence, target_score_max = fit$target_score_max,
                       clipping = nuisance$clipping, valid_B = bootstrap$valid_B, progress_log = bootstrap$progress_log),
    bootstrap = bootstrap, plot = NULL), class = "dynamic_door")
  if (plot) {
    object$plot <- plot.dynamic_door(object)
    print(object$plot)
  }
  object
}

#' Plot one or more dynamic_door objects, optionally with both types of intervals.
#' Corrected intervals are used by default when requested during fitting.
#' Ellipsoid projections are simultaneous bounds for a single trajectory;
#' joint tests use the full covariance in x$joint_tests.
plot.dynamic_door <- function(x, ..., corrected = x$settings$varc,
                              ci_trajectory = x$settings$ci_trajectory,
                              xlab = "Time", ylab = "Win probability") {
  ci_trajectory <- match.arg(ci_trajectory, c("pointwise", "simul", "both"))
  objects <- c(list(x), list(...))
  if (!all(vapply(objects, inherits, logical(1), "dynamic_door"))) stop("Every object must be a dynamic_door fit")
  labels <- names(objects)
  if (is.null(labels)) labels <- rep("", length(objects))
  for (i in seq_along(objects)) {
    object <- objects[[i]]
    if (is.na(labels[i]) || labels[i] == "")
      labels[i] <- paste0(object$settings$method, if (object$settings$cf && object$settings$method != "TMLE-P") "-CF" else "",
                          " / ", object$settings$outcome)
  }
  labels <- make.unique(labels)
  data <- vector("list", length(objects))
  for (i in seq_along(objects)) {
    object <- objects[[i]]
    if (corrected && !object$settings$varc) stop("Corrected intervals were not computed for every inference method")
    if (!identical(object$settings$outcome_levels, x$settings$outcome_levels) || object$settings$conf_level != x$settings$conf_level)
      stop("Overlaid fits must have the same outcome levels and confidence level")
    if (object$settings$ci_trajectory != "both" && ci_trajectory != object$settings$ci_trajectory)
      stop("Requested intervals were not computed; fit with ci_trajectory = 'both' to plot either type")
    d <- object$results[object$results$target == "time", , drop = FALSE]
    suffix <- if (corrected) "_corrected" else ""
    d$point_lower <- d[[paste0("ci_lower", suffix)]]; d$point_upper <- d[[paste0("ci_upper", suffix)]]
    d$simul_lower_plot <- d[[paste0("simul_lower", suffix)]]; d$simul_upper_plot <- d[[paste0("simul_upper", suffix)]]
    d$model <- labels[i]
    data[[i]] <- d
  }
  data <- do.call(rbind, data)
  spacing <- if (length(unique(data$time)) > 1L) min(diff(sort(unique(data$time)))) else 1
  dodge <- ggplot2::position_dodge(width = spacing * 0.2)
  p <- ggplot2::ggplot(data, ggplot2::aes(x = time, y = estimate, color = model, group = model)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60")
  if (ci_trajectory %in% c("simul", "both"))
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = simul_lower_plot, ymax = simul_upper_plot),
                                   position = dodge, width = spacing * 0.15, linewidth = 0.5, alpha = 0.4, na.rm = TRUE)
  if (ci_trajectory %in% c("pointwise", "both"))
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = point_lower, ymax = point_upper),
                                   position = dodge, width = spacing * 0.07, linewidth = 0.7, na.rm = TRUE)
  interval_label <- switch(ci_trajectory, pointwise = "pointwise Wald intervals",
                            simul = "joint Wald ellipsoid projections",
                            both = "dark: pointwise Wald intervals; light: joint Wald ellipsoid projections")
  p + ggplot2::geom_line(position = dodge, linewidth = 0.6) + ggplot2::geom_point(position = dodge, size = 2) +
    ggplot2::scale_x_continuous(breaks = sort(unique(data$time))) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = ceiling(length(objects) / 3))) +
    ggplot2::labs(x = xlab, y = ylab, color = NULL,
                  caption = paste0(100 * x$settings$conf_level, "% ", interval_label,
                    " (", if (corrected) "corrected" else "uncorrected", ").")) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = "top", panel.grid = ggplot2::element_blank())
}
