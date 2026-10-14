\name{fit_dynamic_door}
\alias{fit_dynamic_door}
\title{Fit longitudinal ordinal win probabilities}
\description{ Larger categories are more desirable. The target at each visit is P(Y(1) > Y(0)) + 0.5 P(Y(1) = Y(0)), using independent marginal draws.  This function fits the nuisance models, estimates the trajectory, computes ordinary and corrected inference, and optionally displays the trajectory. }
\usage{fit_dynamic_door(obs_data, treatment, outcome_columns,
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
                             fixed_propensity = NULL, folds = NULL, null_trajectory = 0.5)}
\arguments{
\item{obs_data}{One row per independent participant, with complete analysis columns.}
\item{treatment}{Name of the treatment column, coded 0/1.}
\item{outcome_columns}{Outcome columns in visit order.}
\item{trt_predictors,outcome_predictors}{Character vectors of baseline covariate names.}
\item{K,n_timepoints}{Category and visit counts; inferred when omitted. K must be >= 3.}
\item{method}{"AIPW", "TMLE", or "TMLE-P". With cf = FALSE, AIPW and TMLE use full-sample fitting; cf = TRUE gives AIPW-CF and TMLE-CF. TMLE-P requires cf = TRUE and pools targeting across evaluation folds.}
\item{outcome}{Required outcome learner: "polr", "multinom", or "tree".}
\item{varc}{Required TRUE/FALSE: whether to compute corrected inference.}
\item{correction}{Required: NULL with varc = FALSE; otherwise "proposition", "corollary", or "bootstrap". "proposition" uses the parametric correction for polr/multinom; "corollary" uses the projection under a correctly specified fitted propensity; "bootstrap" uses bootstrap covariance. Tree outcomes cannot use "proposition". The correction is never selected automatically.}
\item{cf,V,folds}{Use cross-fitting, with V folds and optional supplied assignments.}
\item{seed}{Seed for fold assignment and bootstrap resampling.}
\item{times}{Numeric visit times; defaults to 1,...,n_timepoints.}
\item{outcome_levels}{Category labels from least to most desirable. Defaults to 1:K, ordered-factor levels, or sorted observed numeric values.}
\item{weights}{Nonnegative temporal weights; normalized to sum to one.}
\item{ci_trajectory}{"pointwise", "simul", or "both". "simul" returns the chi-square joint Wald region and its coordinate projections for plotting. Joint membership uses the full inverse covariance, not the projected bounds. A singular covariance makes the joint region unavailable; no fallback is used. The weighted average always has a pointwise normal Wald interval.}
\item{conf_level}{Confidence level (default 0.95).}
\item{null_trajectory}{Scalar or vector for the joint Wald test; default 0.5 at every visit.}
\item{plot}{Display and retain a trajectory plot. Corrected intervals are used when varc = TRUE.}
\item{B}{Bootstrap replicates. Bootstrap estimates supply covariance only; intervals remain centered at the original estimate, with normal/chi-square critical values.}
\item{fixed_propensity}{Optional fixed working P(A = 1); omits propensity fitting and all beta-parameter correction terms. Incompatible with "corollary".}
\item{scale}{Standardize numeric predictors before fitting (default TRUE).}
\item{tree_control}{Named rpart.control overrides.}
\item{tol,maxit,tmle_scale}{TMLE controls; tmle_scale is "original" or "log".}
}
\value{A dynamic_door object: results, covariance, joint regions, joint_tests, fit, nuisance, settings, diagnostics, bootstrap, and plot. Results contain ordinary and corrected SEs, ci_lower/upper for pointwise inference, and simul_lower/upper for ellipsoid projections. Unrequested bounds are NA.}
\examples{
data <- gen_toy_data(n = 300, K = 3, seed = 52)$observed_data
fit <- fit_dynamic_door(
  data, treatment = "Treatment",
  outcome_columns = "PrimaryDOOROutcome_collapsed_5",
  trt_predictors = paste0("X", 1:5), outcome_predictors = paste0("X", 1:5),
  outcome_levels = as.character(1:3), method = "AIPW", cf = FALSE,
  outcome = "polr", varc = FALSE, correction = NULL, plot = FALSE
)
fit$results
}
