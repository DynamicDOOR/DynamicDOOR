\name{gen_toy_data}
\alias{gen_toy_data}
\alias{true_psi_p_gen}
\title{Generate synthetic repeated ordinal outcomes and reference probabilities}
\description{Generates six ordinal outcomes at days 5 through 30 with baseline confounding. Larger categories are preferable. Reference probabilities integrate the known outcome probabilities over simulated covariates.}
\usage{
gen_toy_data(n, K = 4, eif = FALSE, seed = 0, true_values = NULL,
  Sigma_time = NULL, propensity_scale = 1)
true_psi_p_gen(M, K = 4, seed = 0, propensity_scale = 1)
}
\arguments{
\item{n}{Number of participants.}
\item{M}{Monte Carlo sample size for reference values.}
\item{K}{Number of categories, at least two. The fitting API requires at least three.}
\item{eif}{Whether to calculate influence contributions at the supplied truth.}
\item{seed}{Random seed.}
\item{true_values}{Reference values from \code{true_psi_p_gen}, required when \code{eif = TRUE}.}
\item{Sigma_time}{Optional six-by-six latent covariance matrix.}
\item{propensity_scale}{Positive multiplier of the treatment linear predictor.}
}
\value{\code{gen_toy_data} returns observed wide and long data, generating probabilities, and optionally influence contributions. \code{true_psi_p_gen} returns marginal category probabilities and visit-specific and average win probabilities.}
