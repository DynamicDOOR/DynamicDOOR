\name{plot.dynamic_door}
\alias{plot.dynamic_door}
\title{Plot a fitted win-probability trajectory}
\description{Shows pointwise intervals, projected joint-region bounds, or both.}
\usage{\method{plot}{dynamic_door}(x, ..., corrected = x$settings$varc,
  ci_trajectory = x$settings$ci_trajectory,
  xlab = "Time", ylab = "Win probability")}
\arguments{
\item{x}{A fitted object from \code{fit_dynamic_door}.}
\item{...}{Additional fitted objects to overlay.}
\item{corrected}{Use corrected rather than ordinary covariance estimates.}
\item{ci_trajectory}{One of \code{"pointwise"}, \code{"simul"}, or \code{"both"}. Requested intervals must have been calculated during fitting.}
\item{xlab, ylab}{Axis labels.}
}
\value{A \code{ggplot} object.}
