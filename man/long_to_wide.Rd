\name{long_to_wide}
\alias{long_to_wide}
\title{Convert participant-visit data to one row per participant}
\description{Preserves outcome classes and checks that baseline variables are constant within participants.}
\usage{long_to_wide(data, id, time, outcome, baseline_columns,
  times = sort(unique(data[[time]])),
  outcome_columns = paste0(outcome, "_", times))}
\arguments{
\item{data}{Long-format data frame.}
\item{id, time, outcome}{Column names for participant, visit, and outcome.}
\item{baseline_columns}{Names of baseline columns to retain, including treatment.}
\item{times}{Visit values to select in the desired order.}
\item{outcome_columns}{Names of the resulting outcome columns.}
}
\value{A data frame with one row per participant. Missing visits become missing values; the fitting API requires complete analysis data.}
