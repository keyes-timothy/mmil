% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{mmil_adjust_predictions_case_control}
\alias{mmil_adjust_predictions_case_control}
\title{Update the predictions of an MMIL model using the case-control (i.e. sampling
bias) adjustment}
\usage{
mmil_adjust_predictions_case_control(
  base_predictions,
  z,
  rho,
  zeta,
  prediction_type = c("logit", "prob")
)
}
\arguments{
\item{base_predictions}{Unadjusted predictions (either logit or probability values)
from an MMIL model.}

\item{z}{A numeric vector of inherited patient labels for each cell. 1 if the
cell was sampled from a patient with the disease and 0 otherwise.}

\item{rho}{A scalar value between 0 and 1 indicating the probability that
that a cell is NOT associated with the outcome-of-interest, given that the
cell was sampled from a patient who DOES have the outcome of interest. In other
words, the expected number of non-disease-associated cells in each sick patient.}

\item{zeta}{A scalar value between 0 and 1 indicating the marginal probability that an
inherited patient label will be 1. In other words, the probability
that a random cell from your population-of-interest will be sampled from a patient with
the condition-of-interest.}

\item{prediction_type}{A string indicating how to interpret `base_predictions`.
Valid options are "logit" (the default), which interprets `base_predictions` as a
vector of logit values, and "prob", which interprets `base_predictions` as a
vector of probability values. Note that `prediction_type` also specifies the
output of this function (either adjusted logit or probability scores, respectively).}
}
\value{
If `prediction_type` == "logit", a vector of adjusted logit values. If
`prediction_type` == "prob", a vector of adjusted probability values.
}
\description{
Update the predictions of an MMIL model using the case-control (i.e. sampling
bias) adjustment
}
\examples{
NULL

}
