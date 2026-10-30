\name{complete_lorps}
\alias{complete_lorps}
\title{Construct Completed Matrices from a LoRPS Fit}
\description{
Reconstructs completed datasets on demand, preserving observed values and
the joint alignment of the imputed entries in each retained draw.
}
\usage{complete_lorps(object, action = 1L)}
\arguments{
  \item{object}{An object returned by \code{\link{lorps_mi}}.}
  \item{action}{One or more retained-draw numbers from 1 to
  \code{object$n_imputations}; 0 for the original incomplete matrix;
  \code{"all"} for all completed matrices; or \code{"mean"} for a single
  matrix with missing entries replaced by their posterior means.}
}
\value{
A numeric matrix if one draw, 0, or \code{"mean"} is requested.
A named list of numeric matrices for multiple draw numbers or \code{"all"}.
}
\details{
Each draw number identifies one entire column of \code{object$imputations}.
Do not independently permute draw columns across missing entries: that
would destroy their joint alignment. Requesting \code{"all"} materializes
all complete matrices and may need much more memory than the compact fit.
For memory-efficient downstream analysis, analyze one completed draw at a
time and retain only its estimates and covariance matrix.

The \code{"mean"} completion is a point summary and does not propagate
imputation uncertainty. For nonlinear quantities, calculate the quantity
within each completed draw and then summarize across draws.

For compatible downstream models, an ordinary list of completed data frames
can be passed to \code{mitools::imputationList}; complete-data estimates and
variances can then be combined with \code{mitools::MIcombine}.
This function does not fit downstream models or pool estimates itself.
}
\examples{
set.seed(21)
X <- matrix(rnorm(48), 8, 6)
X[c(1, 9, 23)] <- NA
fit <- lorps_mi(X, R = 1, K = 2, M = 30, burn_in = 10,
                seed = 22, verbose = FALSE)
one <- complete_lorps(fit, 1)
two <- complete_lorps(fit, c(1, 5))
point_estimate <- complete_lorps(fit, "mean")
original <- complete_lorps(fit, 0)
stopifnot(identical(original, fit$X), length(two) == 2L)
}
\seealso{\code{\link{lorps_mi}}}
