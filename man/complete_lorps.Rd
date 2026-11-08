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
# Generate a complete matrix with a low-rank underlying signal.
set.seed(10)

N <- 30
Q <- 20
signal_rank <- 2

U <- matrix(rnorm(N * signal_rank), nrow = N)
V <- matrix(rnorm(Q * signal_rank), nrow = Q)

signal <- tcrossprod(U, V)
X_true <- signal + matrix(rnorm(N * Q, sd = 0.3), nrow = N)

# Randomly mask 20 percent of the entries.
missing_index <- sort(
  sample.int(length(X_true), size = floor(0.20 * length(X_true)))
)

X <- X_true
X[missing_index] <- NA_real_

# Fit LoRPS with a fixed tensor rank.
# These are illustrative settings; convergence should be assessed.
fit <- lorps_mi(
  X,
  R = 2,
  K = 10,
  M = 2000,
  burn_in = 1000,
  thin = 5,
  seed = 12,
  verbose = FALSE
)

# Construct one completed matrix.
one <- complete_lorps(fit, action = 1)

# Construct a list containing the first and fifth completed matrices.
two <- complete_lorps(fit, action = c(1, 5))

# Replace missing entries with their means across imputations.
point_estimate <- complete_lorps(fit, action = "mean")

# Retrieve the original incomplete matrix.
original <- complete_lorps(fit, action = 0)

# Evaluate point accuracy using only the artificially masked entries.
true_missing <- X_true[fit$missing_index]
imputed_mean <- point_estimate[fit$missing_index]

nmse <- mean((imputed_mean - true_missing)^2) / var(true_missing)
print(c(NMSE = nmse))

# Check that the original matrix and observed values are preserved.
stopifnot(
  identical(original, fit$X),
  length(two) == 2L,
  identical(one[!is.na(X)], X[!is.na(X)]),
  identical(point_estimate[!is.na(X)], X[!is.na(X)])
)

\dontrun{
# Construct all completed matrices.
# For large datasets, processing one completion at a time uses less memory.
all_completed <- complete_lorps(fit, action = "all")
length(all_completed)
}
}
\seealso{\code{\link{lorps_mi}}}
