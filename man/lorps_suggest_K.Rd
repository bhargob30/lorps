\name{lorps_suggest_K}
\alias{lorps_suggest_K}
\title{Choose a Truncation Level from the LoRPS Prior Bound}
\description{
Calculates a sufficient finite truncation level for a specified upper bound
on the L1 discrepancy between the finite and infinite prior joint marginal
densities in Theorem 1 of the LoRPS manuscript.
}
\usage{lorps_suggest_K(N, Q, tolerance = 1e-3)}
\arguments{
  \item{N}{Positive integer number of rows.}
  \item{Q}{Positive integer number of columns.}
  \item{tolerance}{Target L1 bound, strictly between 0 and 2.}
}
\value{An integer \code{K}, at least 2.}
\details{
The theorem gives
\deqn{\|m^{(K)} - m^{(\infty)}\|_1 \le 2NQ\,2^{-(K-1)}.}
Consequently, the helper returns
\deqn{\max\{2,\ 1 + \lceil\log_2(2NQ/\epsilon)\rceil\},}
where \eqn{\epsilon} is \code{tolerance}. The final stick is forced to one.
The bound concerns prior joint marginal densities under the stated model.

The package default \code{K = 20} follows the manuscript's practical
recommendation. It is independent of this helper;
the helper does not silently alter the fit. In a 348 by 574 matrix,
\code{K = 30} makes the raw bound approximately 0.000744.
}
\references{
Kakoty, B., Chekouo, T., and Lock, E.
\emph{Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices}.
Manuscript.
}
\examples{
lorps_suggest_K(348, 574)  # 30
K <- lorps_suggest_K(100, 50, tolerance = 0.01)
2 * 100 * 50 * 2^(-(K - 1))
}
