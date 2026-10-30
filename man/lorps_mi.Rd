\name{lorps_mi}
\alias{lorps_mi}
\title{Multiple Imputation with a Low-Rank Probit Stick-Breaking Model}
\description{
Fits a normal mixture model with entry-specific probit stick-breaking weights
whose latent variables have a low-rank tensor representation. Missing entries
are stored as joint retained sampler draws. The tensor rank can be specified
or selected by held-out posterior predictive likelihood.
}
\usage{
lorps_mi(X, R = NULL, K = 20, CV = is.null(R),
         cand_rank_list = 1:15,
         M = 2000, burn_in = 1000, thin = 1,
         mus_init = NULL, sigma2s_init = NULL, C_ij_init = NULL,
         tau2 = 10000, cv_control = list(),
         seed = NULL, verbose = TRUE)
}
\arguments{
  \item{X}{A nonempty real numeric matrix, with individuals in rows and
  variables in columns. Missing values are \code{NA} or \code{NaN}.
  Infinite values and entirely missing matrices are not supported.
  A complete matrix is also accepted. No automatic scaling, transformation,
  or restriction of imputations to the observed range is applied.}
  \item{R}{A positive integer tensor rank, or \code{NULL}. With the defaults,
  \code{NULL} requests CV and an explicit value requests a fixed-rank fit.
  If \code{CV = FALSE} and \code{R = NULL}, rank 5 is used.
  Explicit \code{CV = TRUE} overrides a supplied \code{R}, with a warning.}
  \item{K}{Number of mixture components in the finite stick-breaking
  representation, at least 2. The default 20 is a computational starting
  choice. It is not a universal bound on imputation error. See
  \code{\link{lorps_suggest_K}} for a prior truncation bound.}
  \item{CV}{Logical. Select the rank by LPPL CV? The default is \code{is.null(R)}.}
  \item{cand_rank_list}{Positive integer candidate ranks. Used only for CV;
  duplicate values are removed and the values are sorted. Larger lists
  require more sampler fits. Exact ties select the smaller rank.}
  \item{M, burn_in, thin}{Total sampler iterations, burn-in setting, and
  thinning interval. \code{M} and \code{thin} must be positive integers;
  \code{burn_in} must be a nonnegative integer. Both \code{M} and
  \code{burn_in} must be multiples of \code{thin}. The sampler retains
  \code{(M - burn_in) / thin - 1} draws, which must be at least 2.
  Defaults give 999 draws. These settings do not guarantee convergence.}
  \item{mus_init, sigma2s_init, C_ij_init}{Optional matched starting values.
  Supply all three together or leave all three \code{NULL}. The first two
  are vectors of length \code{K}, containing finite means and positive
  finite variances. \code{C_ij_init} is a matrix with \code{dim(X)} and
  zero-based integer labels \code{0, ..., K - 1}, including positions missing
  in \code{X}. Custom starts apply to the final fit only; CV always initializes
  from its own training data.}
  \item{tau2}{Positive finite prior variance for component means. The sampler
  uses independent priors \eqn{\mu_k \sim N(0,\tau^2)} and
  \eqn{\sigma_k^2 \sim IG(2,1)}; the inverse-gamma parameters are fixed.}
  \item{cv_control}{A named list overriding any CV defaults:
  \describe{
    \item{\code{nfold = 5}}{Number of disjoint validation subsets, at least 2.}
    \item{\code{holdout_fraction = 0.05}}{Fraction of originally observed
    entries in each validation subset, rounded down to an integer count.}
    \item{\code{M = 1000, burn_in = 500, thin = 1}}{CV sampler settings,
    giving 499 retained draws by default. The same divisibility and
    retained-draw checks apply as for the final chain.}
    \item{\code{cores = NULL}}{Automatic use of up to five workers,
    capped by the number of subsets and \code{parallelly::availableCores()}.
    An integer supplies a different upper limit; 1 forces serial CV.
    A PSOCK cluster is used across platforms. If cluster startup fails,
    CV falls back to serial execution with a warning.}
  }}
  \item{seed}{Optional nonnegative integer. An explicit seed is local to this
  call and the previous R random-number state is restored on exit.
  With \code{NULL}, the current R random-number stream is used.}
  \item{verbose}{Show progress messages and the final sampler's iteration
  output? CV sampler output is suppressed.}
}
\details{
Initialization mean-fills missing entries, then applies Lloyd k-means with
25 starts and at most 500 iterations. Initial component variances that are
zero or undefined use the observed-data variance (or 1 if that is unavailable).
If fewer than \code{K} distinct starting values exist, k-means uses the
available centers and the remaining, unoccupied components are initialized
from their priors. This does not change \code{K} in the fitted model.

CV subsets are sampled without replacement from observed entries and are
shared across ranks. Each subset is hidden in its own fit; the other
validation subsets remain observed in that fit. The default therefore hides
5 percent in each of five fits and uses about 25 percent of observed entries
for validation overall. It is not an exhaustive five-fold partition.
Ranks are fitted sequentially within each worker. Validation subsets are
distributed across available workers. Under the same software and numerical
environment, the stored task seeds make results independent of worker count.

For subset \eqn{A} and rank \eqn{r}, the score is
\deqn{\log\{B^{-1}\sum_{b=1}^B \exp(\sum_{(i,j)\in A}
 \log p(x_{ij}\mid\theta^{(b)},r))\}.}
This integrates the joint held-out likelihood by Monte Carlo. It is not
the sum of separately integrated pointwise log predictive densities, nor
the mean of log-likelihood draws. Scores are added across subsets and the
rank with the largest total is selected. The final model is refitted using
all entries observed in the original input. Genuine missing entries are
imputed but never scored against invented truth.

The backend uses zero-based iteration labels, runs from 0 through
\code{M - 1}, and saves when \code{t > burn_in} and \code{t \%\% thin == 0}.
The wrapper's divisibility requirement keeps its allocation and saving loop
consistent without modifying the backend. Draws within a chain can be
autocorrelated. Inspect convergence and assess Monte Carlo error for the
scientific summaries of interest. The current backend does not poll for an
R interrupt inside its iteration loop.

The normal kernels have support on the real line. For bounded data, decide
on a suitable analysis scale before fitting and back-transform each completed
draw afterwards. The fixed variance prior makes units relevant to inference.
}
\value{
An object of class \code{lorps_mi}, with:
\describe{
  \item{\code{X}}{The numeric input matrix, preserving row and column names.}
  \item{\code{missing_index}}{Linear indices of missing cells in R's
  column-major order: the same order as \code{which(is.na(X))}.}
  \item{\code{imputations}}{A matrix with one row per missing entry and
  one column per retained joint draw. With complete input it has zero rows.}
  \item{\code{mu_draws}}{Component-mean draws, with \code{K} rows and the
  same draw columns. These are not all the model parameters.}
  \item{\code{R, K, n_imputations}}{Fitted rank, truncation level, and actual
  number of retained draws.}
  \item{\code{iterations}}{Retained zero-based backend iteration labels.}
  \item{\code{cv}}{\code{NULL} for fixed rank, otherwise a list:
    \code{results} gives one row per rank, fold scores, total and mean log
    LPPL, and a selection flag; \code{fold_scores} is the rank-by-subset
    score matrix; \code{folds} contains the held-out linear indices;
    \code{n_observed} and \code{fold_size} record sizes; \code{controls}
    contains resolved settings; \code{cores_used} is the actual worker
    count; and \code{seeds} contains initialization and sampler task seeds.}
  \item{\code{settings}}{Final chain settings, prior variance, seed,
  rank-selection method, and a prior joint-density truncation upper bound
  capped at the maximum possible L1 distance of 2.}
  \item{\code{seeds, rng_kind, call}}{Final initialization and sampler seeds,
  R RNG configuration, and the call.}
}
No accuracy statistic is computed for genuinely missing values. Complete
datasets can be reconstructed with \code{\link{complete_lorps}}.
}
\references{
Kakoty, B., Chekouo, T., and Lock, E.
\emph{Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices}.
Manuscript, Theorem 1 and Section 3.4.
}
\seealso{\code{\link{complete_lorps}}, \code{\link{lorps_suggest_K}}}
\examples{
# Short run to illustrate the interface, not to establish convergence.
set.seed(10)
X <- matrix(rnorm(48), 8, 6)
X[c(2, 17, 28)] <- NA
fit <- lorps_mi(X, R = 1, K = 2, M = 30, burn_in = 10,
                seed = 12, verbose = FALSE)
fit
dim(fit$imputations)  # 3 missing entries by 19 retained draws
completed <- complete_lorps(fit, 1)
stopifnot(identical(completed[!is.na(X)], X[!is.na(X)]))

\dontrun{
# Defaults: select R from 1:15 using five disjoint 5-percent holdouts.
fit_cv <- lorps_mi(X, seed = 12)
fit_cv$cv$results

# Smaller rank search and a serial run.
fit_cv <- lorps_mi(X, K = 3, cand_rank_list = 1:3,
                   cv_control = list(cores = 1), seed = 12)
}
}
