\name{lorps_mi}
\alias{lorps_mi}
\title{Multiple Imputation with a Low-Rank Probit Stick-Breaking Model}
\description{
Performs Bayesian multiple imputation for incomplete numeric matrices
using a Low-Rank Probit Stick-Breaking Process (LoRPS). Entry-specific
normal mixtures accommodate flexible distributions, while a low-rank
latent structure shares information across rows and columns. The tensor
rank can be fixed or selected using held-out predictive likelihood.
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

  \item{R}{A positive integer tensor rank parameter, or \code{NULL}.
  With the defaults, \code{NULL} requests cross-validation (CV), and
  an explicit value requests a fixed-rank fit.
  If \code{CV = FALSE} and \code{R = NULL}, rank 5 is used.
  Explicit \code{CV = TRUE} overrides a supplied \code{R}, with a warning.}

  \item{K}{Number of normal mixture components in the finite stick-breaking
  representation, at least 2. The default is 20. The final component
  receives the remaining stick mass. See \code{\link{lorps_suggest_K}} for
  guidance based on a prior truncation bound.}

  \item{CV}{Logical. Whether to select \code{R} by maximizing the
  cross-validated log posterior predictive likelihood (LPPL).
  The default is \code{is.null(R)}.}

  \item{cand_rank_list}{Positive integer candidate ranks. Used only for CV;
  duplicate values are removed and the values are sorted. Larger lists
  require more sampler fits. Exact ties select the smaller rank.}

  \item{M, burn_in, thin}{Total sampler iterations, burn-in setting, and
  thinning interval. \code{M} and \code{thin} must be positive integers;
  \code{burn_in} must be a nonnegative integer. Both \code{M} and
  \code{burn_in} must be multiples of \code{thin}. Under the sampler's
  saving convention, \code{(M - burn_in) / thin - 1} draws are retained,
  and this number must be at least 2. The defaults give 999 draws.
  These settings do not guarantee convergence.}

  \item{mus_init, sigma2s_init, C_ij_init}{Optional starting values for
  the component means, variances, and assignments. Supply all three
  together or leave all three as \code{NULL}. By default, k-means-based
  initialization is used, with missing entries temporarily mean-filled
  to construct starting values. The first two arguments must be vectors
  of length \code{K}, containing finite means and positive finite
  variances. \code{C_ij_init} must be a matrix with \code{dim(X)} and
  zero-based integer labels \code{0:(K - 1)}, including positions missing
  in \code{X}. Custom starting values apply only to the final fit;
  CV fits are initialized using their own training data.}

  \item{tau2}{Positive finite variance of the normal prior
  \eqn{\mu_k \sim \mathcal{N}(0,\tau^2)}{mu[k] ~ Normal(0, tau2)}
  on each component mean. The default is 10000. The remaining prior
  distributions are fixed and described under Model.}

  \item{cv_control}{A named list overriding any CV defaults:
  \describe{
    \item{\code{nfold = 5}}{Number of disjoint validation subsets, at least 2.}

    \item{\code{holdout_fraction = 0.05}}{Fraction of originally observed
    entries in each validation subset, rounded down to an integer count.}

    \item{\code{M = 1000, burn_in = 500, thin = 1}}{CV sampler settings,
    giving 499 retained draws by default. The same divisibility and
    retained-draw checks apply as for the final fit.}

    \item{\code{cores = NULL}}{Automatic use of up to five workers,
    capped by the number of subsets and \code{parallelly::availableCores()}.
    A positive integer supplies a different upper limit; 1 forces serial
    CV. A PSOCK cluster is used across platforms. If cluster startup
    fails, CV falls back to serial execution with a warning.}
  }}

  \item{seed}{Optional nonnegative integer. An explicit seed is local to
  this call, and the previous R random-number state is restored on exit.
  With \code{NULL}, the current R random-number stream is used.}

  \item{verbose}{Logical. Whether to show progress messages and the final
  sampler's iteration output. CV sampler output is suppressed.}
}
\details{
\subsection{Model}{
Let \eqn{\Theta}{Theta} denote the model parameters. Conditional on these
parameters, the matrix entries are independent, with
\deqn{
X_{ij} \mid \Theta \sim
\sum_{k=1}^{K} w_{ijk}\,\mathcal{N}(\mu_k,\sigma_k^2).
}{
X[i,j] | Theta ~ sum over k of w[i,j,k] Normal(mu[k], sigma2[k]).
}
The component means and variances are shared across the matrix,
whereas the mixture weights vary by row and column. The weights follow
\eqn{w_{ijk}=V_{ijk}\prod_{\ell<k}(1-V_{ij\ell})}{w[i,j,k] = V[i,j,k] times the product of (1 - V[i,j,l]) for l < k},
where \eqn{V_{ijk}=\Phi(Z_{ijk})}{V[i,j,k] = Phi(Z[i,j,k])} for
\eqn{k<K}{k < K}, and \eqn{V_{ijK}=1}{V[i,j,K] = 1}.
Here, \eqn{\Phi}{Phi} is the standard normal cumulative distribution
function. The latent tensor has a CP/PARAFAC representation with
\eqn{R} terms:
\deqn{
Z_{ijk}=\sum_{r=1}^{R}Z_1[i,r]Z_2[j,r]Z_3[k,r].
}{
Z[i,j,k] = sum over r of Z1[i,r] Z2[j,r] Z3[k,r].
}
The factor matrices link rows, columns, and mixture components.
Low-rank structure is therefore imposed on the latent variables governing
the weights, rather than directly on the data matrix.

Factor entries have independent standard normal priors. Independently,
component means have normal priors with mean zero and variance
\code{tau2}, and component variances have inverse-gamma priors with
fixed shape 2 and scale 1.
}

\subsection{Mixture truncation: K}{
\code{K} is the number of normal components used to approximate the
infinite mixture. The helper \code{\link{lorps_suggest_K}()} provides
guidance using a bound on the discrepancy between the finite and infinite
prior joint marginal densities. Its result must be supplied explicitly
as \code{K}; \code{lorps_mi()} does not select \code{K} automatically.
}

\subsection{Tensor rank: R}{
\code{R} controls the complexity of the latent factorization. When
\code{CV = TRUE}, the model is fitted at each candidate rank with
selected observed entries temporarily held out. For each validation
set, its joint held-out likelihood is averaged over posterior draws,
and the logarithms of these averages are summed across validation sets.
This gives the cross-validated log posterior predictive likelihood (LPPL) score. The rank with the largest
score is selected, and the model is refitted using all originally
observed entries.

By default, five disjoint validation sets each contain approximately
5\% of the observed entries; they do not partition all observed entries
into five equal folds. The same sets are used for every candidate rank,
and genuinely missing entries are never scored. A larger rank search
requires additional fits; a fixed rank avoids this search.
}

\subsection{Multiple imputations}{
Posterior computation uses Markov chain Monte Carlo (MCMC). Each
retained iteration supplies one joint draw of all missing entries,
which, together with the unchanged observed entries, forms a completed
dataset. The draws are stored compactly in \code{imputations}; completed
matrices are constructed using \code{\link{complete_lorps}()}.
The number retained is reported in \code{n_imputations} and depends
on \code{M}, \code{burn_in}, and \code{thin}.
}
}
\value{
An object of class \code{lorps_mi}, with:
\describe{
  \item{\code{X}}{The numeric input matrix, preserving row and column names.}

  \item{\code{missing_index}}{Linear indices of missing entries in R's
  column-major order, matching \code{which(is.na(X))}.}

  \item{\code{imputations}}{A matrix with one row per missing entry and
  one column per retained joint draw. With complete input it has zero rows.}

  \item{\code{mu_draws}}{Component-mean draws, with \code{K} rows and
  columns corresponding to the same retained iterations as
  \code{imputations}. Other model-parameter draws are not returned.}

  \item{\code{R, K, n_imputations}}{Selected or supplied rank, truncation
  level, and actual number of retained draws.}

  \item{\code{iterations}}{Retained zero-based sampler iteration labels.}

  \item{\code{cv}}{\code{NULL} for a fixed-rank fit; otherwise, a list:
    \code{results} gives one row per rank, validation scores, their total
    and mean, and a selection flag; \code{fold_scores} is the rank-by-subset
    score matrix; \code{folds} contains the held-out linear indices;
    \code{n_observed} and \code{fold_size} record sizes; \code{controls}
    contains resolved settings; \code{cores_used} is the actual worker
    count; and \code{seeds} contains initialization and sampler task seeds.}

  \item{\code{settings}}{Final sampler settings, prior variance, seed,
  rank-selection method, and a prior joint-density truncation upper bound
  capped at the maximum possible L1 distance of 2.}

  \item{\code{seeds, rng_kind, call}}{Final initialization and sampler seeds,
  R random-number generator configuration, and the function call.}
}
No accuracy statistic is computed for genuinely missing values. Completed
datasets can be reconstructed with \code{\link{complete_lorps}}.
}
\references{
Kakoty, B., Chekouo, T., and Lock, E. F.
\emph{Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices}.
Manuscript.
}
\seealso{\code{\link{complete_lorps}}, \code{\link{lorps_suggest_K}}}
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

fit
dim(fit$imputations)

# Average the imputation draws for each missing entry.
imputed_mean <- rowMeans(fit$imputations)

# Match the true values to the ordering of the stored imputations.
true_missing <- X_true[fit$missing_index]

# NMSE: mean squared prediction error divided by
# the sample variance of the true masked values.
nmse <- mean((imputed_mean - true_missing)^2) / var(true_missing)

print(c(NMSE = nmse))

# Obtain one completed dataset or a posterior-mean completion.
completed <- complete_lorps(fit, action = 1)
completed_mean <- complete_lorps(fit, action = "mean")

# Observed entries remain unchanged.
stopifnot(identical(completed[!is.na(X)], X[!is.na(X)]))

\dontrun{
# Alternatively, select the tensor rank by predictive-likelihood CV.
fit_cv <- lorps_mi(
  X,
  K = 10,
  cand_rank_list = 1:3,
  cv_control = list(cores = 1),
  M = 2000,
  burn_in = 1000,
  thin = 5,
  seed = 12,
  verbose = FALSE
)

fit_cv$R
fit_cv$cv$results

imputed_mean_cv <- rowMeans(fit_cv$imputations)
true_missing_cv <- X_true[fit_cv$missing_index]

nmse_cv <- mean((imputed_mean_cv - true_missing_cv)^2) /
  var(true_missing_cv)

print(c(NMSE = nmse_cv))
}
}