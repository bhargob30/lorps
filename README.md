
# lorps

**Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices**

`lorps` implements Bayesian multiple imputation using a Low-Rank Probit Stick-Breaking Process (LoRPS). The method models matrix entries through flexible normal mixtures with entry-specific weights, while a low-rank tensor structure allows information to be shared across rows and columns. Missing entries are sampled from their posterior predictive distributions to generate multiple completed datasets.

## Installation

```r
install.packages("remotes")
remotes::install_github("bhargob30/lorps")
```

## Usage

The following example generates a matrix with an underlying low-rank signal, randomly masks 20% of its entries, and evaluates imputation accuracy using normalized mean squared error (NMSE).

```r
library(lorps)

# Generate a complete matrix with a low-rank signal.
set.seed(10)

N <- 30
Q <- 20
signal_rank <- 2

U <- matrix(rnorm(N * signal_rank), nrow = N)
V <- matrix(rnorm(Q * signal_rank), nrow = Q)

signal <- tcrossprod(U, V)
X_true <- signal + matrix(rnorm(N * Q, sd = 0.3), nrow = N)

# Randomly mask 20% of the entries.
missing_index <- sort(
  sample.int(length(X_true), size = floor(0.20 * length(X_true)))
)

X <- X_true
X[missing_index] <- NA_real_

# Fit LoRPS with a fixed tensor rank.
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

# Retrieve completed datasets.
X1 <- complete_lorps(fit, action = 1)
X2 <- complete_lorps(fit, action = 2)

# Posterior mean imputation.
X_mean <- complete_lorps(fit, action = "mean")

# Evaluate imputation accuracy at the masked entries.
true_missing <- X_true[fit$missing_index]
imputed_mean <- X_mean[fit$missing_index]

nmse <- mean((imputed_mean - true_missing)^2) /
  var(true_missing)

print(c(NMSE = nmse))
```

The NMSE compares the mean imputed values with the true values at the artificially masked entries, normalized by the sample variance of those true values. Smaller values indicate better imputation accuracy. For genuinely missing data, NMSE cannot be computed without known reference values.

### Automatic rank selection

The tensor rank `R` can also be selected by maximizing the cross-validated log posterior predictive likelihood (LPPL):

```r
fit_cv <- lorps_mi(
  X,
  K = 10,
  cand_rank_list = 1:3,
  cv_control = list(cores = 2),
  M = 2000,
  burn_in = 1000,
  thin = 5,
  seed = 12,
  verbose = FALSE
)

fit_cv$R
fit_cv$cv$results
```

Here, `R` controls the rank of the latent mixture-weight tensor, whereas `K` is the number of normal mixture components. The helper function `lorps_suggest_K()` provides guidance for choosing `K` based on a prior truncation bound.

The examples illustrate the package interface; MCMC convergence should be assessed before using imputations for inference.

For additional options and documentation, see `?lorps_mi`, `?complete_lorps`, and `?lorps_suggest_K`.

## Reference

Kakoty, B., Chekouo, T., and Lock, E. F. *Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices*. Manuscript.
