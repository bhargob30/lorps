# lorps

**Multiple imputation for numeric matrices using a low-rank probit stick-breaking mixture.**

`lorps` models each matrix entry with a normal mixture. Its mixture weights share information across rows and columns through a low-rank latent representation. The package saves draws of missing entries and constructs completed datasets when requested.

The main function is `lorps_mi()`. Supply a numeric matrix with missing entries represented by `NA`. You can fix the tensor rank or select it using held-out predictive likelihood. No values for genuinely missing entries are needed.

## Installation

The package requires R 4.1 or later and a C++17 compiler. Source installation uses `Rcpp`, `RcppArmadillo` (at least 14.6.0-1), `RcppTN`, and `parallelly`; the installer obtains the R dependencies. Windows source installation requires an appropriate [Rtools](https://cran.r-project.org/bin/windows/Rtools/) toolchain.

```r
install.packages("remotes")
remotes::install_github(bhargob30/lorps", upgrade = "never")
```

For local installation and repository setup, see [GITHUB_SETUP.md](GITHUB_SETUP.md).

## A first fit

```r
library(lorps)

set.seed(12)
X <- matrix(rnorm(30 * 20), nrow = 30, ncol = 20)
X[sample.int(length(X), size = 120)] <- NA_real_

# Providing R uses a fixed rank by default. K defaults to 20.
# This short chain illustrates the interface.
fit <- lorps_mi(
  X,
  R = 2,
  M = 600,
  burn_in = 300,
  thin = 3,
  seed = 2026,
  verbose = FALSE
)

fit$R
fit$K
fit$n_imputations

completed_1 <- complete_lorps(fit, action = 1)
completed_mean <- complete_lorps(fit, action = "mean")
completed_first_three <- complete_lorps(fit, action = 1:3)
```

The completed matrices retain the observed entries. `action = "mean"` replaces each missing entry with its average across saved draws. It is a point summary; use separate completed datasets when you need to propagate imputation uncertainty into a subsequent analysis.

## Fixed rank or automatic rank selection

`R` is the tensor rank. It is different from `K`, the mixture truncation level.

| Call | Rank behavior |
|---|---|
| `lorps_mi(X)` | Selects a rank from `1:15` using validation holdouts. |
| `lorps_mi(X, R = 3)` | Fits rank 3 without rank selection. |
| `lorps_mi(X, CV = FALSE)` | Fits the fallback fixed rank 5. |
| `lorps_mi(X, R = 3, CV = TRUE)` | Warns that the supplied rank is ignored and selects from `cand_rank_list`. |

For example, this uses a smaller candidate set and caps validation parallelism at two workers:

```r
fit_cv <- lorps_mi(
  X,
  K = 20,
  cand_rank_list = 1:5,
  cv_control = list(cores = 2),
  seed = 2026
)

fit_cv$R
fit_cv$cv$results
fit_cv$cv$cores_used
```

Rank selection fits several models before the final fit. It can therefore take substantially longer than using a fixed rank.

### What the validation procedure does

By default, the observed entries supply **five disjoint validation sets, each containing approximately 5% of the observed entries**. Each validation fit hides one set and keeps the other observed entries available. Approximately 25% of observed entries are used for validation across the five sets; the remaining entries stay in every training set. This follows the manuscript's validation design, rather than a conventional five-fold partition of all observed entries.

The same validation sets are used for every candidate rank. Genuine `NA` entries remain missing and are never scored. Once a rank is selected, the model is fitted again using all originally observed entries.

The default validation controls are:

```r
list(
  nfold = 5,
  holdout_fraction = 0.05,
  M = 1000,
  burn_in = 500,
  thin = 1,
  cores = NULL
)
```

When `cores = NULL`, the worker count is limited by five workers, the number of validation sets, and `parallelly::availableCores()`. Supply `cores = 1` for serial execution, or a positive integer to cap the requested workers. Parallel execution uses PSOCK workers on Windows, Linux, and macOS. Validation sets run in parallel; candidate ranks run sequentially within each worker. Each worker holds a separate model fit, so reducing the worker count also reduces memory use.

### The predictive score

For validation set A and posterior draw t, let

$$
\ell_t(A)=\sum_{(i,j)\in A}\log p(X_{ij}\mid\Theta^{(t)}).
$$

The validation score is

$$
\log\left\{\frac{1}{T}\sum_{t=1}^T\exp\bigl(\ell_t(A)\bigr)\right\}.
$$

The implementation evaluates this expression with a stable log-mean-exp calculation. Scores are added across validation sets, and the rank with the largest total is selected.

This is a **joint held-out predictive-density score**. It is different from adding separately averaged marginal log predictive densities for individual cells. The saved validation results report the criterion used for selection.

## Choosing the number of mixture components

`K` includes the final component, which receives the remaining stick mass. The package default is `K = 20`, a practical starting value also used in the manuscript's simulation study. It is not a universal accuracy guarantee.

Theorem 1 in the manuscript gives the prior joint-density bound

$$
\lVert m^{(K)}-m^{(\infty)}\rVert_1
\leq 2NQ\,2^{-(K-1)}.
$$

For a requested tolerance, `lorps_suggest_K()` returns

$$
K=\max\left\{2,\;1+\left\lceil\log_2\left(\frac{2NQ}{\text{tolerance}}\right)\right\rceil\right\}.
$$

```r
lorps_suggest_K(N = 348, Q = 574, tolerance = 0.001)
# 30

K_for_X <- lorps_suggest_K(nrow(X), ncol(X), tolerance = 0.001)
```

For a 348 by 574 matrix, `K = 20` gives a bound of approximately 0.762, while `K = 30` gives approximately 0.000744. The manuscript uses `K = 30` for that application.

The bound compares **prior joint marginal densities under the finite and infinite models**. It does not bound posterior imputation error, determine the optimal number of occupied components, or guarantee a particular predictive score. The `0.001` tolerance is a numerical choice. The theorem assumes the centered latent-factor priors and terminal-stick construction used by the model. The recorded `fit$settings$truncation_L1_bound` is capped at 2, the largest possible L1 distance between probability densities.

## Sampling controls and saved draws

`M` is the total number of sampler iterations, not the requested number of imputations. The sampler uses zero-based iteration indices and saves an iteration when its index is strictly greater than `burn_in` and divisible by `thin`.

The wrapper requires both `M` and `burn_in` to be divisible by `thin`. Under this convention, the number saved is

```text
n_imputations = (M - burn_in) / thin - 1
```

At least two saved draws are required. In particular:

| Controls | Saved draws |
|---|---:|
| Default final fit: `M = 2000`, `burn_in = 1000`, `thin = 1` | 999 |
| Default validation fit: `M = 1000`, `burn_in = 500`, `thin = 1` | 499 |
| Short example: `M = 600`, `burn_in = 300`, `thin = 3` | 99 |

The actual retained source indices are recorded in `fit$iterations`. These are MCMC draws and need not be independent. A short example chain demonstrates the interface; it does not establish convergence or adequate sampling for another dataset.

## Initialization, priors, and input scale

With no custom initialization, missing entries are mean-filled temporarily for initialization and the resulting values are clustered. If there are fewer distinct values than `K`, the available centers are used and the remaining initial means and variances are drawn from their priors. The fitted model still has `K` components. Mean filling is only an initialization step; missing entries are subsequently sampled by the model.

To supply initial values, provide all three arguments together:

| Argument | Required shape and interpretation |
|---|---|
| `mus_init` | Numeric vector of length `K`, containing initial component means. |
| `sigma2s_init` | Positive numeric vector of length `K`, containing initial variances. |
| `C_ij_init` | Matrix with the same dimensions as `X`, containing **zero-based** component labels `0:(K - 1)` for every entry. |

Custom initial values apply only to the final fit. Every validation fit receives fresh initialization based on its own training data.

The component-mean prior is normal with mean zero and variance `tau2`, whose default is 10000. Component variances have an inverse-gamma IG(2, 1) prior. The latent-factor entries have independent standard-normal priors. Only `tau2` is exposed among these prior settings.

The package fits the scale of the matrix you supply. It does not automatically standardize, transform, or impose value bounds. For example, if you fit logit-transformed methylation values, transform each completed dataset back with `plogis()` before calculating quantities on the methylation scale. Observed proportions must be strictly between zero and one before using `qlogis()`. In general, transforming the mean of draws differs from averaging their transformed values.

Supply `seed` to reproduce the validation and sampler random streams. With an explicit seed, the call preserves the caller's outside RNG state. Used seeds are retained in the result.

## Returned object

`lorps_mi()` returns an object of class `lorps_mi` with these main components:

| Component | Contents |
|---|---|
| `X` | Original matrix, including its `NA` entries. |
| `missing_index` | Indices of the missing entries in R's column-major order. |
| `imputations` | Matrix with one row per missing entry and one column per retained joint draw. |
| `mu_draws` | Saved component-mean draws. |
| `R`, `K` | Fitted tensor rank and mixture truncation level. |
| `n_imputations` | Actual number of saved imputations. |
| `iterations` | Retained zero-based sampler iteration indices. |
| `cv` | Validation results when selection is used: `results`, `fold_scores`, `folds`, `n_observed`, `fold_size`, `controls`, `cores_used`, and `seeds`. |
| `settings` | Sampling and prior settings, rank-selection information, and the truncation bound. |
| `seeds` | Recorded random seeds. |
| `call` | Original function call. |

An imputation column is a coherent draw of all missing entries. Keeping each column together preserves the dependence represented by that draw.

`complete_lorps()` provides these forms:

| `action` | Return value |
|---|---|
| `0` | Original incomplete matrix. |
| One positive integer | One completed matrix. |
| Several positive integers | A list of the corresponding completed matrices. |
| `"all"` | A list of all completed matrices. |
| `"mean"` | One matrix with mean imputations. |

Creating every completed matrix can use much more memory than retaining the compact `imputations` matrix.

## Downstream multiple-imputation analysis

Completed matrices can be converted to the ordinary data-frame lists used
by established MI tools. For the first example, whose unnamed columns become
`V1`, `V2`, and so on:

```r
install.packages("mitools")  # Optional; not needed for imputation itself.
completed <- lapply(complete_lorps(fit, "all"), as.data.frame)
mi_data <- mitools::imputationList(completed)
analyses <- with(mi_data, lm(V1 ~ V2 + V3))
mitools::MIcombine(analyses)
```

Use the model and column names appropriate to the scientific question.
Pooling complete-data estimates also needs their within-imputation
variances; the variation between imputations alone is not the full
population sampling uncertainty. The package does not automatically select
or validate a downstream model.

When complete matrices are large, avoid making the entire list at once:

```r
analyses <- lapply(seq_len(fit$n_imputations), function(m) {
  dat <- as.data.frame(complete_lorps(fit, m))
  lm(V1 ~ V2 + V3, data = dat)
})
mitools::MIcombine(analyses)
```

Save and reload the compact result with `saveRDS(fit, "lorps_fit.rds")`
and `readRDS("lorps_fit.rds")`. No refitting is needed to reconstruct those
same completed datasets or compute new summaries from them.

## Simulation example and help

A complete example simulates data from a small finite low-rank probit stick-breaking mixture, masks 20% of the entries, and evaluates reconstruction of those hidden values:

```r
source(system.file("examples", "simulated_example.R", package = "lorps"))
```

The example reports entry-level normalized mean squared error and empirical interval coverage because the artificially hidden values are known. These accuracy quantities cannot be computed for genuinely missing entries without a reference truth. Automatic rank selection is included as an optional commented block.

```r
?lorps_mi
?complete_lorps
?lorps_suggest_K
citation("lorps")
```

## Citation and license

Bhargob Kakoty, Thierry Chekouo, and Eric Lock. *Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices*. Unpublished manuscript.

The package is distributed under GPL-3.

Before an inferential release, review the augmented-state scan question and
other preserved-core properties in [development notes](inst/DEVELOPMENT_NOTES.md).
