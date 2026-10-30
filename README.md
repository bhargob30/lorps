# lorps

**Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices**

`lorps` implements a Bayesian nonparametric multiple imputation method based on a low-rank probit stick-breaking process. The method accommodates flexible data distributions while borrowing information across rows and columns of a matrix to impute missing values and quantify uncertainty.

## Installation

```r
install.packages("remotes")
remotes::install_github("bhargob30/lorps")
```

## Usage

```r
library(lorps)

# Example matrix with missing entries
set.seed(123)
X <- matrix(rnorm(300), nrow = 20)
X[sample(length(X), 60)] <- NA

# Multiple imputation with automatic rank selection
fit <- lorps_mi(X)

# Alternatively, specify the rank
fit <- lorps_mi(X, R = 3)

# Retrieve completed datasets
X1 <- complete_lorps(fit, 1)
X2 <- complete_lorps(fit, 2)

# Posterior mean imputation
X_mean <- complete_lorps(fit, "mean")
```

For additional options and documentation, see `?lorps_mi`.

## Reference

Kakoty, B., Chekouo, T., and Lock, E. *Bayesian Nonparametric Multiple Imputation for High-Dimensional Matrices*. Manuscript.
