# A small illustration of the model and the lorps API.
# The short chain demonstrates usage; it is not a convergence assessment.

library(lorps)

set.seed(121)
N <- 30L
Q <- 20L
K_true <- 3L
R_true <- 2L

# Latent row, column, and stick factors.
U <- matrix(rnorm(N * R_true), nrow = N)
V <- matrix(rnorm(Q * R_true), nrow = Q)
B <- matrix(rnorm((K_true - 1L) * R_true), nrow = K_true - 1L)

# Probit stick-breaking weights, with the last component taking the remainder.
weights <- array(0, dim = c(N, Q, K_true))
remaining <- matrix(1, nrow = N, ncol = Q)
for (k in seq_len(K_true - 1L)) {
  eta <- tcrossprod(sweep(U, 2L, B[k, ], FUN = "*"), V)
  stick <- pnorm(eta)
  weights[, , k] <- remaining * stick
  remaining <- remaining * (1 - stick)
}
weights[, , K_true] <- remaining

component_means <- c(-2, 0, 2)
component_sds <- c(0.35, 0.5, 0.4)
X_complete <- matrix(NA_real_, nrow = N, ncol = Q,
                      dimnames = list(paste0("sample", seq_len(N)),
                                      paste0("feature", seq_len(Q))))
for (j in seq_len(Q)) {
  for (i in seq_len(N)) {
    component <- sample.int(K_true, size = 1L, prob = weights[i, j, ])
    X_complete[i, j] <- rnorm(1L, component_means[component], component_sds[component])
  }
}

# Mask 20% of entries. These reference values are used only for evaluation.
heldout_index <- sort(sample.int(length(X_complete), floor(0.20 * length(X_complete))))
X <- X_complete
X[heldout_index] <- NA_real_

# Fixed rank and the small known generating truncation keep this example light.
# K is a modeling choice; for a prior-density tolerance use lorps_suggest_K().
fit <- lorps_mi(
  X,
  R = R_true,
  K = K_true,
  M = 600L,
  burn_in = 300L,
  thin = 3L,
  seed = 2026L,
  verbose = FALSE
)

completed_1 <- complete_lorps(fit, action = 1L)
completed_mean <- complete_lorps(fit, action = "mean")
completed_first_three <- complete_lorps(fit, action = 1:3)

# Missing entries are ordered by fit$missing_index, in R column-major order.
truth <- X_complete[fit$missing_index]
prediction <- completed_mean[fit$missing_index]
nmse <- mean((prediction - truth)^2) / var(truth)

# Entry-level empirical reconstruction intervals from the saved MI draws.
limits <- apply(fit$imputations, 1L, quantile,
                 probs = c(0.025, 0.975), type = 1L, names = FALSE)
coverage95 <- mean(truth >= limits[1L, ] & truth <= limits[2L, ])
mean_width95 <- mean(limits[2L, ] - limits[1L, ])

cat("Retained imputations:", fit$n_imputations, "\n")
cat("Fitted rank:", fit$R, "  Truncation:", fit$K, "\n")
print(c(NMSE = nmse, coverage95 = coverage95, mean_interval_width95 = mean_width95))

# These are illustrative results on artificially hidden entries. Genuine NA
# entries have no known reference truth, so their NMSE and coverage are unknown.
# Reconstruction intervals do not by themselves include population sampling
# uncertainty for a later scientific estimand.

# Optional: automatically select the rank using validation holdouts.
# This adds candidate-rank fits and is intentionally not run by this example.
# fit_cv <- lorps_mi(
#   X,
#   K = K_true,
#   cand_rank_list = 1:4,
#   M = 1000L,
#   burn_in = 500L,
#   thin = 5L,
#   cv_control = list(nfold = 5L, holdout_fraction = 0.05, cores = 1L),
#   seed = 2027L
# )
# fit_cv$R
# fit_cv$cv$results

# For the manuscript's joint prior-density bound on this matrix:
# lorps_suggest_K(N, Q, tolerance = 0.001)
