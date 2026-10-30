.lorps_initial <- function(X, K, tau2) {
  observed <- X[!is.na(X)]
  X_init <- X
  X_init[is.na(X_init)] <- mean(observed)
  x <- as.vector(X_init)
  n_centers <- min(K, length(unique(x)))
  fallback_variance <- stats::var(observed)
  if (!is.finite(fallback_variance) || fallback_variance <= 0) fallback_variance <- 1

  if (n_centers == 1L) {
    centers <- mean(x)
    cluster <- rep.int(1L, length(x))
  } else {
    km <- stats::kmeans(x, centers = n_centers, nstart = 25L,
                       iter.max = 500L, algorithm = "Lloyd")
    centers <- as.numeric(km$centers[, 1L])
    cluster <- km$cluster
  }
  variances <- vapply(seq_len(n_centers), function(k) stats::var(x[cluster == k]),
                     numeric(1L))
  variances[!is.finite(variances) | variances <= 0] <- fallback_variance
  if (n_centers < K) {
    # Unoccupied initial components are legitimate: initialize them from the model prior.
    centers <- c(centers, stats::rnorm(K - n_centers, sd = sqrt(tau2)))
    variances <- c(variances, 1 / stats::rgamma(K - n_centers, shape = 2, rate = 1))
  }
  list(mus = centers, sigma2 = variances,
       C = matrix(as.integer(cluster - 1L), nrow = nrow(X), ncol = ncol(X)))
}

.lorps_run_sampler <- function(X, K, R, initial, score_values, chain, tau2, verbose) {
  run <- function() {
    mcmcDraws_fsbp(X = X, K = K, R = R,
                   mus_init = initial$mus, sigma2s_init = initial$sigma2,
                   C_ij_init = initial$C, all_mis_elem = score_values,
                   M = chain$M, burn_in = chain$burn_in, thin = chain$thin,
                   tau2 = tau2)
  }
  if (verbose) {
    raw <- run()
  } else {
    invisible(utils::capture.output(raw <- run()))
  }
  expected <- c(sum(is.na(X)), chain$n_save)
  if (length(raw) != 3L || !is.matrix(raw[[1L]]) ||
      !identical(as.integer(dim(raw[[1L]])), as.integer(expected)) ||
      !identical(as.integer(dim(raw[[2L]])), c(K, chain$n_save)) ||
      length(raw[[3L]]) != chain$n_save) {
    stop("The sampler returned unexpected dimensions; inspect the compiled sampler version.",
         call. = FALSE)
  }
  if (any(!is.finite(raw[[1L]])) || any(!is.finite(raw[[2L]])) ||
      any(!is.finite(raw[[3L]]))) {
    stop("The sampler returned non-finite values. Check the data scale, starting values, and model settings.",
         call. = FALSE)
  }
  raw
}

.lorps_log_mean_exp <- function(x) {
  if (!length(x) || any(!is.finite(x))) {
    stop("Validation log-likelihood draws must be finite and nonempty.", call. = FALSE)
  }
  maximum <- max(x)
  maximum + log(mean(exp(x - maximum)))
}
