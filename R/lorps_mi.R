# The public interface keeps validation, initialization, and CV in small helpers.
# The sampler itself is compiled from src/finite_stickbreak_eff.cpp.

lorps_mi <- function(X, R = NULL, K = 20, CV = is.null(R),
                     cand_rank_list = 1:15,
                     M = 2000, burn_in = 1000, thin = 1,
                     mus_init = NULL, sigma2s_init = NULL, C_ij_init = NULL,
                     tau2 = 10000, cv_control = list(),
                     seed = NULL, verbose = TRUE) {
  call <- match.call()
  X <- .lorps_validate_X(X)
  K <- .lorps_integer(K, "K", lower = 2L)
  .lorps_flag(CV, "CV")
  .lorps_flag(verbose, "verbose")
  .lorps_positive(tau2, "tau2")
  chain <- .lorps_chain(M, burn_in, thin)
  initial <- .lorps_validate_initial(X, K, mus_init, sigma2s_init, C_ij_init)
  if (!is.null(R)) R <- .lorps_integer(R, "R")
  if (!is.null(seed)) seed <- .lorps_integer(seed, "seed", lower = 0L)

  if (CV && !is.null(R)) {
    warning("CV = TRUE: selecting R from cand_rank_list; the supplied R is ignored.",
            call. = FALSE)
  }
  if (CV && !is.null(initial)) {
    warning("Custom starting values are used for the final fit only; CV uses training-only k-means.",
            call. = FALSE)
  }
  if (CV) {
    ranks <- .lorps_ranks(cand_rank_list)
    control <- .lorps_cv_control(cv_control, sum(!is.na(X)))
  } else {
    if (length(cv_control)) {
      warning("cv_control is ignored because CV = FALSE.", call. = FALSE)
    }
    if (is.null(R)) R <- 5L
  }

  # An explicit seed is local to this call; set.seed() also works when seed = NULL.
  restore_rng <- .lorps_local_seed(seed)
  if (!is.null(restore_rng)) on.exit(restore_rng(), add = TRUE)
  rng_kind <- RNGkind()
  cv <- NULL
  seeds <- list(final_init = sample.int(.Machine$integer.max, 1L),
                final_sampler = sample.int(.Machine$integer.max, 1L))

  if (CV) {
    cv <- .lorps_cross_validate(X, K, ranks, tau2, control, rng_kind, verbose)
    R <- cv$results$R[which.max(cv$results$total_log_lppl)]
    if (verbose) message("Selected R = ", R, "; refitting all observed entries.")
  }

  .lorps_set_seed(seeds$final_init, rng_kind)
  if (is.null(initial)) initial <- .lorps_initial(X, K, tau2)
  missing_index <- which(is.na(X))
  if (verbose) {
    message("Final fit: R = ", R, ", K = ", K, ", ", chain$n_save,
            " retained imputations; ", length(missing_index), " missing entries.")
  }
  .lorps_set_seed(seeds$final_sampler, rng_kind)
  raw <- .lorps_run_sampler(X, K, R, initial,
                           score_values = rep(-Inf, length(missing_index)),
                           chain = chain, tau2 = tau2, verbose = verbose)
  imputations <- raw[[1L]]
  mu_draws <- raw[[2L]]
  draw_names <- paste0("iteration_", chain$iterations)
  colnames(imputations) <- colnames(mu_draws) <- draw_names
  rownames(mu_draws) <- paste0("component_", seq_len(K))

  # With no validation truth the backend score is identically zero, not a fit metric.
  # It is deliberately not reported as the likelihood of genuinely missing values.
  structure(list(call = call, X = X, missing_index = missing_index,
                 imputations = imputations, mu_draws = mu_draws,
                 R = R, K = K, n_imputations = ncol(imputations),
                 iterations = chain$iterations, cv = cv,
                 settings = list(M = chain$M, burn_in = chain$burn_in,
                                 thin = chain$thin, tau2 = tau2, seed = seed,
                                 rank_selection = if (CV) "LPPL CV" else "fixed",
                                 truncation_L1_bound = min(2, 2 * length(X) * 2^(1 - K))),
                 seeds = seeds, rng_kind = rng_kind),
            class = "lorps_mi")
}
