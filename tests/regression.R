# Base-R regression checks for the public wrapper and its safety boundaries.
# These short chains check software contracts, not posterior convergence.
library(lorps)

run_regression_checks <- function() {
  namespace <- asNamespace("lorps")
  internal <- function(name) get(name, envir = namespace, inherits = FALSE)
  assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
  }
  expect_error <- function(expr, pattern = NULL) {
    condition <- tryCatch({ force(expr); NULL }, error = identity)
    assert(inherits(condition, "error"), "Expected an error, but the call succeeded.")
    if (!is.null(pattern)) {
      assert(grepl(pattern, conditionMessage(condition), fixed = TRUE),
             paste("Unexpected error:", conditionMessage(condition)))
    }
    invisible(condition)
  }
  collect_warnings <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(force(expr), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, messages = messages)
  }

  # Leave the caller's random-number settings unchanged, including on failure.
  original_kind <- RNGkind()
  original_had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  original_seed <- if (original_had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(original_kind))
    if (original_had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")

  truth <- matrix(sin(seq_len(48L) * 0.31) + seq_len(48L) / 73,
                  nrow = 8L, ncol = 6L,
                  dimnames = list(paste0("row_", 1:8), paste0("col_", 1:6)))
  missing <- c(2L, 9L, 18L, 29L, 35L, 47L)
  X <- truth
  X[missing] <- NA_real_
  X_original <- X
  observed <- which(!is.na(X))
  chain <- internal(".lorps_chain")(20L, 10L, 1L)

  # The core's strict t > burn_in rule retains nine, not ten, draws here.
  assert(identical(chain$iterations, 11:19) && identical(chain$n_save, 9L),
         "The saved-iteration convention changed.")
  large_chain <- internal(".lorps_chain")(10000L, 5000L, 50L)
  assert(identical(large_chain$n_save, 99L) &&
           identical(large_chain$iterations, seq.int(5050L, 9950L, by = 50L)),
         "The historical final-chain saving convention was misreported.")
  expect_error(internal(".lorps_chain")(101, 50, 10), "multiples of thin")
  expect_error(internal(".lorps_chain")(100, 51, 10), "multiples of thin")
  expect_error(internal(".lorps_chain")(20, 19, 1), "at least two draws")
  expect_error(internal(".lorps_chain")(20, 20, 1), "at least two draws")
  expect_error(internal(".lorps_chain")(20.5, 10, 1), "M must be")
  expect_error(internal(".lorps_chain")(20, -1, 1), "burn_in must be")
  expect_error(internal(".lorps_chain")(20, 10, 0), "thin must be")

  # Fail before entering C++ on structurally unsafe input.
  expect_error(lorps_mi(as.vector(X), R = 1, K = 2), "numeric matrix")
  expect_error(lorps_mi(as.data.frame(X), R = 1, K = 2), "numeric matrix")
  expect_error(lorps_mi(matrix(NA_real_, 2, 2), R = 1, K = 2), "at least one observed")
  X_infinite <- X
  X_infinite[1L] <- Inf
  expect_error(lorps_mi(X_infinite, R = 1, K = 2), "not Inf")
  expect_error(lorps_mi(X, R = 0, K = 2), "R must be")
  expect_error(lorps_mi(X, R = 1, K = 1), "K must be")
  expect_error(lorps_mi(X, R = 1, K = 2, tau2 = 0), "tau2 must be")
  expect_error(lorps_mi(X, R = 1, K = 2, CV = NA), "CV must be")
  expect_error(lorps_mi(X, R = 1, K = 2, mus_init = c(-1, 1)), "together")
  C_valid <- matrix(0L, nrow(X), ncol(X))
  expect_error(lorps_mi(X, R = 1, K = 2, mus_init = c(-1, 1),
                       sigma2s_init = c(1, 0), C_ij_init = C_valid),
               "positive finite variances")
  C_bad <- C_valid
  C_bad[1L] <- 2L
  expect_error(lorps_mi(X, R = 1, K = 2, mus_init = c(-1, 1),
                       sigma2s_init = c(1, 1), C_ij_init = C_bad),
               "zero-based labels")
  expect_error(internal(".lorps_ranks")(c(1, 1.5)), "positive integer ranks")
  expect_error(internal(".lorps_cv_control")(list(unknown = 1), length(observed)),
               "named list")
  expect_error(internal(".lorps_cv_control")(list(holdout_fraction = .001), length(observed)),
               "Cannot form")
  X_empty_row <- X
  X_empty_row[1L, ] <- NA_real_
  empty_row_check <- collect_warnings(internal(".lorps_validate_X")(X_empty_row))
  assert(any(grepl("entirely missing", empty_row_check$messages, fixed = TRUE)),
         "Entirely missing rows or columns should be identified.")

  # Explicit seeds are reproducible and local to the call.
  set.seed(911)
  caller_seed <- .Random.seed
  caller_kind <- RNGkind()
  fit <- lorps_mi(X, R = 1, K = 2, M = 20, burn_in = 10, thin = 1,
                  seed = 314, verbose = FALSE)
  assert(identical(.Random.seed, caller_seed) && identical(RNGkind(), caller_kind),
         "An explicitly seeded fit changed the caller's RNG state.")
  fit_again <- lorps_mi(X, R = 1, K = 2, M = 20, burn_in = 10, thin = 1,
                        seed = 314, verbose = FALSE)
  assert(identical(fit$imputations, fit_again$imputations) &&
           identical(fit$mu_draws, fit_again$mu_draws),
         "Repeated seeded fits are not identical.")
  assert(identical(X, X_original), "The sampler mutated the input matrix.")
  assert(identical(fit$missing_index, missing) &&
           identical(dim(fit$imputations), c(length(missing), 9L)) &&
           identical(dim(fit$mu_draws), c(2L, 9L)) &&
           identical(fit$iterations, 11:19) && identical(fit$n_imputations, 9L),
         "The fit did not retain the expected ordered draws and dimensions.")
  assert(is.null(fit$cv) && fit$R == 1L, "Supplying R should disable CV by default.")
  assert(!any(grepl("lppl|likelihood", names(fit), ignore.case = TRUE)),
         "A genuine-missing-data final fit should not report a likelihood score.")
  assert(all(is.finite(fit$imputations)) && all(is.finite(fit$mu_draws)),
         "A small valid fit returned nonfinite draws.")

  # Reconstruction keeps observed values and matrix names exactly; a column is a draw.
  completed <- complete_lorps(fit, "all")
  assert(is.list(completed) && length(completed) == fit$n_imputations,
         "Completing all draws should return one matrix per imputation.")
  for (i in seq_along(completed)) {
    assert(identical(dimnames(completed[[i]]), dimnames(X)) &&
             identical(completed[[i]][observed], X[observed]) &&
             identical(unname(completed[[i]][missing]), unname(fit$imputations[, i])),
           "Completion changed an observed entry, a name, or the draw order.")
  }
  assert(identical(complete_lorps(fit, 0), X), "Action 0 must return the original incomplete X.")
  selected <- complete_lorps(fit, c(3L, 1L))
  assert(identical(selected[[1L]], completed[[3L]]) &&
           identical(selected[[2L]], completed[[1L]]),
         "Requested nonconsecutive draw indices were reordered.")
  mean_completed <- complete_lorps(fit, "mean")
  assert(identical(mean_completed[observed], X[observed]) &&
           identical(unname(mean_completed[missing]), unname(rowMeans(fit$imputations))),
         "Mean completion does not use the missing-cell posterior means.")
  expect_error(complete_lorps(fit, 10L), "integer indices")
  expect_error(complete_lorps(fit, -1L), "integer indices")
  expect_error(complete_lorps(fit, 1.5), "integer indices")
  expect_error(complete_lorps(fit, NA_real_), "integer indices")
  expect_error(complete_lorps(fit, "median"), "action must be")
  expect_error(complete_lorps(list(), 1), "lorps_mi fit")

  # No preexisting .Random.seed must remain no seed after an explicitly seeded call.
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(lorps_mi(X, R = 1, K = 2, M = 20, burn_in = 10, thin = 1,
                     seed = 314, verbose = FALSE))
  assert(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
         "An explicit seed left .Random.seed behind when it did not exist before.")
  # set.seed() remains usable when the wrapper's seed argument is omitted.
  set.seed(432)
  from_global_seed <- lorps_mi(X, R = 1, K = 2, M = 20, burn_in = 10, verbose = FALSE)
  set.seed(432)
  from_global_seed_again <- lorps_mi(X, R = 1, K = 2, M = 20, burn_in = 10, verbose = FALSE)
  assert(identical(from_global_seed$imputations, from_global_seed_again$imputations),
         "set.seed() did not reproduce a fit with seed = NULL.")

  # Directly exercise the unchanged core's score sentinel. Scores must not feed back
  # into fitting: only S_Pt changes when held-out truths change at the same seed.
  set.seed(77)
  initial <- internal(".lorps_initial")(X, 2L, 10000)
  run_scored <- function(values) {
    set.seed(78)
    internal(".lorps_run_sampler")(X, 2L, 1L, initial, values, chain, 10000, FALSE)
  }
  score_values <- rep(-Inf, length(missing))
  score_values[c(2L, 5L)] <- truth[missing[c(2L, 5L)]]
  raw_none <- run_scored(rep(-Inf, length(missing)))
  raw_scored <- run_scored(score_values)
  changed_score_values <- score_values
  changed_score_values[c(2L, 5L)] <- changed_score_values[c(2L, 5L)] + 4
  raw_changed <- run_scored(changed_score_values)
  assert(all(raw_none[[3L]] == 0), "Genuinely unknown values contributed a backend score.")
  assert(identical(raw_none[[1L]], raw_scored[[1L]]) &&
           identical(raw_scored[[1L]], raw_changed[[1L]]) &&
           identical(raw_none[[2L]], raw_scored[[2L]]) &&
           identical(raw_scored[[2L]], raw_changed[[2L]]),
         "Held-out truths changed the fitted draws instead of only their scores.")
  assert(all(is.finite(raw_scored[[3L]])) &&
           !identical(raw_scored[[3L]], raw_changed[[3L]]),
         "Finite held-out truths were not reflected in likelihood scores.")

  # K is the fitted model truncation, even when the data contain fewer unique values.
  constant_X <- matrix(0.4, 4L, 3L)
  constant_X[c(2L, 9L)] <- NA_real_
  set.seed(93)
  padded_initial <- internal(".lorps_initial")(constant_X, 20L, 10000)
  assert(length(padded_initial$mus) == 20L && length(padded_initial$sigma2) == 20L &&
           all(is.finite(padded_initial$mus)) &&
           all(is.finite(padded_initial$sigma2) & padded_initial$sigma2 > 0) &&
           identical(unique(as.vector(padded_initial$C)), 0L),
         "Unoccupied starting components did not preserve K or valid initial values.")
  padded_fit <- lorps_mi(constant_X, R = 1, K = 20, M = 20, burn_in = 10,
                         seed = 94, verbose = FALSE)
  assert(padded_fit$K == 20L && nrow(padded_fit$mu_draws) == 20L &&
           all(is.finite(padded_fit$imputations)),
         "Fitting a low-uniqueness matrix silently reduced K or failed.")
  tiny_complete <- matrix(c(-1, 0, 1, 2), 2L, 2L)
  set.seed(95)
  tiny_initial <- internal(".lorps_initial")(tiny_complete, 20L, 10000)
  assert(length(tiny_initial$mus) == 20L && length(tiny_initial$sigma2) == 20L &&
           length(unique(as.vector(tiny_initial$C))) == 4L &&
           all(is.finite(tiny_initial$mus)) &&
           all(is.finite(tiny_initial$sigma2) & tiny_initial$sigma2 > 0),
         "Four distinct cells could not initialize four occupied and sixteen empty components.")
  assert(identical(lorps_suggest_K(348, 574, tolerance = 1e-3), 30L),
         "The known BRCA truncation recommendation changed.")
  recommended_K <- lorps_suggest_K(8, 6, tolerance = .01)
  bound <- function(K) 2 * 8 * 6 * 2^(1 - K)
  assert(bound(recommended_K) <= .01 && bound(recommended_K - 1L) > .01,
         "The truncation helper did not return the minimal sufficient integer K.")

  # Five disjoint 5% observed-cell holdouts; a final refit retains every original
  # observed value. CV is deliberately serial here so checks need no socket access.
  cv_control <- list(M = 16L, burn_in = 8L, thin = 2L, cores = 1L)
  fit_cv <- lorps_mi(X, K = 2, cand_rank_list = c(2L, 1L, 2L),
                     M = 20, burn_in = 10, cv_control = cv_control,
                     seed = 560, verbose = FALSE)
  fold_indices <- unlist(fit_cv$cv$folds, use.names = FALSE)
  assert(length(fit_cv$cv$folds) == 5L &&
           all(lengths(fit_cv$cv$folds) == floor(.05 * length(observed))) &&
           !anyDuplicated(fold_indices) && all(fold_indices %in% observed),
         "CV folds overlap, include genuine missing cells, or have the wrong size.")
  assert(identical(fit_cv$missing_index, missing) && identical(fit_cv$X, X) &&
           identical(complete_lorps(fit_cv, 1)[observed], X[observed]),
         "CV cells remained hidden in the final refit.")
  assert(identical(as.integer(fit_cv$cv$results$R), 1:2) &&
           fit_cv$R == fit_cv$cv$results$R[which.max(fit_cv$cv$results$total_log_lppl)],
         "Rank candidates or rank selection are inconsistent with saved CV scores.")
  assert(all(is.finite(fit_cv$cv$fold_scores)) &&
           isTRUE(all.equal(unname(rowSums(fit_cv$cv$fold_scores)),
                            unname(fit_cv$cv$results$total_log_lppl))),
         "Saved CV totals do not match their fold scores.")

  # Complete matrices can still select rank by held-out scoring. The final fit
  # has zero missing rows in its draw matrix and completion returns X unchanged.
  complete_input_fit <- lorps_mi(truth, K = 2, cand_rank_list = 1:2,
                                 M = 20, burn_in = 10, cv_control = cv_control,
                                 seed = 561, verbose = FALSE)
  assert(identical(dim(complete_input_fit$imputations), c(0L, 9L)) &&
           identical(complete_input_fit$missing_index, integer()) &&
           identical(complete_lorps(complete_input_fit, 1L), truth) &&
           identical(complete_lorps(complete_input_fit, "mean"), truth),
         "Complete input failed the zero-missing-entry fit or reconstruction contract.")
  complete_folds <- unlist(complete_input_fit$cv$folds, use.names = FALSE)
  assert(length(complete_folds) == 10L && !anyDuplicated(complete_folds) &&
           all(complete_folds %in% seq_along(truth)) &&
           all(is.finite(complete_input_fit$cv$fold_scores)) &&
           any(complete_input_fit$cv$fold_scores != 0),
         "Complete input did not produce valid held-out CV likelihood scores.")

  # Custom starts are restricted to the final fit. They must neither determine
  # validation masks nor enter the CV chains or their scores.
  custom_cv <- collect_warnings(lorps_mi(
    X, R = 9, CV = TRUE, K = 2, cand_rank_list = c(1L, 2L),
    M = 20, burn_in = 10, cv_control = cv_control,
    mus_init = c(-7, 7), sigma2s_init = c(2, 3), C_ij_init = C_valid,
    seed = 560, verbose = FALSE
  ))
  assert(any(grepl("final fit only", custom_cv$messages, fixed = TRUE)) &&
           any(grepl("supplied R is ignored", custom_cv$messages, fixed = TRUE)),
         "The custom-initialization and rank-override rules were not made explicit.")
  assert(identical(fit_cv$cv$folds, custom_cv$value$cv$folds) &&
           identical(fit_cv$cv$fold_scores, custom_cv$value$cv$fold_scores),
         "Custom final-fit initialization leaked into CV or reproducible fold scores.")

  # Test the fold boundary independently: change only held-out truths while
  # recording the training matrix and initialized parameters passed to the sampler.
  fold_function <- internal(".lorps_cv_fold")
  recording_env <- new.env(parent = namespace)
  recorded <- list()
  real_initializer <- internal(".lorps_initial")
  recording_env$.lorps_initial <- function(X, K, tau2) real_initializer(X, K, tau2)
  recording_env$.lorps_run_sampler <- function(X, K, R, initial, score_values,
                                               chain, tau2, verbose) {
    recorded[[length(recorded) + 1L]] <<- list(X = X, initial = initial,
                                             score_values = score_values)
    list(matrix(0, sum(is.na(X)), chain$n_save),
         matrix(0, K, chain$n_save), rep(-1, chain$n_save))
  }
  environment(fold_function) <- recording_env
  heldout <- observed[c(3L, 8L)]
  changed_data <- X
  changed_data[heldout] <- changed_data[heldout] + 1000
  fold_args <- list(fold_id = 1L, K = 2L, ranks = 1L, tau2 = 10000,
                    control = list(chain = chain), folds = list(heldout),
                    init_seeds = 701L, sampler_seeds = matrix(702L, 1L, 1L),
                    rng_kind = RNGkind())
  invisible(do.call(fold_function, c(list(data = X), fold_args)))
  invisible(do.call(fold_function, c(list(data = changed_data), fold_args)))
  assert(length(recorded) == 2L && identical(recorded[[1L]]$X, recorded[[2L]]$X) &&
           identical(recorded[[1L]]$initial, recorded[[2L]]$initial),
         "Held-out truths entered fold training data or initialization.")
  fold_missing <- which(is.na(recorded[[1L]]$X))
  assert(all(is.infinite(recorded[[1L]]$score_values[match(missing, fold_missing)])) &&
           all(recorded[[1L]]$score_values[match(missing, fold_missing)] < 0) &&
           identical(unname(recorded[[1L]]$score_values[match(heldout, fold_missing)]),
                     unname(X[heldout])),
         "The fold score payload failed to separate genuine NA from validation truth.")

  fallback_rank <- lorps_mi(X, CV = FALSE, K = 2, M = 20, burn_in = 10,
                            seed = 888, verbose = FALSE)
  assert(fallback_rank$R == 5L && is.null(fallback_rank$cv),
         "CV = FALSE with R omitted should use rank five.")
  cat("lorps base-R regression checks passed.\n")
}

run_regression_checks()
