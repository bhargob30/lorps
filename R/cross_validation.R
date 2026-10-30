.lorps_cv_fold <- function(fold_id, data, K, ranks, tau2, control,
                           folds, init_seeds, sampler_seeds, rng_kind) {
  X <- data
  heldout <- folds[[fold_id]]
  truth <- X[heldout]
  X_cv <- X
  X_cv[heldout] <- NA_real_
  missing_index <- which(is.na(X_cv))
  score_values <- rep(-Inf, length(missing_index))
  score_values[match(heldout, missing_index)] <- truth
  .lorps_set_seed(init_seeds[fold_id], rng_kind)
  initial <- .lorps_initial(X_cv, K, tau2)
  scores <- numeric(length(ranks))
  for (r in seq_along(ranks)) {
    .lorps_set_seed(sampler_seeds[r, fold_id], rng_kind)
    raw <- tryCatch(
      .lorps_run_sampler(X_cv, K, ranks[r], initial, score_values,
                        control$chain, tau2, verbose = FALSE),
      error = function(e) stop("CV subset ", fold_id, ", R = ", ranks[r], ": ",
                               conditionMessage(e), call. = FALSE)
    )
    scores[r] <- .lorps_log_mean_exp(as.numeric(raw[[3L]]))
    # Keep only the fold score, not every CV imputation matrix.
    rm(raw)
  }
  scores
}

.lorps_start_cluster <- function(cores) {
  # Each worker runs one serial sampler. Set limits before starting worker R sessions.
  thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
  old <- Sys.getenv(thread_vars, unset = NA_character_)
  on.exit({
    Sys.unsetenv(thread_vars[is.na(old)])
    if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
  }, add = TRUE)
  one <- stats::setNames(rep("1", length(thread_vars)), thread_vars)
  do.call(Sys.setenv, as.list(one))
  cl <- parallel::makePSOCKcluster(cores)
  ready <- FALSE
  on.exit(if (!ready) parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(paths) {
    .libPaths(paths)
    loadNamespace("lorps")
    NULL
  }, .libPaths())
  ready <- TRUE
  cl
}

.lorps_cross_validate <- function(X, K, ranks, tau2, control, rng_kind, verbose) {
  observed <- which(!is.na(X))
  chosen <- observed[sample.int(length(observed), control$nfold * control$fold_size)]
  folds <- split(chosen, rep(seq_len(control$nfold), each = control$fold_size))
  folds <- lapply(folds, sort)
  init_seeds <- sample.int(.Machine$integer.max, control$nfold)
  sampler_seeds <- matrix(sample.int(.Machine$integer.max,
                                    length(ranks) * control$nfold),
                          nrow = length(ranks), ncol = control$nfold)
  if (verbose) {
    message("LPPL CV: ", length(ranks), " ranks x ", control$nfold,
            " disjoint subsets of ", control$fold_size, " observed entries.")
    message("CV workers: ", control$cores,
            ". Each worker evaluates ranks in sequence; this can take time.")
  }
  args <- list(data = X, K = K, ranks = ranks, tau2 = tau2, control = control,
               folds = folds, init_seeds = init_seeds,
               sampler_seeds = sampler_seeds, rng_kind = rng_kind)
  workers <- control$cores
  cl <- NULL
  if (workers > 1L) {
    cl <- tryCatch(.lorps_start_cluster(workers), error = function(e) {
      warning("Parallel CV could not start; using serial CV: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  }
  if (is.null(cl)) {
    workers <- 1L
    result <- do.call(lapply, c(list(X = seq_len(control$nfold),
                                    FUN = .lorps_cv_fold), args))
  } else {
    on.exit(parallel::stopCluster(cl), add = TRUE)
    result <- do.call(parallel::parLapply,
                      c(list(cl = cl, X = seq_len(control$nfold),
                             fun = .lorps_cv_fold), args))
  }
  fold_scores <- do.call(cbind, result)
  rownames(fold_scores) <- paste0("R_", ranks)
  colnames(fold_scores) <- paste0("fold_", seq_len(control$nfold))
  total <- rowSums(fold_scores)
  table <- data.frame(R = ranks, fold_scores, total_log_lppl = total,
                      mean_log_lppl = rowMeans(fold_scores),
                      selected = seq_along(ranks) == which.max(total),
                      row.names = NULL)
  list(results = table, fold_scores = fold_scores, folds = folds,
       n_observed = length(observed), fold_size = control$fold_size,
       controls = control, cores_used = workers,
       seeds = list(initialization = init_seeds, sampler = sampler_seeds))
}
