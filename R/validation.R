.lorps_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max) {
    stop(name, " must be one integer between ", lower, " and ",
         .Machine$integer.max, ".", call. = FALSE)
  }
  as.integer(x)
}

.lorps_positive <- function(x, name) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be one positive finite number.", call. = FALSE)
  }
  invisible(x)
}

.lorps_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

.lorps_validate_X <- function(X) {
  if (!is.matrix(X) || !is.numeric(X) || is.complex(X) ||
      any(dim(X) < 1L) || length(X) > .Machine$integer.max) {
    stop("X must be a nonempty real numeric matrix (rows by columns).", call. = FALSE)
  }
  if (any(is.infinite(X))) stop("X may contain NA, but not Inf or -Inf.", call. = FALSE)
  if (all(is.na(X))) stop("X must contain at least one observed value.", call. = FALSE)
  if (any(rowSums(!is.na(X)) == 0L) || any(colSums(!is.na(X)) == 0L)) {
    warning("X has entirely missing rows or columns; their factors have no direct observations.",
            call. = FALSE)
  }
  storage.mode(X) <- "double"
  X
}

.lorps_chain <- function(M, burn_in, thin) {
  M <- .lorps_integer(M, "M")
  burn_in <- .lorps_integer(burn_in, "burn_in", lower = 0L)
  thin <- .lorps_integer(thin, "thin")
  if (M %% thin != 0L || burn_in %% thin != 0L) {
    stop("M and burn_in must both be multiples of thin for this sampler's saving rule.",
         call. = FALSE)
  }
  n_save <- (M - burn_in) / thin - 1L
  if (n_save < 2L) {
    stop("Choose M, burn_in, and thin to retain at least two draws: (M - burn_in)/thin - 1 >= 2.",
         call. = FALSE)
  }
  list(M = M, burn_in = burn_in, thin = thin, n_save = as.integer(n_save),
       iterations = seq.int(burn_in + thin, M - thin, by = thin))
}

.lorps_validate_initial <- function(X, K, mus, sigma2, C) {
  supplied <- !vapply(list(mus, sigma2, C), is.null, logical(1L))
  if (!any(supplied)) return(NULL)
  if (!all(supplied)) {
    stop("Supply mus_init, sigma2s_init, and C_ij_init together, or leave all three NULL.",
         call. = FALSE)
  }
  if (!is.numeric(mus) || is.complex(mus) || length(mus) != K || any(!is.finite(mus))) {
    stop("mus_init must contain K finite means.", call. = FALSE)
  }
  if (!is.numeric(sigma2) || is.complex(sigma2) || length(sigma2) != K ||
      any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    stop("sigma2s_init must contain K positive finite variances.", call. = FALSE)
  }
  if (!is.matrix(C) || !is.numeric(C) || is.complex(C) || !identical(dim(C), dim(X)) ||
      any(!is.finite(C)) || any(C != floor(C)) || any(C < 0 | C >= K)) {
    stop("C_ij_init must match dim(X) and contain zero-based labels 0, ..., K - 1.",
         call. = FALSE)
  }
  storage.mode(C) <- "integer"
  list(mus = as.numeric(mus), sigma2 = as.numeric(sigma2), C = C)
}

.lorps_ranks <- function(x) {
  if (!is.numeric(x) || is.complex(x) || !length(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < 1 | x > .Machine$integer.max)) {
    stop("cand_rank_list must contain positive integer ranks.", call. = FALSE)
  }
  sort(unique(as.integer(x)))
}

.lorps_cv_control <- function(control, n_observed) {
  defaults <- list(nfold = 5L, holdout_fraction = 0.05,
                   M = 1000L, burn_in = 500L, thin = 1L, cores = NULL)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) ||
       anyNA(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% names(defaults))))) {
    stop("cv_control must be a named list using nfold, holdout_fraction, M, burn_in, thin, or cores.",
         call. = FALSE)
  }
  defaults[names(control)] <- control
  control <- defaults
  control$nfold <- .lorps_integer(control$nfold, "cv_control$nfold", lower = 2L)
  .lorps_positive(control$holdout_fraction, "cv_control$holdout_fraction")
  if (control$holdout_fraction >= 1) {
    stop("cv_control$holdout_fraction must be smaller than one.", call. = FALSE)
  }
  control$chain <- .lorps_chain(control$M, control$burn_in, control$thin)
  control$fold_size <- floor(control$holdout_fraction * n_observed)
  if (control$fold_size < 1L || control$nfold * control$fold_size > n_observed ||
      control$fold_size >= n_observed) {
    stop("Cannot form these disjoint validation subsets. Change nfold or holdout_fraction, or supply R to skip CV.",
         call. = FALSE)
  }
  available <- as.integer(parallelly::availableCores(constraints = "connections-16"))
  if (length(available) != 1L || is.na(available) || available < 1L) available <- 1L
  if (!is.null(control$cores)) {
    requested <- .lorps_integer(control$cores, "cv_control$cores")
  } else {
    requested <- 5L
  }
  control$cores <- min(requested, control$nfold, available)
  control$available_cores <- available
  control
}

.lorps_local_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    do.call(RNGkind, as.list(kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

.lorps_set_seed <- function(seed, kind) {
  do.call(RNGkind, as.list(kind))
  set.seed(seed)
  invisible(NULL)
}
