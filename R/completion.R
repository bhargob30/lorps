complete_lorps <- function(object, action = 1L) {
  if (!inherits(object, "lorps_mi")) stop("object must be a lorps_mi fit.", call. = FALSE)
  if (is.character(action)) {
    if (length(action) != 1L || is.na(action) || !action %in% c("all", "mean")) {
      stop("action must be 'all', 'mean', or integer draw indices (0 returns X).", call. = FALSE)
    }
    if (action == "mean") {
      X <- object$X
      X[object$missing_index] <- rowMeans(object$imputations)
      return(X)
    }
    action <- seq_len(object$n_imputations)
  }
  if (!is.numeric(action) || is.complex(action) || !length(action) ||
      any(!is.finite(action)) || any(action != floor(action)) ||
      any(action < 0 | action > object$n_imputations)) {
    stop("action must contain integer indices from 0 to n_imputations.", call. = FALSE)
  }
  one <- function(draw) {
    X <- object$X
    if (draw > 0L) X[object$missing_index] <- object$imputations[, draw]
    X
  }
  if (length(action) == 1L) return(one(action))
  stats::setNames(lapply(action, one), paste0("imputation_", action))
}

print.lorps_mi <- function(x, ...) {
  cat("LoRPS multiple imputation\n")
  cat("  Matrix:", nrow(x$X), "rows x", ncol(x$X), "columns\n")
  cat("  Missing entries:", length(x$missing_index), "\n")
  cat("  Rank R: ", x$R, " (", x$settings$rank_selection, ")\n", sep = "")
  cat("  Mixture components K:", x$K, "\n")
  cat("  Retained imputations:", x$n_imputations, "\n")
  cat("  Use complete_lorps(x, 1) for one completed matrix.\n")
  invisible(x)
}

lorps_suggest_K <- function(N, Q, tolerance = 1e-3) {
  N <- .lorps_integer(N, "N")
  Q <- .lorps_integer(Q, "Q")
  .lorps_positive(tolerance, "tolerance")
  if (tolerance >= 2) stop("tolerance must be smaller than 2 (the maximum L1 distance).", call. = FALSE)
  as.integer(max(2, 1 + ceiling(1 + log2(N) + log2(Q) - log2(tolerance))))
}
