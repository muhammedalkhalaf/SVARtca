#' @title Transmission Effect and Path Control Functions
#' @description Functions for computing transmission effects by controlling
#'   paths through the systems form DAG.
#' @name transmission
NULL

#' Map Variable Index to Systems Form Index
#'
#' Converts an original variable index and time period to the
#' corresponding index in the systems form vector.
#' Equivalent to \code{mapY2XInt.m} in the MATLAB toolbox.
#'
#' @param var_idx Variable number in original ordering (1-based).
#' @param time    Time period (0 = contemporaneous).
#' @param K       Number of variables.
#' @param order   Transmission ordering vector.
#' @return 1-based index in the systems form.
#' @keywords internal
mapY2X <- function(var_idx, time, K, order) {
  pos <- which(order == var_idx)
  K * time + pos
}

#' Apply AND Condition (Force Paths Through Variable)
#'
#' Zeros entries in B and Omega that bypass a variable, forcing
#' all paths to pass through it.
#' Equivalent to \code{applyAndToB.m} in the MATLAB toolbox.
#'
#' @param B     Systems form B matrix.
#' @param Omega Systems form Omega matrix.
#' @param from  Shock column index.
#' @param var   Variable index in systems form.
#' @return List with modified \code{B} and \code{Omega}.
#' @keywords internal
applyAndToB <- function(B, Omega, from, var) {
  n <- nrow(B)
  if (var < n) {
    Omega[(var + 1):n, from] <- 0
    if (var > 1) {
      B[(var + 1):n, 1:(var - 1)] <- 0
    }
  }
  list(B = B, Omega = Omega)
}

#' Apply NOT Condition (Block Paths Through Variable)
#'
#' Zeros entries in B and Omega leading into a variable, blocking
#' all paths through it.
#' Equivalent to \code{applyNotToB.m} in the MATLAB toolbox.
#'
#' @param B     Systems form B matrix.
#' @param Omega Systems form Omega matrix.
#' @param from  Shock column index.
#' @param var   Variable index in systems form.
#' @return List with modified \code{B} and \code{Omega}.
#' @keywords internal
applyNotToB <- function(B, Omega, from, var) {
  Omega[var, from] <- 0
  if (var > 1) {
    B[var, 1:(var - 1)] <- 0
  }
  list(B = B, Omega = Omega)
}

#' Compute Transmission Effect
#'
#' Computes the transmission effect by applying AND/NOT conditions
#' to the systems form matrices and solving the resulting linear system.
#' Equivalent to \code{transmissionBOmega.m} in the MATLAB toolbox.
#'
#' @param from     Shock index (1-based, in systems form ordering).
#' @param B        Systems form B matrix.
#' @param Omega    Systems form Omega matrix.
#' @param and_vars Integer vector of systems form indices to force
#'   paths through (AND conditions).
#' @param not_vars Integer vector of systems form indices to block
#'   paths through (NOT conditions).
#' @return Numeric vector of effects (length K*(h+1)).
#' @keywords internal
transmissionEffect <- function(from, B, Omega,
                                and_vars = integer(0),
                                not_vars = integer(0)) {
  B_mod <- B
  O_mod <- Omega

  for (v in and_vars) {
    result <- applyAndToB(B_mod, O_mod, from, v)
    B_mod <- result$B
    O_mod <- result$Omega
  }
  for (v in not_vars) {
    result <- applyNotToB(B_mod, O_mod, from, v)
    B_mod <- result$B
    O_mod <- result$Omega
  }

  n <- nrow(B_mod)
  effects <- solve(diag(n) - B_mod, O_mod[, from])

  if (length(and_vars) > 0) {
    max_and <- max(and_vars)
    effects[1:max_and] <- 0
  }
  effects
}

#' Build NOT Variable Indices for All Time Periods
#'
#' Creates a vector of systems form indices for blocking specified
#' variables across all time periods 0, ..., h.
#'
#' @param var_indices Integer vector of original variable numbers.
#' @param K           Number of variables.
#' @param h           Maximum horizon.
#' @param order       Transmission ordering.
#' @return Integer vector of systems form indices.
#' @keywords internal
not_vars_for <- function(var_indices, K, h, order) {
  nv <- integer(0)
  for (v in var_indices) {
    for (t in 0:h) {
      nv <- c(nv, mapY2X(v, t, K, order))
    }
  }
  nv
}

#' Convert Systems Form Vector to IRF Matrix
#'
#' Transforms a systems form effect vector back to the original
#' variable ordering as an (h+1) x K matrix.
#'
#' @param vec   Numeric vector of length K*(h+1).
#' @param K     Number of variables.
#' @param h     Maximum horizon.
#' @param order Transmission ordering.
#' @return Matrix of dimension (h+1) x K.
#' @keywords internal
vec_to_irf <- function(vec, K, h, order) {
  TT_inv <- t(permmatrix(order))
  mat <- matrix(0, h + 1, K)
  for (t in 0:h) {
    idx <- (t * K + 1):((t + 1) * K)
    mat[t + 1, ] <- as.numeric(TT_inv %*% vec[idx])
  }
  mat
}
