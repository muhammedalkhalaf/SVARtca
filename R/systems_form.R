#' @title Systems Form Construction Functions
#' @description Core functions for building the systems form representation
#'   x = Bx + Omega*epsilon following Wegner, Lieb, Smeekes (2025).
#' @references
#'   Wegner, E., Lieb, L., Smeekes, S. (2025).
#'   \emph{Transmission Channel Analysis in Dynamic Models.}
#'   arXiv:2405.18987.
#'   \url{https://github.com/enweg/tca-matlab-toolbox}
#' @name systems_form
NULL

#' LD Decomposition (Cholesky-inverse method)
#'
#' Computes the LD decomposition of a symmetric positive-definite matrix
#' following the MATLAB TCA toolbox: \code{L = inv(chol(Sigma, 'lower'))},
#' \code{D = diag(1 / diag(L))}. This ensures that \code{D \%*\% L} has
#' unit diagonal, making \code{I - D \%*\% L} strictly lower-triangular.
#' Equivalent to \code{makeLD.m} in the MATLAB toolbox.
#'
#' @param Sigma Symmetric positive-definite matrix.
#' @return A list with components \code{L} and \code{D}.
#' @keywords internal
makeLD <- function(Sigma) {
  Linv <- t(chol(Sigma))   # lower-triangular Cholesky factor
  L <- solve(Linv)          # L = inv(chol(Sigma))
  D <- diag(1 / diag(L))   # D rescales so D*L has unit diagonal
  list(L = L, D = D)
}

#' Permutation Matrix
#'
#' Creates a permutation matrix from an ordering vector.
#' Equivalent to \code{permmatrix.m} in the MATLAB toolbox.
#'
#' @param order Integer vector specifying the permutation.
#' @return Permutation matrix.
#' @keywords internal
permmatrix <- function(order) {
  P <- diag(length(order))
  P[order, , drop = FALSE]
}

#' Moving Average Coefficients from VAR
#'
#' Computes MA coefficients C_1, ..., C_h from VAR coefficient matrices
#' using the recursion C_h = sum_{j=1}^{min(h,p)} A_j C_{h-j}, with C_0 = I.
#'
#' @param As List of VAR coefficient matrices (A_1, A_2, ..., A_p).
#' @param h  Maximum horizon.
#' @return List of MA coefficient matrices (C_1, ..., C_h).
#' @keywords internal
computeMA <- function(As, h) {
  K <- nrow(As[[1]])
  p <- length(As)
  Cs <- list()

  for (j in seq_len(h)) {
    C_j <- matrix(0, K, K)
    for (i in seq_len(min(j, p))) {
      if (i == j) {
        C_j <- C_j + As[[i]]
      } else {
        C_j <- C_j + As[[i]] %*% Cs[[j - i]]
      }
    }
    Cs[[j]] <- C_j
  }
  Cs
}

#' Slide-In Block Construction
#'
#' Inserts a row block into a block lower-triangular matrix.
#' Equivalent to \code{slideIn.m} in the MATLAB toolbox.
#'
#' @param B Square matrix.
#' @param A Row block matrix (K x ncols).
#' @return Modified matrix B.
#' @keywords internal
slideIn <- function(B, A) {
  K <- nrow(A)
  nBlocks <- nrow(B) %/% K
  nA_cols <- ncol(A)

  for (i in seq_len(nBlocks)) {
    n_cols <- min(i * K, nA_cols)
    block <- A[, (nA_cols - n_cols + 1):nA_cols, drop = FALSE]
    r_start <- (i - 1) * K + 1
    r_end   <- i * K
    c_end   <- i * K
    c_start <- c_end - n_cols + 1
    B[r_start:r_end, c_start:c_end] <- block
  }
  B
}

#' Build Systems Form B Matrix
#'
#' Constructs the B matrix of dimension K*(h+1) x K*(h+1) for the
#' systems form x = Bx + Omega*epsilon.
#' Equivalent to \code{makeB.m} in the MATLAB toolbox.
#'
#' @param As    List of VAR coefficient matrices.
#' @param Sigma Residual covariance matrix.
#' @param order Transmission ordering vector.
#' @param h     Maximum horizon.
#' @return Systems form B matrix.
#' @keywords internal
makeB_systems <- function(As, Sigma, order, h) {
  K <- nrow(Sigma)
  TT <- permmatrix(order)
  As_star <- lapply(As, function(A) TT %*% A %*% t(TT))
  Sigma_star <- TT %*% Sigma %*% t(TT)
  ld <- makeLD(Sigma_star)
  DL <- ld$D %*% ld$L
  As_trans <- lapply(As_star, function(A) DL %*% A)
  As_flipped <- rev(As_trans)
  rowBlock <- do.call(cbind, As_flipped)
  rowBlock <- cbind(rowBlock, diag(K) - DL)
  dim <- K * (h + 1)
  B <- matrix(0, dim, dim)
  slideIn(B, rowBlock)
}

#' Build Systems Form Omega Matrix
#'
#' Constructs the Omega matrix for the systems form.
#' Equivalent to \code{makeOmega.m} in the MATLAB toolbox.
#'
#' @param Phi0  Structural impact matrix.
#' @param Psis  List of reduced-form MA coefficients.
#' @param order Transmission ordering vector.
#' @param h     Maximum horizon.
#' @return Systems form Omega matrix.
#' @keywords internal
makeOmega_systems <- function(Phi0, Psis, order, h) {
  Sigma <- Phi0 %*% t(Phi0)
  K <- nrow(Sigma)
  TT <- permmatrix(order)
  Sigma_star <- TT %*% Sigma %*% t(TT)
  ld <- makeLD(Sigma_star)
  DL <- ld$D %*% ld$L
  Qt <- ld$L %*% TT %*% Phi0
  Psis_trans <- lapply(Psis, function(Psi) DL %*% TT %*% Psi %*% Phi0)
  Psis_flipped <- rev(Psis_trans)
  if (length(Psis_flipped) > 0) {
    rowBlock <- do.call(cbind, Psis_flipped)
    rowBlock <- cbind(rowBlock, ld$D %*% Qt)
  } else {
    rowBlock <- ld$D %*% Qt
  }
  dim <- K * (h + 1)
  Omega <- matrix(0, dim, dim)
  slideIn(Omega, rowBlock)
}

#' Build Complete Systems Form
#'
#' Constructs both B and Omega matrices for the systems form
#' representation x = Bx + Omega*epsilon.
#'
#' @param Phi0  Structural impact matrix (K x K). For Cholesky
#'   identification, use \code{t(chol(Sigma))}.
#' @param As    List of VAR coefficient matrices (A_1, A_2, ..., A_p).
#' @param h     Maximum IRF horizon.
#' @param order Transmission ordering vector (default: \code{1:K}).
#' @return A list with components:
#'   \describe{
#'     \item{B}{Systems form B matrix (K*(h+1) x K*(h+1)).}
#'     \item{Omega}{Systems form Omega matrix (K*(h+1) x K*(h+1)).}
#'   }
#' @export
#' @examples
#' # 2-variable VAR(1)
#' Phi0 <- matrix(c(1, 0.3, 0, 0.95), 2, 2)
#' As <- list(matrix(c(0.5, -0.1, 0.2, 0.4), 2, 2))
#' sf <- tca_systems_form(Phi0, As, h = 10)
#' dim(sf$B)  # 22 x 22
tca_systems_form <- function(Phi0, As, h, order = NULL, Psis = NULL) {
  K <- nrow(Phi0)
  if (is.null(order)) order <- seq_len(K)
  Sigma <- Phi0 %*% t(Phi0)
  # For pure VAR/SVAR models: Psis should be empty (list()).
  # The MA dynamics emerge from (I-B)^{-1}. Only DSGE/VARMA models
  # pass non-empty Psis. This matches MATLAB's SVAR.m and Recursive.m
  # which pass cell(0,0) (empty) for Psis.
  if (is.null(Psis)) Psis <- list()
  B     <- makeB_systems(As, Sigma, order, h)
  Omega <- makeOmega_systems(Phi0, Psis, order, h)
  list(B = B, Omega = Omega)
}
