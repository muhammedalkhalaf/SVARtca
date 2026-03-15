#' @title Interface with VAR Estimation
#' @description Convenience function to run TCA directly from
#'   estimated VAR model objects.
#' @name tca_var_interface
NULL

#' Run TCA from a VAR Estimation Object
#'
#' Extracts coefficient matrices and residual covariance from a fitted
#' VAR model (from the \pkg{vars} package) and runs TCA.
#'
#' @param var_model A fitted VAR model object from \code{\link[vars]{VAR}}.
#' @param from Shock variable (integer or name).
#' @param intermediates Integer vector or character vector of intermediate
#'   variable names.
#' @param h Maximum horizon (default: 20).
#' @param order Transmission ordering (default: variable ordering in the VAR).
#' @param mode Decomposition mode: \code{"overlapping"},
#'   \code{"exhaustive_3way"}, or \code{"exhaustive_4way"}.
#' @param identification Identification scheme: \code{"cholesky"} (default)
#'   or \code{"manual"}.
#' @param Phi0 Manual impact matrix. Required if \code{identification = "manual"}.
#' @return A \code{tca_result} object (see \code{\link{tca_analyze}}).
#' @export
#' @examples
#' \dontrun{
#' library(vars)
#' data(Canada)
#' var_est <- VAR(Canada, p = 2, type = "const")
#' result <- tca_from_var(var_est, from = "e",
#'                         intermediates = c("prod", "rw"),
#'                         h = 20, mode = "exhaustive_4way")
#' plot_tca(result, target = "U")
#' }
tca_from_var <- function(var_model, from, intermediates, h = 20,
                          order = NULL, mode = "overlapping",
                          identification = "cholesky", Phi0 = NULL) {
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("Package 'vars' is required. Install with install.packages('vars').")
  }

  K <- var_model$K
  p <- var_model$p
  var_names <- colnames(var_model$y)

  # Extract coefficient matrices A1, A2, ..., Ap
  coef_list <- vars::Acoef(var_model)
  As <- lapply(seq_len(p), function(i) coef_list[[i]])

  # Residual covariance
  Sigma <- summary(var_model)$covres

  # Identification
  if (identification == "cholesky") {
    Phi0 <- t(chol(Sigma))
  } else if (identification == "manual") {
    if (is.null(Phi0)) stop("Phi0 must be provided for manual identification.")
    if (!all(dim(Phi0) == c(K, K))) stop("Phi0 must be K x K.")
  } else {
    stop("identification must be 'cholesky' or 'manual'.")
  }

  # Default ordering
  if (is.null(order)) order <- seq_len(K)

  # Resolve character names to indices
  if (is.character(from)) {
    from <- match(from, var_names)
    if (is.na(from)) stop("'from' variable not found in VAR model.")
  }
  if (is.character(intermediates)) {
    intermediates <- match(intermediates, var_names)
    if (any(is.na(intermediates))) {
      stop("One or more 'intermediates' not found in VAR model.")
    }
  }
  if (is.character(order)) {
    order <- match(order, var_names)
    if (any(is.na(order))) stop("Order contains names not found in VAR model.")
  }

  # Systems form
  sf <- tca_systems_form(Phi0, As, h = h, order = order)

  # Run TCA
  tca_analyze(from = from, B = sf$B, Omega = sf$Omega,
              intermediates = intermediates, K = K, h = h,
              order = order, mode = mode, var_names = var_names)
}
