#' @title Main TCA Analysis Functions
#' @description Primary user-facing functions for running Transmission
#'   Channel Analysis with different decomposition modes.
#' @name tca_main
NULL

#' Transmission Channel Analysis
#'
#' Decomposes impulse response functions into transmission channel
#' contributions using the methodology of Wegner, Lieb, Smeekes (2025).
#'
#' Three decomposition modes are supported:
#' \describe{
#'   \item{\code{"overlapping"}}{Each channel is through(j) = total -
#'     not_through(j). Channels may overlap, so their sum may differ
#'     from the total.}
#'   \item{\code{"exhaustive_3way"}}{(2 intermediates only) Non-overlapping:
#'     (1) through var1 inclusive, (2) through var2 only,
#'     (3) direct. Sum equals total.}
#'   \item{\code{"exhaustive_4way"}}{(2 intermediates only) Full
#'     inclusion-exclusion: (1) var1 only, (2) var2 only,
#'     (3) both, (4) direct. Sum equals total.}
#' }
#'
#' @param from          Shock variable number (1-based).
#' @param B             Systems form B matrix (from \code{\link{tca_systems_form}}).
#' @param Omega         Systems form Omega matrix.
#' @param intermediates Integer vector of intermediate variable numbers
#'   (1-based, original ordering).
#' @param K             Number of variables.
#' @param h             Maximum horizon.
#' @param order         Transmission ordering vector.
#' @param mode          Decomposition mode: \code{"overlapping"},
#'   \code{"exhaustive_3way"}, or \code{"exhaustive_4way"}.
#' @param var_names     Character vector of variable names (optional).
#' @return A list of class \code{"tca_result"} with components:
#'   \describe{
#'     \item{irf_total}{Matrix (h+1) x K of total IRFs.}
#'     \item{irf_channels}{Named list of channel IRF matrices, each (h+1) x K.}
#'     \item{channel_names}{Character vector of channel names.}
#'     \item{mode}{Decomposition mode used.}
#'     \item{from}{Shock variable number.}
#'     \item{K}{Number of variables.}
#'     \item{h}{Maximum horizon.}
#'     \item{order}{Transmission ordering.}
#'     \item{var_names}{Variable names.}
#'   }
#' @export
#' @examples
#' # Monetary policy model
#' K <- 4
#' A1 <- matrix(c(0.7,-0.1,0.05,-0.05, -0.3,0.6,0.10,-0.10,
#'                 -0.2,0.1,0.70,0.05, -0.1,0.2,0.05,0.65), K, K, byrow=TRUE)
#' Sigma <- matrix(c(1,0.3,0.2,0.1, 0.3,1.5,0.25,0.15,
#'                    0.2,0.25,0.8,0.1, 0.1,0.15,0.1,0.6), K, K, byrow=TRUE)
#' Phi0 <- t(chol(Sigma))
#' sf <- tca_systems_form(Phi0, list(A1), h = 20)
#' result <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
#'                        intermediates = c(2, 4), K = K, h = 20,
#'                        order = 1:K, mode = "exhaustive_4way",
#'                        var_names = c("IntRate","GDP","Inflation","Wages"))
#' print(result)
tca_analyze <- function(from, B, Omega, intermediates, K, h, order,
                         mode = "overlapping", var_names = NULL) {
  if (is.null(var_names)) var_names <- paste0("Var", seq_len(K))

  mode <- match.arg(mode, c("overlapping", "exhaustive_3way", "exhaustive_4way"))

  if (mode %in% c("exhaustive_3way", "exhaustive_4way") && length(intermediates) != 2) {
    stop(mode, " requires exactly 2 intermediate variables.")
  }

     # Total effect
   total_vec <- transmissionEffect(from, B, Omega)
   total_mat <- vec_to_irf(total_vec, K, h, order)

     # not_through for each intermediate
   nt_vecs <- list()
   th_vecs <- list()
   for (v in intermediates) {
     nv <- not_vars_for(v, K, h, order)
     nt_vecs[[as.character(v)]] <- transmissionEffect(from, B, Omega, not_vars = nv)
     th_vecs[[as.character(v)]] <- total_vec - nt_vecs[[as.character(v)]]
  }

     # not_through all intermediates together
   nv_all <- not_vars_for(intermediates, K, h, order)
   nt_all_vec <- transmissionEffect(from, B, Omega, not_vars = nv_all)

  channel_results <- list()
  channel_names <- character(0)

  if (mode == "overlapping") {
    for (v in intermediates) {
      nm <- paste0("Through ", var_names[v])
      channel_results[[nm]] <- vec_to_irf(th_vecs[[as.character(v)]], K, h, order)
      channel_names <- c(channel_names, nm)
    }
    channel_results[["Direct"]] <- vec_to_irf(nt_all_vec, K, h, order)
    channel_names <- c(channel_names, "Direct")

  } else if (mode == "exhaustive_3way") {
    v1 <- intermediates[1]; v2 <- intermediates[2]
    nm1 <- paste0("Through ", var_names[v1], " (incl.)")
    channel_results[[nm1]] <- vec_to_irf(th_vecs[[as.character(v1)]], K, h, order)
    nm2 <- paste0("Through ", var_names[v2], " only")
    channel_results[[nm2]] <- vec_to_irf(nt_vecs[[as.character(v1)]] - nt_all_vec, K, h, order)
    channel_results[["Direct"]] <- vec_to_irf(nt_all_vec, K, h, order)
    channel_names <- c(nm1, nm2, "Direct")

  } else if (mode == "exhaustive_4way") {
    v1 <- intermediates[1]; v2 <- intermediates[2]
    th_v1 <- th_vecs[[as.character(v1)]]
    th_v2 <- th_vecs[[as.character(v2)]]
    th_or  <- total_vec - nt_all_vec
    th_and <- th_v1 + th_v2 - th_or 
    ch_v1_only <- th_v1 - th_and
    ch_v2_only <- th_v2 - th_and

    nm1 <- paste0(var_names[v1], " only")
    nm2 <- paste0(var_names[v2], " only")
    nm3 <- paste0(var_names[v1], " & ", var_names[v2])
    channel_results[[nm1]] <- vec_to_irf(ch_v1_only, K, h, order)
    channel_results[[nm2]] <- vec_to_irf(ch_v2_only, K, h, order)
    channel_results[[nm3]] <- vec_to_irf(th_and, K, h, order)
    channel_results[["Direct"]] <- vec_to_irf(nt_all_vec, K, h, order)
    channel_names <- c(nm1, nm2, nm3, "Direct")
  }

  result <- list(
    irf_total      = total_mat,
    irf_channels  = channel_results,
    channel_names = channel_names,
    mode          = mode,
    from          = from,
    K             = K,
    h             = h,
    order         = order,
    var_names     = var_names
  )
  class(result) <- "tca_result"
  result
}


#' Binary Decomposition: Total = Through + Not-Through
#'
#' Decomposes the total IRF into the effect passing through a variable
#' and the effect not passing through it. This is an exact decomposition
#' (residual is zero at machine precision).
#'
#' @param from   Shock variable (1-based).
#' @param B      Systems form B matrix.
#' @param Omega   Systems form Omega matrix.
#' @param var_idx Variable to decompose through (1-based).
#' @param K       Number of variables.
#' @param h       Maximum horizon.
#' @param order   Transmission ordering.
#' @return A list with matrices \code{total}, \code{"through},
#'   \code{not_through} (each (h+1) x K).
#' @export
tca_decompose_binary <- function(from, B, Omega, var_idx, K, h, order) {
  total_vec <- transmissionEffect(from, B, Omega)
  nv <- not_vars_for(var_idx, K, h, order)
  nt_vec <- transmissionEffect(from, B, Omega, not_vars = nv)
  th_vec <- total_vec - nt_vec

  list(
    total       = vec_to_irf(total_vec, K, h, order),
    through     = vec_to_irf(th_vec, K, h, order),
    not_through = vec_to_irf(nt_vec, K, h, order)
  )
}

#' Validate Binary Additivity
#'
#' Tests that total = through(j) + not_through(j) holds for all
#' variables at machine precision. This is a diagnostic to verify
#' correct implementation.
#'
#' @param from      Shock variable (1-based).
#' @param B         Systems form B matrix.
#' @param Omega     Systems form Omega matrix.
#' @param K         Number of variables.
#' @param h         Maximum horizon.
#' @param order     Transmission ordering.
#' @param var_names Character vector of variable names (optional).
#' @param verbose   Logical; print results? Default \code{TRUE}.
#' @return Invisibly returns \code{TRUE} if all tests pass.
#' @export
tca_validate_additivity <- function(from, B, Omega, K, h, order,
                                     var_names = NULL, verbose = TRUE) {
  if (is.null(var_names)) var_names <- paste0("Var", seq_len(K))
  total_vec <- transmissionEffect(from, B, Omega)
  max_resid <- 0

  if (verbose) cat("===== Binary Additivity Test =====\n")
  for (v in seq_len(K)) {
    nv <- not_vars_for(v, K, h, order)
    nt <- transmissionEffect(from, B, Omega, not_vars = nv)
    th <- total_vec - nt
    max_r <- max(abs(total_vec - (th + nt)))
    max_resid <- max(max_resid, max_r)
    if (verbose) cat(sprintf("  %s: max |residual| = %.2e\n", var_names[v], max_r))
  }

  passed <- max_resid < 1e-12
  if (verbose) {
    cat(sprintf("\nOverall max |residual| = %.2e\n", max_resid))
    if (passed) cat("PASSED\n") else cat("WARNING: Additivity violation\n")
  }
  invisible(passed)
}

#' Print TCA Result
#'
#' @param x   A \code{tca_result} object.
#' @param target Target variable to display (default: first intermediate).
#' @param ...   Additional arguments (ignored).
#' @export
print.tca_result <- function(x, target = NULL, ...) {
  if (is.null(target)) target <- 1
  cat(sprintf("\nTCA Results (mode: %s)\n", x$mode))
  cat(sprintf("Shock from: %s | Horizon: %d\n", x$var_names[x$from], x$h))
  cat(sprintf("Response variable: %s\n", x$var_names[target]))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  cat(sprintf("%4s | %12s", "h", "Total"))
  for (nm in x$channel_names) {
    cat(sprintf(" | %12s", substr(nm, 1, 12)))
  }
  cat("\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  horizons <- unique(c(0, 1, 2, 4, 8, 12, 16, 20, x$h))
  horizons <- horizons[horizons <= x$h]

  for (t in horizons) {
    cat(sprintf("%4d | %12.6f", t, x$irf_total[t + 1, target]))
    for (nm in x$channel_names) {
      cat(sprintf(" | %12.6f", x$irf_channels[[nm]][t + 1, target]))
    }
    cat("\n")
  }
  cat(paste(rep("-", 70), collapse = ""), "\n")
  invisible(x)
}