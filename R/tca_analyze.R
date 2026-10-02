#' @title Main TCA Analysis Functions
#' @description Primary user-facing functions for running Transmission
#'   Channel Analysis with different decomposition modes.
#' @name tca_main
NULL

#' Transmission Channel Analysis
#'
#' Decomposes impulse response functions into transmission channel
#' contributions using the methodology of Wegner, Lieb, Smeekes and Wilms
#' (2025), \doi{10.48550/arXiv.2405.18987}.
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
#' and the effect not passing through it. The not-through effect is
#' obtained by blocking every node of the variable (NOT condition) and the
#' through effect is defined as \code{total - not_through}, so the three
#' matrices satisfy \code{total = through + not_through} by construction.
#' Use \code{\link{tca_validate_additivity}} to check the decomposition
#' against an independent computation of the through effect.
#'
#' @param from   Shock variable (1-based).
#' @param B      Systems form B matrix.
#' @param Omega   Systems form Omega matrix.
#' @param var_idx Variable to decompose through (1-based).
#' @param K       Number of variables.
#' @param h       Maximum horizon.
#' @param order   Transmission ordering.
#' @return A list with matrices \code{total}, \code{through},
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
#' Checks Theorem 2(ii) of Wegner, Lieb, Smeekes and Wilms (2025): the
#' effects of disjoint transmission channels sum to the total effect.
#' For every variable \code{j}, the effect of the paths passing through
#' \code{j} (at any horizon) is computed independently of the identity
#' \code{through = total - not_through}, by partitioning those paths by the
#' first node of \code{j} they visit and summing the AND-conditioned
#' effects of the parts (see \code{\link{tca_analyze}} for the AND and NOT
#' conditions). This independent through effect plus the NOT-conditioned
#' not-through effect is then compared with the total effect. The residual
#' is evaluated at all response variables other than \code{j} itself,
#' because the AND condition sets the effect on the conditioning node to
#' zero while the effect on \code{j} counts entirely as "through \code{j}"
#' in the \code{total - not_through} convention.
#'
#' If \code{pair} is given, the inclusion-exclusion identity used by the
#' \code{"exhaustive_4way"} mode of \code{\link{tca_analyze}},
#' \code{through(v1 and v2) = through(v1) + through(v2) - through(v1 or v2)},
#' is also checked: the left-hand side is computed directly with AND
#' conditions on both variables (a double first-passage sum of
#' \code{(h+1)^2} linear solves, so keep \code{h} moderate), the right-hand
#' side from the independent single-variable through effects and
#' \code{total - not_through(v1, v2)}. The residual is evaluated at all
#' response variables other than \code{v1} and \code{v2}.
#'
#' The function also reports whether \code{B} is strictly lower triangular
#' (entries on and above the diagonal at most \code{1e-12} times the
#' largest absolute entry), which the systems form requires for the graph
#' to be acyclic. Before
#' version 1.0.3 the residual was computed as
#' \code{total - ((total - not_through) + not_through)}, which is zero by
#' construction and could not detect any error.
#'
#' @param from      Shock variable (1-based).
#' @param B         Systems form B matrix.
#' @param Omega     Systems form Omega matrix.
#' @param K         Number of variables.
#' @param h         Maximum horizon.
#' @param order     Transmission ordering.
#' @param var_names Character vector of variable names (optional).
#' @param verbose   Logical; print results? Default \code{TRUE}.
#' @param pair      Optional integer vector of two distinct variable numbers
#'   (1-based) for which the inclusion-exclusion identity is checked.
#' @param tol       Tolerance for the maximum absolute residual. Default
#'   \code{1e-10}.
#' @return Invisibly returns a logical, \code{TRUE} if \code{B} is strictly
#'   lower triangular and every residual is below \code{tol}, with
#'   attributes \code{residuals} (named numeric vector, the maximum absolute
#'   additivity residual for each variable), \code{pair_residual} (the
#'   maximum absolute inclusion-exclusion residual, or \code{NA} if
#'   \code{pair} is \code{NULL}), \code{max_residual} (the largest of these)
#'   and \code{lower_triangular} (logical).
#' @export
#' @examples
#' Phi0 <- matrix(c(1, 0.3, 0, 0.95), 2, 2)
#' As <- list(matrix(c(0.5, -0.1, 0.2, 0.4), 2, 2))
#' sf <- tca_systems_form(Phi0, As, h = 5)
#' ok <- tca_validate_additivity(1, sf$B, sf$Omega, K = 2, h = 5,
#'                               order = 1:2, pair = c(1, 2))
#' attr(ok, "max_residual")
#' # A backward edge breaks the acyclic structure and is detected
#' Bc <- sf$B; Bc[1, 4] <- 0.2
#' tca_validate_additivity(1, Bc, sf$Omega, K = 2, h = 5, order = 1:2)
tca_validate_additivity <- function(from, B, Omega, K, h, order,
                                     var_names = NULL, verbose = TRUE,
                                     pair = NULL, tol = 1e-10) {
  if (is.null(var_names)) var_names <- paste0("Var", seq_len(K))
  total_vec <- transmissionEffect(from, B, Omega)
  n <- length(total_vec)

  # strictly lower triangular up to rounding noise from the LD decomposition
  upper_max <- max(abs(B[upper.tri(B, diag = TRUE)]))
  lower_tri <- upper_max <= 1e-12 * max(1, max(abs(B)))

  if (verbose) cat("===== Binary Additivity Test =====\n")
  residuals <- numeric(K)
  names(residuals) <- var_names
  th_indep <- vector("list", K)
  for (v in seq_len(K)) {
    nv <- not_vars_for(v, K, h, order)
    nt <- transmissionEffect(from, B, Omega, not_vars = nv)
    th <- through_any_effect(from, B, Omega, nv)
    th_indep[[v]] <- th
    keep <- setdiff(seq_len(n), nv)
    residuals[v] <- max(abs((total_vec - (th + nt))[keep]))
    if (verbose) cat(sprintf("  %s: max |total - (through + not_through)| = %.2e\n",
                             var_names[v], residuals[v]))
  }

  pair_resid <- NA_real_
  if (!is.null(pair)) {
    if (length(pair) != 2 || pair[1] == pair[2] || any(!pair %in% seq_len(K))) {
      stop("'pair' must contain two distinct variable numbers in 1:K.")
    }
    n1 <- not_vars_for(pair[1], K, h, order)
    n2 <- not_vars_for(pair[2], K, h, order)
    th_and <- through_both_effect(from, B, Omega, n1, n2)
    nt_both <- transmissionEffect(from, B, Omega, not_vars = c(n1, n2))
    th_or <- total_vec - nt_both
    rhs <- th_indep[[pair[1]]] + th_indep[[pair[2]]] - th_or
    keep <- setdiff(seq_len(n), c(n1, n2))
    pair_resid <- max(abs((th_and - rhs)[keep]))
    if (verbose) cat(sprintf("  %s and %s: max |inclusion-exclusion residual| = %.2e\n",
                             var_names[pair[1]], var_names[pair[2]], pair_resid))
  }

  max_resid <- max(c(residuals, pair_resid), na.rm = TRUE)
  passed <- lower_tri && max_resid < tol
  if (verbose) {
    cat(sprintf("\nB strictly lower triangular: %s\n", if (lower_tri) "yes" else "NO"))
    cat(sprintf("Overall max |residual| = %.2e (tolerance %.1e)\n", max_resid, tol))
    if (passed) cat("PASSED\n") else cat("WARNING: Additivity violation\n")
  }
  attr(passed, "residuals") <- residuals
  attr(passed, "pair_residual") <- pair_resid
  attr(passed, "max_residual") <- max_resid
  attr(passed, "lower_triangular") <- lower_tri
  invisible(passed)
}

#' Print TCA Result
#'
#' S3 print method for \code{tca_result} objects. Displays a formatted
#' table of transmission channel contributions across horizons.
#'
#' @param x   A \code{tca_result} object.
#' @param target Integer index (original ordering) of the response variable
#'   to display. Default \code{NULL}, which displays variable 1 (the first
#'   variable), not the first intermediate.
#' @param ...   Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
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