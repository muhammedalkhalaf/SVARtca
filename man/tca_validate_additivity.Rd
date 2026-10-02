\name{tca_validate_additivity}
\alias{tca_validate_additivity}
\title{Validate Binary Additivity}
\usage{
tca_validate_additivity(from, B, Omega, K, h, order,
  var_names = NULL, verbose = TRUE, pair = NULL, tol = 1e-10)
}
\arguments{
\item{from}{Shock variable (1-based).}

\item{B}{Systems form B matrix.}

\item{Omega}{Systems form Omega matrix.}

\item{K}{Number of variables.}

\item{h}{Maximum horizon.}

\item{order}{Transmission ordering.}

\item{var_names}{Character vector of variable names (optional).}

\item{verbose}{Logical; print results? Default \code{TRUE}.}

\item{pair}{Optional integer vector of two distinct variable numbers
(1-based) for which the inclusion-exclusion identity is checked.}

\item{tol}{Tolerance for the maximum absolute residual. Default
\code{1e-10}.}
}
\value{
Invisibly returns a logical, \code{TRUE} if \code{B} is strictly
  lower triangular and every residual is below \code{tol}, with
  attributes \code{residuals} (named numeric vector, the maximum absolute
  additivity residual for each variable), \code{pair_residual} (the
  maximum absolute inclusion-exclusion residual, or \code{NA} if
  \code{pair} is \code{NULL}), \code{max_residual} (the largest of these)
  and \code{lower_triangular} (logical).
}
\description{
Checks Theorem 2(ii) of Wegner, Lieb, Smeekes and Wilms (2025): the
effects of disjoint transmission channels sum to the total effect.
For every variable \code{j}, the effect of the paths passing through
\code{j} (at any horizon) is computed independently of the identity
\code{through = total - not_through}, by partitioning those paths by the
first node of \code{j} they visit and summing the AND-conditioned
effects of the parts (see \code{\link{tca_analyze}} for the AND and NOT
conditions). This independent through effect plus the NOT-conditioned
not-through effect is then compared with the total effect. The residual
is evaluated at all response variables other than \code{j} itself,
because the AND condition sets the effect on the conditioning node to
zero while the effect on \code{j} counts entirely as "through \code{j}"
in the \code{total - not_through} convention.

If \code{pair} is given, the inclusion-exclusion identity used by the
\code{"exhaustive_4way"} mode of \code{\link{tca_analyze}},
\code{through(v1 and v2) = through(v1) + through(v2) - through(v1 or v2)},
is also checked: the left-hand side is computed directly with AND
conditions on both variables (a double first-passage sum of
\code{(h+1)^2} linear solves, so keep \code{h} moderate), the right-hand
side from the independent single-variable through effects and
\code{total - not_through(v1, v2)}. The residual is evaluated at all
response variables other than \code{v1} and \code{v2}.

The function also reports whether \code{B} is strictly lower triangular
(entries on and above the diagonal at most \code{1e-12} times the
largest absolute entry), which the systems form requires for the graph
to be acyclic. Before
version 1.0.3 the residual was computed as
\code{total - ((total - not_through) + not_through)}, which is zero by
construction and could not detect any error.
}
\examples{
Phi0 <- matrix(c(1, 0.3, 0, 0.95), 2, 2)
As <- list(matrix(c(0.5, -0.1, 0.2, 0.4), 2, 2))
sf <- tca_systems_form(Phi0, As, h = 5)
ok <- tca_validate_additivity(1, sf$B, sf$Omega, K = 2, h = 5,
                              order = 1:2, pair = c(1, 2))
attr(ok, "max_residual")
# A backward edge breaks the acyclic structure and is detected
Bc <- sf$B; Bc[1, 4] <- 0.2
tca_validate_additivity(1, Bc, sf$Omega, K = 2, h = 5, order = 1:2)
}
\references{
Wegner, E., Lieb, L., Smeekes, S. and Wilms, I. (2025).
\emph{Transmission Channel Analysis in Dynamic Models.}
arXiv:2405.18987, \doi{10.48550/arXiv.2405.18987}.
}
