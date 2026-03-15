\name{tca_systems_form}
\alias{tca_systems_form}
\title{Build Complete Systems Form}
\usage{
tca_systems_form(Phi0, As, h, order = NULL, Psis = NULL)
}
\arguments{
\item{Phi0}{Structural impact matrix (K x K). For Cholesky
identification, use \code{t(chol(Sigma))}.}

\item{As}{List of VAR coefficient matrices (A_1, A_2, ..., A_p).}

\item{h}{Maximum IRF horizon.}

\item{order}{Transmission ordering vector (default: \code{1:K}).}

\item{Psis}{List of reduced-form MA coefficients (optional, for DSGE/VARMA models).}
}
\value{
A list with components:
  \describe{
    \item{B}{Systems form B matrix (K*(h+1) x K*(h+1)).}
    \item{Omega}{Systems form Omega matrix (K*(h+1) x K*(h+1)).}
  }
}
\description{
Constructs both B and Omega matrices for the systems form
representation x = Bx + Omega*epsilon, following the methodology of
Wegner, Lieb, Smeekes (2025).
}
\examples{
# 2-variable VAR(1)
Phi0 <- matrix(c(1, 0.3, 0, 0.95), 2, 2)
As <- list(matrix(c(0.5, -0.1, 0.2, 0.4), 2, 2))
sf <- tca_systems_form(Phi0, As, h = 10)
dim(sf$B)  # 22 x 22
}
\references{
Wegner, E., Lieb, L., Smeekes, S. (2025).
\emph{Transmission Channel Analysis in Dynamic Models.}
arXiv:2405.18987.
\url{https://github.com/enweg/tca-matlab-toolbox}
}
