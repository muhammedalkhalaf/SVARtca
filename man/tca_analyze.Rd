\name{tca_analyze}
\alias{tca_analyze}
\title{Transmission Channel Analysis}
\usage{
tca_analyze(from, B, Omega, intermediates, K, h, order,
  mode = "overlapping", var_names = NULL)
}
\arguments{
\item{from}{Shock variable number (1-based).}

\item{B}{Systems form B matrix (from \code{\link{tca_systems_form}}).}

\item{Omega}{Systems form Omega matrix.}

\item{intermediates}{Integer vector of intermediate variable numbers
(1-based, original ordering).}

\item{K}{Number of variables.}

\item{h}{Maximum horizon.}

\item{order}{Transmission ordering vector.}

\item{mode}{Decomposition mode: \code{"overlapping"},
\code{"exhaustive_3way"}, or \code{"exhaustive_4way"}.}

\item{var_names}{Character vector of variable names (optional).}
}
\value{
A list of class \code{"tca_result"} with components:
  \describe{
    \item{irf_total}{Matrix (h+1) x K of total IRFs.}
    \item{irf_channels}{Named list of channel IRF matrices, each (h+1) x K.}
    \item{channel_names}{Character vector of channel names.}
    \item{mode}{Decomposition mode used.}
    \item{from}{Shock variable number.}
    \item{K}{Number of variables.}
    \item{h}{Maximum horizon.}
    \item{order}{Transmission ordering.}
    \item{var_names}{Variable names.}
  }
}
\description{
Decomposes impulse response functions into transmission channel
contributions using the methodology of Wegner, Lieb, Smeekes (2025).

Three decomposition modes are supported:
\describe{
  \item{\code{"overlapping"}}{Each channel is through(j) = total -
    not_through(j). Channels may overlap, so their sum may differ
    from the total.}
  \item{\code{"exhaustive_3way"}}{(2 intermediates only) Non-overlapping:
    (1) through var1 inclusive, (2) through var2 only,
    (3) direct. Sum equals total.}
  \item{\code{"exhaustive_4way"}}{(2 intermediates only) Full
    inclusion-exclusion: (1) var1 only, (2) var2 only,
    (3) both, (4) direct. Sum equals total.}
}
}
\examples{
# Monetary policy model
K <- 4
A1 <- matrix(c(0.7,-0.1,0.05,-0.05, -0.3,0.6,0.10,-0.10,
                -0.2,0.1,0.70,0.05, -0.1,0.2,0.05,0.65), K, K, byrow=TRUE)
Sigma <- matrix(c(1,0.3,0.2,0.1, 0.3,1.5,0.25,0.15,
                   0.2,0.25,0.8,0.1, 0.1,0.15,0.1,0.6), K, K, byrow=TRUE)
Phi0 <- t(chol(Sigma))
sf <- tca_systems_form(Phi0, list(A1), h = 20)
result <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                       intermediates = c(2, 4), K = K, h = 20,
                       order = 1:K, mode = "exhaustive_4way",
                       var_names = c("IntRate","GDP","Inflation","Wages"))
print(result)
}
\references{
Wegner, E., Lieb, L., Smeekes, S. (2025).
\emph{Transmission Channel Analysis in Dynamic Models.}
arXiv:2405.18987.
}
