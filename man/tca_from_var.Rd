\name{tca_from_var}
\alias{tca_from_var}
\title{Run TCA from a VAR Estimation Object}
\usage{
tca_from_var(var_model, from, intermediates, h = 20,
  order = NULL, mode = "overlapping",
  identification = "cholesky", Phi0 = NULL)
}
\arguments{
\item{var_model}{A fitted VAR model object from \code{\link[vars]{VAR}}.}

\item{from}{Shock variable (integer or name).}

\item{intermediates}{Integer vector or character vector of intermediate
variable names.}

\item{h}{Maximum horizon (default: 20).}

\item{order}{Transmission ordering (default: variable ordering in the VAR).}

\item{mode}{Decomposition mode: \code{"overlapping"},
\code{"exhaustive_3way"}, or \code{"exhaustive_4way"}.}

\item{identification}{Identification scheme: \code{"cholesky"} (default)
or \code{"manual"}.}

\item{Phi0}{Manual impact matrix. Required if \code{identification = "manual"}.}
}
\value{
A \code{tca_result} object (see \code{\link{tca_analyze}}).
}
\description{
Extracts coefficient matrices and residual covariance from a fitted
VAR model (from the \pkg{vars} package) and runs TCA.
}
\examples{
\dontrun{
library(vars)
data(Canada)
var_est <- VAR(Canada, p = 2, type = "const")
result <- tca_from_var(var_est, from = "e",
                        intermediates = c("prod", "rw"),
                        h = 20, mode = "exhaustive_4way")
plot_tca(result, target = "U")
}
}
