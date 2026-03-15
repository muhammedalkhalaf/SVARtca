\name{tca_validate_additivity}
\alias{tca_validate_additivity}
\title{Validate Binary Additivity}
\usage{
tca_validate_additivity(from, B, Omega, K, h, order,
  var_names = NULL, verbose = TRUE)
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
}
\value{
Invisibly returns \code{TRUE} if all tests pass.
}
\description{
Tests that total = through(j) + not_through(j) holds for all
variables at machine precision. This is a diagnostic to verify
correct implementation.
}
