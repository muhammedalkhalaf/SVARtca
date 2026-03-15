\name{tca_decompose_binary}
\alias{tca_decompose_binary}
\title{Binary Decomposition: Total = Through + Not-Through}
\usage{
tca_decompose_binary(from, B, Omega, var_idx, K, h, order)
}
\arguments{
\item{from}{Shock variable (1-based).}

\item{B}{Systems form B matrix.}

\item{Omega}{Systems form Omega matrix.}

\item{var_idx}{Variable to decompose through (1-based).}

\item{K}{Number of variables.}

\item{h}{Maximum horizon.}

\item{order}{Transmission ordering.}
}
\value{
A list with matrices \code{total}, \code{through},
  \code{not_through} (each (h+1) x K).
}
\description{
Decomposes the total IRF into the effect passing through a variable
and the effect not passing through it. This is an exact decomposition
(residual is zero at machine precision).
}
