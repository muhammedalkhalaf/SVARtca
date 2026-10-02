\name{print.tca_result}
\alias{print.tca_result}
\title{Print TCA Result}
\usage{
\method{print}{tca_result}(x, target = NULL, ...)
}
\arguments{
\item{x}{A \code{tca_result} object.}

\item{target}{Integer index (original ordering) of the response variable
to display. Default \code{NULL}, which displays variable 1 (the first
variable), not the first intermediate.}

\item{...}{Additional arguments (ignored).}
}
\value{
Invisibly returns \code{x}.
}
\description{
S3 print method for \code{tca_result} objects. Displays a formatted
table of transmission channel contributions across horizons.
}
