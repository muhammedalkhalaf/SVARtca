\name{SVARtca-package}
\alias{SVARtca}
\alias{SVARtca-package}
\docType{package}
\title{Transmission Channel Analysis in Structural VAR Models}
\description{
The SVARtca package implements Transmission Channel Analysis (TCA) for
decomposing impulse response functions in structural vector autoregressive
(SVAR) models. The package follows the methodology presented in Wegner,
Lieb, Smeekes (2025).

TCA enables researchers to decompose the total impulse response to a
shock into contributions from different transmission channels. Supported
channels include effects passing through specific intermediate variables,
allowing analysts to distinguish direct effects from indirect effects
that operate through identified transmission mechanisms.

Key features:
\itemize{
  \item Multiple decomposition modes: overlapping, exhaustive 3-way, and
        exhaustive 4-way.
  \item Support for VAR and SVAR models via integration with the
        \pkg{vars} package.
  \item Visualization of channel contributions using \pkg{ggplot2}.
  \item Diagnostic tools for validating decomposition accuracy.
}
}
\section{Main Functions}{

\describe{
  \item{\code{\link{tca_systems_form}}}{Build the systems form representation
    (B and Omega matrices) from VAR coefficients and structural impact matrix.}

  \item{\code{\link{tca_analyze}}}{Run the main TCA analysis with specified
    decomposition mode.}

  \item{\code{\link{tca_decompose_binary}}}{Decompose IRF into "through" and
    "not_through" components for a single variable.}

  \item{\code{\link{tca_validate_additivity}}}{Validate that the binary
    decomposition satisfies exact additivity.}

  \item{\code{\link{tca_from_var}}}{Convenience wrapper to run TCA directly
    from a fitted VAR model object.}

  \item{\code{\link{plot_tca}}}{Create publication-quality visualizations
    of channel contributions.}
}
}
\section{Identification}{

The package supports two identification schemes:
\itemize{
  \item \strong{Cholesky}: Lower-triangular Cholesky decomposition of the
        residual covariance matrix. Set \code{Phi0 = t(chol(Sigma))}.

  \item \strong{Manual}: User-supplied structural impact matrix for
        alternative identifications (e.g., sign restrictions, zero
        restrictions).
}
}
\section{Citation}{

If you use the SVARtca package in your research, please cite the original
methodology paper:

Wegner, E., Lieb, L., Smeekes, S. (2025).
\emph{Transmission Channel Analysis in Dynamic Models.}
arXiv:2405.18987.
\url{https://github.com/enweg/tca-matlab-toolbox}
}
\author{
Muhammad Alkhalaf
}
\references{
Wegner, E., Lieb, L., Smeekes, S. (2025).
Transmission Channel Analysis in Dynamic Models.
arXiv:2405.18987.
}
\keyword{package}
