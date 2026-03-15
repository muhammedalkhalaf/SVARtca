\name{plot_tca}
\alias{plot_tca}
\title{Plot TCA Channel Decomposition}
\usage{
plot_tca(x, target = NULL, type = NULL, title = NULL,
  colors = NULL)
}
\arguments{
\item{x}{A \code{tca_result} object from \code{\link{tca_analyze}}.}

\item{target}{Integer index of the response variable to plot
(default: first variable after shock).}

\item{type}{Plot type: \code{"bar"} for stacked bar chart (default
for exhaustive modes), \code{"line"} for overlaid lines (default
for overlapping mode). If \code{NULL}, chosen automatically.}

\item{title}{Custom plot title (optional).}

\item{colors}{Named character vector of colours for channels (optional).}
}
\value{
A \code{ggplot} object.
}
\description{
Creates a stacked bar chart (for exhaustive modes) or a line chart
(for overlapping mode) showing channel contributions over horizons.
}
\examples{
K <- 4
A1 <- matrix(c(0.7,-0.1,0.05,-0.05, -0.3,0.6,0.10,-0.10,
                -0.2,0.1,0.70,0.05, -0.1,0.2,0.05,0.65), K, K, byrow=TRUE)
Sigma <- matrix(c(1,0.3,0.2,0.1, 0.3,1.5,0.25,0.15,
                   0.2,0.25,0.8,0.1, 0.1,0.15,0.1,0.6), K, K, byrow=TRUE)
Phi0 <- t(chol(Sigma))
sf <- tca_systems_form(Phi0, list(A1), h = 20)
res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                    intermediates = c(2, 4), K = K, h = 20,
                    order = 1:K, mode = "exhaustive_4way",
                    var_names = c("IntRate","GDP","Inflation","Wages"))
plot_tca(res, target = 3)
}
