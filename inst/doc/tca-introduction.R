### code from vignettes/tca-introduction.Rmd ###
library(SVARtca)

# Define a 4-variable VAR(1) model
K <- 4
A1 <- matrix(c( 0.7, -0.1,  0.05, -0.05,
                -0.3,  0.6,  0.10, -0.10,
                -0.2,  0.1,  0.70,  0.05,
                -0.1,  0.2,  0.05,  0.65), K, K, byrow = TRUE)

Sigma <- matrix(c(1.00, 0.30, 0.20, 0.10,
                  0.30, 1.50, 0.25, 0.15,
                  0.20, 0.25, 0.80, 0.10,
                  0.10, 0.15, 0.10, 0.60), K, K, byrow = TRUE)

Phi0 <- t(chol(Sigma))
var_names <- c("IntRate", "GDP", "Inflation", "Wages")

sf <- tca_systems_form(Phi0, list(A1), h = 20)

result <- tca_analyze(
  from          = 1,
  B             = sf\,
  Omega         = sf\,
  intermediates = c(2, 4),
  K             = K,
  h             = 20,
  order         = 1:K,
  mode          = "exhaustive_4way",
  var_names     = var_names
)