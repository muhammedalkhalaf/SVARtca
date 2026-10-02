## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)


## ----quickstart---------------------------------------------------------------
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

# Step 1: Build systems form
sf <- tca_systems_form(Phi0, list(A1), h = 20)

# Step 2: Run TCA
result <- tca_analyze(
  from          = 1,
  B             = sf$B,
  Omega         = sf$Omega,
  intermediates = c(2, 4),
  K             = K,
  h             = 20,
  order         = 1:K,
  mode          = "exhaustive_4way",
  var_names     = var_names
)

print(result, target = 3)


## ----overlapping--------------------------------------------------------------
res_ov <- tca_analyze(
  from = 1, B = sf$B, Omega = sf$Omega,
  intermediates = c(2, 4), K = K, h = 20,
  order = 1:K, mode = "overlapping", var_names = var_names
)
print(res_ov, target = 3)


## ----three_way----------------------------------------------------------------
res_3w <- tca_analyze(
  from = 1, B = sf$B, Omega = sf$Omega,
  intermediates = c(2, 4), K = K, h = 20,
  order = 1:K, mode = "exhaustive_3way", var_names = var_names
)
print(res_3w, target = 3)


## ----four_way-----------------------------------------------------------------
res_4w <- tca_analyze(
  from = 1, B = sf$B, Omega = sf$Omega,
  intermediates = c(2, 4), K = K, h = 20,
  order = 1:K, mode = "exhaustive_4way", var_names = var_names
)
print(res_4w, target = 3)


## ----plot_bar-----------------------------------------------------------------
plot_tca(res_4w, target = 3, type = "bar")


## ----plot_line----------------------------------------------------------------
plot_tca(res_ov, target = 3, type = "line")


## ----validate-----------------------------------------------------------------
tca_validate_additivity(
  from = 1, B = sf$B, Omega = sf$Omega,
  K = K, h = 20, order = 1:K, var_names = var_names
)
sf8 <- tca_systems_form(Phi0, list(A1), h = 8, order = 1:K)
ok <- tca_validate_additivity(
  from = 1, B = sf8$B, Omega = sf8$Omega,
  K = K, h = 8, order = 1:K, var_names = var_names, pair = c(2, 4),
  verbose = FALSE
)
attr(ok, "pair_residual")


## ----vars_example, eval = FALSE-----------------------------------------------
# library(vars)
# data(Canada)
# var_est <- VAR(Canada, p = 2, type = "const")
# result <- tca_from_var(var_est, from = "e",
#                         intermediates = c("prod", "rw"),
#                         h = 20, mode = "exhaustive_4way")
# plot_tca(result, target = "U")

