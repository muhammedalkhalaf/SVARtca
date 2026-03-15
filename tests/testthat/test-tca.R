# ============================================================
# Unit tests for the TCA package
# Run with: testthat::test_dir("tests/testthat")
# ============================================================

# ---- Setup ----
K <- 4
A1 <- matrix(c(0.7, -0.1, 0.05, -0.05,
               -0.3, 0.6, 0.10, -0.10,
               -0.2, 0.1, 0.70, 0.05,
               -0.1, 0.2, 0.05, 0.65), K, K, byrow = TRUE)

Sigma <- matrix(c(1.00, 0.30, 0.20, 0.10,
                  0.30, 1.50, 0.25, 0.15,
                  0.20, 0.25, 0.80, 0.10,
                  0.10, 0.15, 0.10, 0.60), K, K, byrow = TRUE)

Phi0 <- t(chol(Sigma))
h <- 20
order_vec <- 1:K
var_names <- c("IntRate", "GDP", "Inflation", "Wages")

# ---- Systems form ----
test_that("tca_systems_form produces correct dimensions", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  n <- K * (h + 1)
  expect_equal(nrow(sf$B), n)
  expect_equal(ncol(sf$B), n)
  expect_equal(nrow(sf$Omega), n)
  expect_equal(ncol(sf$Omega), n)
})

test_that("tca_systems_form with two lags", {
  A2 <- matrix(c(0.1, 0, 0, 0,
                 -0.1, 0.1, 0, 0,
                 -0.1, 0, 0.1, 0,
                  0, 0, 0, 0.1), K, K, byrow = TRUE)
  sf <- tca_systems_form(Phi0, list(A1, A2), h = h, order = order_vec)
  n <- K * (h + 1)
  expect_equal(nrow(sf$B), n)
  expect_equal(ncol(sf$B), n)
})

# ---- Binary additivity ----
test_that("binary additivity holds for all variables", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  passed <- tca_validate_additivity(
    from = 1, B = sf$B, Omega = sf$Omega,
    K = K, h = h, order = order_vec,
    var_names = var_names, verbose = FALSE
  )
  expect_true(passed)
})

# ---- Binary decomposition ----
test_that("tca_decompose_binary sums to total", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  dec <- tca_decompose_binary(1, sf$B, sf$Omega, var_idx = 2,
                               K = K, h = h, order = order_vec)
  residual <- max(abs(dec$total - (dec$through + dec$not_through)))
  expect_lt(residual, 1e-12)
})

# ---- Overlapping mode ----
test_that("overlapping mode produces correct number of channels", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h,
                      order = order_vec, mode = "overlapping",
                      var_names = var_names)
  expect_s3_class(res, "tca_result")
  expect_equal(length(res$channel_names), 3)  # ThruGDP, ThruWages, Direct
  expect_equal(res$mode, "overlapping")
})

# ---- Exhaustive 3-way ----
test_that("exhaustive_3way sums to total", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h,
                      order = order_vec, mode = "exhaustive_3way",
                      var_names = var_names)
  expect_equal(length(res$channel_names), 3)

  # Check additivity: sum of channels = total for target = 3 (Inflation)
  target <- 3
  ch_sum <- Reduce("+", lapply(res$irf_channels, function(m) m[, target]))
  max_resid <- max(abs(ch_sum - res$irf_total[, target]))
  expect_lt(max_resid, 1e-10)
})

# ---- Exhaustive 4-way ----
test_that("exhaustive_4way sums to total", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h,
                      order = order_vec, mode = "exhaustive_4way",
                      var_names = var_names)
  expect_equal(length(res$channel_names), 4)

  target <- 3
  ch_sum <- Reduce("+", lapply(res$irf_channels, function(m) m[, target]))
  max_resid <- max(abs(ch_sum - res$irf_total[, target]))
  expect_lt(max_resid, 1e-10)
})

# ---- Reference values from verified Python implementation ----
test_that("results match Python reference values", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h,
                      order = order_vec, mode = "overlapping",
                      var_names = var_names)

  tol <- 1e-6

  # h=0 total Inflation (row 1)
  expect_equal(res$irf_total[1, 3], 0.2000000000, tolerance = tol)

  # h=20 total Inflation (row 21)
  expect_equal(res$irf_total[21, 3], -0.0443717030, tolerance = tol)

  # h=8 total Inflation (row 9)
  expect_equal(res$irf_total[9, 3], -0.2539489911, tolerance = tol)

  # h=8 Through GDP -> Inflation
  ch_gdp <- res$irf_channels[["Through GDP"]]
  expect_equal(ch_gdp[9, 3], -0.1542802189, tolerance = tol)

  # h=20 Through GDP -> Inflation
  expect_equal(ch_gdp[21, 3], -0.0413144427, tolerance = tol)

  # h=0 Direct -> Inflation
  ch_direct <- res$irf_channels[["Direct"]]
  expect_equal(ch_direct[1, 3], 0.1595744681, tolerance = tol)

  # h=20 Direct -> Inflation
  expect_equal(ch_direct[21, 3], -0.0008756125, tolerance = tol)
})

# ---- Input validation ----
test_that("exhaustive modes require exactly 2 intermediates", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  expect_error(
    tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                intermediates = c(2, 3, 4), K = K, h = h,
                order = order_vec, mode = "exhaustive_3way"),
    "exactly 2 intermediate"
  )
})

test_that("invalid mode is rejected", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  expect_error(
    tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                intermediates = c(2, 4), K = K, h = h,
                order = order_vec, mode = "invalid_mode"),
    "arg"
  )
})

# ---- Plot function ----
test_that("plot_tca returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h,
                      order = order_vec, mode = "exhaustive_4way",
                      var_names = var_names)
  p <- plot_tca(res, target = 3)
  expect_s3_class(p, "gg")
})
