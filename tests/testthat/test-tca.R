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

# ---- Hand-computed IRFs: C_j = sum_i A_i C_{j-i}, IRF_j = C_j Phi0 ----
hand_irf <- function(As, Phi0, h, shock) {
  K <- nrow(Phi0); p <- length(As)
  C <- list(diag(K))
  for (j in seq_len(h)) {
    Cj <- matrix(0, K, K)
    for (i in seq_len(min(j, p))) Cj <- Cj + As[[i]] %*% C[[j - i + 1]]
    C[[j + 1]] <- Cj
  }
  t(sapply(C, function(Cj) (Cj %*% Phi0)[, shock]))
}

# ---- Binary additivity (Theorem 2(ii)) ----
test_that("binary additivity holds for all variables", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  passed <- tca_validate_additivity(
    from = 1, B = sf$B, Omega = sf$Omega,
    K = K, h = h, order = order_vec,
    var_names = var_names, verbose = FALSE
  )
  expect_true(passed)
  expect_true(attr(passed, "lower_triangular"))
  res <- attr(passed, "residuals")
  expect_length(res, K)
  expect_named(res, var_names)
  expect_true(all(is.finite(res)))
  expect_lt(attr(passed, "max_residual"), 1e-12)
  expect_true(is.na(attr(passed, "pair_residual")))
})

test_that("binary additivity holds under a non-trivial ordering", {
  ord <- c(2, 4, 1, 3)
  sf <- tca_systems_form(Phi0, list(A1), h = 6, order = ord)
  passed <- tca_validate_additivity(3, sf$B, sf$Omega, K = K, h = 6,
                                    order = ord, verbose = FALSE)
  expect_true(passed)
  expect_lt(attr(passed, "max_residual"), 1e-12)
})

test_that("inclusion-exclusion identity holds for a pair", {
  sf <- tca_systems_form(Phi0, list(A1), h = 6, order = order_vec)
  passed <- tca_validate_additivity(1, sf$B, sf$Omega, K = K, h = 6,
                                    order = order_vec, pair = c(2, 4),
                                    verbose = FALSE)
  expect_true(passed)
  expect_false(is.na(attr(passed, "pair_residual")))
  expect_lt(attr(passed, "pair_residual"), 1e-12)
  # the AND-based through(v1 and v2) must equal the 4-way "both" channel
  # at response variables other than v1 and v2
  ns <- asNamespace("SVARtca")
  n2 <- ns$not_vars_for(2, K, 6, order_vec)
  n4 <- ns$not_vars_for(4, K, 6, order_vec)
  th_and <- ns$vec_to_irf(ns$through_both_effect(1, sf$B, sf$Omega, n2, n4),
                          K, 6, order_vec)
  res <- tca_analyze(1, sf$B, sf$Omega, intermediates = c(2, 4), K = K,
                     h = 6, order = order_vec, mode = "exhaustive_4way",
                     var_names = var_names)
  expect_equal(th_and[, c(1, 3)], res$irf_channels[["GDP & Wages"]][, c(1, 3)],
               tolerance = 1e-12)
  expect_error(tca_validate_additivity(1, sf$B, sf$Omega, K = K, h = 6,
                                       order = order_vec, pair = c(2, 2),
                                       verbose = FALSE), "distinct")
})

test_that("a corrupted B is reported with a nonzero residual", {
  sf <- tca_systems_form(Phi0, list(A1), h = 6, order = order_vec)
  # backward edge (above the diagonal): the graph is no longer acyclic
  Bc <- sf$B
  Bc[3, 10] <- 0.3
  bad <- tca_validate_additivity(1, Bc, sf$Omega, K = K, h = 6,
                                 order = order_vec, pair = c(2, 4),
                                 verbose = FALSE)
  expect_false(bad)
  expect_false(attr(bad, "lower_triangular"))
  expect_gt(attr(bad, "max_residual"), 1e-3)
  expect_gt(max(attr(bad, "residuals")), 1e-3)
  expect_gt(attr(bad, "pair_residual"), 1e-3)
  # the printed diagnostic reports the violation
  expect_output(tca_validate_additivity(1, Bc, sf$Omega, K = K, h = 6,
                                        order = order_vec),
                "Additivity violation")
  # a nonzero diagonal entry is caught by the structural check
  Bd <- sf$B
  Bd[7, 7] <- 0.4
  bad2 <- tca_validate_additivity(1, Bd, sf$Omega, K = K, h = 6,
                                  order = order_vec, verbose = FALSE)
  expect_false(bad2)
  expect_false(attr(bad2, "lower_triangular"))
})

# ---- Binary decomposition ----
test_that("tca_decompose_binary matches independent computations", {
  sf <- tca_systems_form(Phi0, list(A1), h = h, order = order_vec)
  dec <- tca_decompose_binary(1, sf$B, sf$Omega, var_idx = 2,
                               K = K, h = h, order = order_vec)
  # total equals the hand-computed MA recursion
  expect_equal(dec$total, hand_irf(list(A1), Phi0, h, 1), tolerance = 1e-12)
  # through equals the AND-based first-passage computation at the other
  # variables, and not_through is the complement
  ns <- asNamespace("SVARtca")
  n2 <- ns$not_vars_for(2, K, h, order_vec)
  th_indep <- ns$vec_to_irf(ns$through_any_effect(1, sf$B, sf$Omega, n2),
                            K, h, order_vec)
  expect_equal(dec$through[, -2], th_indep[, -2], tolerance = 1e-12)
  expect_equal(dec$not_through[, -2], (dec$total - th_indep)[, -2],
               tolerance = 1e-12)
  # effect on the variable itself counts entirely as "through"
  expect_equal(dec$not_through[, 2], rep(0, h + 1))
  expect_equal(dec$through[, 2], dec$total[, 2])
})

# ---- Worked example, Section 2.2 of Wegner, Lieb, Smeekes and Wilms ----
test_that("static three-equation example of Section 2.2 reproduces", {
  a1 <- -0.4; a2 <- 0.3; a3 <- 0.5; a4 <- 1.5
  A <- matrix(c(1, 0, -a1,
                -a2, 1, 0,
                -a3, -a4, 1), 3, 3, byrow = TRUE)
  eta <- 1 - a1 * a3 - a1 * a2 * a4
  sf <- tca_systems_form(solve(A), list(matrix(0, 3, 3)), h = 0, order = 1:3)
  r <- tca_analyze(1, sf$B, sf$Omega, intermediates = 2, K = 3, h = 0,
                   order = 1:3)
  TE <- (a2 * a4 + a3) / eta
  IE <- a2 * a4 / ((1 + a1^2) * eta)
  DE <- (a3 + a1 * (1 - eta)) / ((1 + a1^2) * eta)
  expect_equal(TE, 0.6884057971, tolerance = 1e-9)
  expect_equal(IE, 0.2811094453, tolerance = 1e-9)
  expect_equal(DE, 0.4072963518, tolerance = 1e-9)
  expect_equal(r$irf_total[1, 3], TE, tolerance = 1e-12)
  expect_equal(r$irf_channels[["Through Var2"]][1, 3], IE, tolerance = 1e-12)
  expect_equal(r$irf_channels[["Direct"]][1, 3], DE, tolerance = 1e-12)
  expect_equal(IE + DE, TE, tolerance = 1e-12)
  ok <- tca_validate_additivity(1, sf$B, sf$Omega, K = 3, h = 0,
                                order = 1:3, verbose = FALSE)
  expect_true(ok)
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

# ---- Brute-force path enumeration on the systems form DAG ----
# Returns, for every target node, the sum of path weights (product of the
# B entries along the path times the Omega entry at entry) over all paths
# from the shock whose set of visited variables (excluding the target node
# itself) satisfies `keep(vars_visited)`.
enum_channel <- function(B, Omega, from, K, h, order, keep) {
  n <- nrow(B)
  out <- numeric(n)
  var_of <- function(node) order[((node - 1) %% K) + 1]
  walk <- function(node, weight, visited) {
    out[node] <<- out[node] + if (keep(visited)) weight else 0
    nxt <- which(B[, node] != 0)
    for (m in nxt[nxt > node]) walk(m, weight * B[m, node], c(visited, var_of(node)))
  }
  for (m in which(Omega[, from] != 0)) walk(m, Omega[m, from], integer(0))
  vec_to_irf_local <- function(vec) {
    mat <- matrix(0, h + 1, K)
    for (t in 0:h) {
      idx <- (t * K + 1):((t + 1) * K)
      mat[t + 1, order] <- vec[idx]
    }
    mat
  }
  vec_to_irf_local(out)
}

# ---- Exhaustive 3-way ----
test_that("exhaustive_3way channels match brute-force path enumeration", {
  h2 <- 2
  sf <- tca_systems_form(Phi0, list(A1), h = h2, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h2,
                      order = order_vec, mode = "exhaustive_3way",
                      var_names = var_names)
  expect_equal(length(res$channel_names), 3)
  tot   <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec, function(v) TRUE)
  thr2  <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec, function(v) 2 %in% v)
  only4 <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                        function(v) 4 %in% v && !(2 %in% v))
  direct <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                         function(v) !(2 %in% v) && !(4 %in% v))
  # compare at response variables other than the intermediates
  tg <- c(1, 3)
  expect_equal(res$irf_total[, tg], tot[, tg], tolerance = 1e-12)
  expect_equal(res$irf_channels[["Through GDP (incl.)"]][, tg], thr2[, tg],
               tolerance = 1e-12)
  expect_equal(res$irf_channels[["Through Wages only"]][, tg], only4[, tg],
               tolerance = 1e-12)
  expect_equal(res$irf_channels[["Direct"]][, tg], direct[, tg],
               tolerance = 1e-12)
  # the enumerated channels partition the enumerated total (Theorem 2(ii))
  expect_equal((thr2 + only4 + direct)[, tg], tot[, tg], tolerance = 1e-12)
  # total agrees with the hand MA recursion
  expect_equal(res$irf_total, hand_irf(list(A1), Phi0, h2, 1), tolerance = 1e-12)
})

# ---- Exhaustive 4-way ----
test_that("exhaustive_4way channels match brute-force path enumeration", {
  h2 <- 2
  sf <- tca_systems_form(Phi0, list(A1), h = h2, order = order_vec)
  res <- tca_analyze(from = 1, B = sf$B, Omega = sf$Omega,
                      intermediates = c(2, 4), K = K, h = h2,
                      order = order_vec, mode = "exhaustive_4way",
                      var_names = var_names)
  expect_equal(length(res$channel_names), 4)
  only2 <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                        function(v) 2 %in% v && !(4 %in% v))
  only4 <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                        function(v) 4 %in% v && !(2 %in% v))
  both  <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                        function(v) 2 %in% v && 4 %in% v)
  direct <- enum_channel(sf$B, sf$Omega, 1, K, h2, order_vec,
                         function(v) !(2 %in% v) && !(4 %in% v))
  tg <- c(1, 3)
  expect_equal(res$irf_channels[["GDP only"]][, tg], only2[, tg], tolerance = 1e-12)
  expect_equal(res$irf_channels[["Wages only"]][, tg], only4[, tg], tolerance = 1e-12)
  expect_equal(res$irf_channels[["GDP & Wages"]][, tg], both[, tg], tolerance = 1e-12)
  expect_equal(res$irf_channels[["Direct"]][, tg], direct[, tg], tolerance = 1e-12)
  # the AND-based helper agrees with the enumerated "both" channel
  ns <- asNamespace("SVARtca")
  th_and <- ns$vec_to_irf(
    ns$through_both_effect(1, sf$B, sf$Omega,
                           ns$not_vars_for(2, K, h2, order_vec),
                           ns$not_vars_for(4, K, h2, order_vec)),
    K, h2, order_vec)
  expect_equal(th_and[, tg], both[, tg], tolerance = 1e-12)
  # the four channels partition the total
  ch_sum <- Reduce("+", res$irf_channels)
  expect_equal(ch_sum[, tg], res$irf_total[, tg], tolerance = 1e-12)
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
