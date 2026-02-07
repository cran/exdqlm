# Reference R implementation (legacy) for comparison in tests
.log_g_ref <- function(gam) {
  base::log(2) + stats::pnorm(-abs(gam), log = TRUE) + 0.5 * gam^2
}
.L_ref <- function(p0) {
  stats::uniroot(function(gam) exp(.log_g_ref(gam)) - (1 - p0), c(-1000, 0))$root
}
.U_ref <- function(p0) {
  stats::uniroot(function(gam) exp(.log_g_ref(gam)) - p0, c(0, 1000))$root
}


test_that("get_gamma_bounds sane at p0=0.5", {
  b <- get_gamma_bounds(0.5)
  expect_length(b, 2)
  expect_true(all(is.finite(b)))
  expect_true(b["L"] < b["U"])
  expect_lt(abs(b["L"]), 10)
  expect_lt(abs(b["U"]), 10)
  expect_lt(abs(b["L"] + b["U"]), 1e-3)

  ref <- c(L = .L_ref(0.5), U = .U_ref(0.5))
  expect_equal(unname(b), unname(ref), tolerance = 1e-3)
})

test_that("get_gamma_bounds works at extreme p0 (regression)", {
  for (p0 in c(0.001, 0.994170922668371)) {
    b <- get_gamma_bounds(p0)
    expect_length(b, 2)
    expect_true(all(is.finite(b)))
    expect_true(b["L"] < b["U"])
    expect_true(b["L"] <= 0)
    expect_true(b["U"] >= 0)

    ref <- c(L = .L_ref(p0), U = .U_ref(p0))
    expect_equal(unname(b), unname(ref), tolerance = 1e-3)

    # C++ wrapper should agree with the reference as well.
    cpp <- exdqlm:::get_gamma_bounds_cpp(p0)
    expect_true(all(is.finite(cpp)))
    expect_true(cpp["L"] < cpp["U"])
    expect_equal(unname(cpp), unname(ref), tolerance = 1e-3)
  }
})


test_that("get_gamma_bounds finite and ordered on grid", {
  pgrid <- c(seq(0.01, 0.99, by = 0.01), 0.001, 1 - 0.001)
  for (p0 in pgrid) {
    b <- get_gamma_bounds(p0)
    expect_true(all(is.finite(b)))
    expect_true(b["L"] < b["U"])
    expect_true(b["L"] <= 0)
    expect_true(b["U"] >= 0)

    # Validate against defining equations on log-scale (stable for extremes).
    expect_lt(abs(.log_g_ref(b["L"]) - log1p(-p0)), 1e-3)
    expect_lt(abs(.log_g_ref(b["U"]) - log(p0)), 1e-3)
  }
})


test_that("get_gamma_bounds stable under repeated calls", {
  set.seed(123)
  p0s <- stats::runif(50, min = 0.001, max = 1 - 0.001)
  for (p0 in p0s) {
    b <- get_gamma_bounds(p0)
    expect_true(all(is.finite(b)))
    expect_true(b["L"] < b["U"])
    expect_lt(abs(.log_g_ref(b["L"]) - log1p(-p0)), 1e-3)
    expect_lt(abs(.log_g_ref(b["U"]) - log(p0)), 1e-3)
  }
})


test_that("ISVB/MCMC reach bounds path on tiny model", {
  # Minimal exdqlm model
  model <- as.exdqlm(list(m0 = 0, C0 = matrix(1, 1, 1), FF = 1, GG = 1))
  y <- c(0.1, -0.2)

  old_opts <- options(exdqlm.use_cpp_kf = FALSE, exdqlm.compute_elbo = FALSE)
  on.exit(options(old_opts), add = TRUE)

  expect_silent(
    exdqlmISVB(y, p0 = 0.5, model = model, df = 1, dim.df = 1,
              sig.init = 1, n.IS = 2, n.samp = 2, tol = 1e6, verbose = FALSE)
  )

  expect_silent(
    exdqlmMCMC(y, p0 = 0.5, model = model, df = 1, dim.df = 1,
              Sig.mh = diag(c(0.005, 0.005)), n.burn = 0, n.mcmc = 1,
              init.from.isvb = FALSE, verbose = FALSE)
  )
})
