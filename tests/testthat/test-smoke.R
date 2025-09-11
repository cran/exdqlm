test_that("exdqlm loads and key symbols exist", {
  # Package can be loaded quietly
  expect_true("package:exdqlm" %in% search() || requireNamespace("exdqlm", quietly = TRUE))

  # A couple of exported symbols exist (adjust if names differ)
  ns <- asNamespace("exdqlm")
  expect_true(exists("exdqlmISVB",  where = ns, inherits = FALSE))
  expect_true(exists("exdqlmMCMC",  where = ns, inherits = FALSE))
  expect_true(exists("exdqlmForecast", where = ns, inherits = FALSE))
})
