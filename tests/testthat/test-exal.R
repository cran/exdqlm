test_that("exAL wrappers return numeric/finites and right shapes", {
  expect_true(is.numeric(dexal(0)))
  expect_true(is.numeric(pexal(0)))
  expect_true(is.numeric(qexal(0.5)))
  x <- rexal(5)
  expect_type(x, "double")
  expect_length(x, 5)
  expect_false(any(!is.finite(x)))
})

test_that("qexal and pexal are approximately inverse in AL case", {
  p  <- seq(0.1, 0.9, by = 0.2)
  q  <- qexal(p, p0 = 0.5, mu = 0, sigma = 1, gamma = 0)
  p2 <- pexal(q, p0 = 0.5, mu = 0, sigma = 1, gamma = 0)
  expect_equal(p2, p, tolerance = 1e-4)
})

test_that("dexal works with a valid gamma bound", {
  expect_silent(dexal(1, p0 = 0.75, sigma = 2, gamma = 0.25))
})
