#' Density Function for the Extended Asymmetric Laplace (exAL) Distribution
#'
#' Computes the PDF of the Extended Asymmetric Laplace (exAL) distribution.
#' Vectorized over `x`.
#'
#' @param x Numeric vector of quantiles.
#' @param p0 Probability level used in the quantile parametrization. Scalar in (0, 1). Default `0.5`.
#' @param mu Location parameter (scalar). Default `0`.
#' @param sigma Scale parameter (scalar, strictly positive). Default `1`.
#' @param gamma Skewness parameter controlling asymmetry (scalar). Must be within valid bounds implied by `p0`. Default `0`.
#' @param log Logical scalar; if `TRUE` return log-density. Default `FALSE`.
#'
#' @return Numeric vector of densities (same length as `x`).
#'
#' @examples
#' dexal(0)
#' dexal(1, p0 = 0.75, mu = 0, sigma = 2, gamma = 0.25)
#' dexal(seq(-3, 3, by = 0.1), p0 = 0.3, mu = 0, sigma = 1, gamma = -0.45)
#'
#' @export
dexal <- function(x, p0 = 0.5, mu = 0, sigma = 1, gamma = 0, log = FALSE) {
  x <- as.numeric(x)
  vapply(x, function(xx) dexal_cpp(xx, p0, mu, sigma, gamma, log), numeric(1))
}

#' Cumulative Distribution Function (CDF) for the exAL Distribution
#'
#' Vectorized over `q`.
#'
#' @param q Numeric vector of quantiles.
#' @inheritParams dexal
#' @param lower.tail Logical scalar; if `TRUE` (default) return \eqn{P(X \le q)}, otherwise \eqn{P(X > q)}.
#' @param log.p Logical scalar; if `TRUE`, return log-probabilities.
#'
#' @return Numeric vector of CDF values (same length as `q`).
#'
#' @examples
#' pexal(0)
#' pexal(c(-1, 0, 1), p0 = 0.5, mu = 0, sigma = 1, gamma = 0.1)
#'
#' @export
pexal <- function(q, p0 = 0.5, mu = 0, sigma = 1, gamma = 0, lower.tail = TRUE, log.p = FALSE) {
  q <- as.numeric(q)
  vapply(q, function(qq) pexal_cpp(qq, p0, mu, sigma, gamma, lower.tail, log.p), numeric(1))
}

#' Quantile Function for the exAL Distribution
#'
#' Vectorized over `p`.
#'
#' @param p Numeric vector of probabilities in (0, 1).
#' @inheritParams pexal
#'
#' @return Numeric vector of quantiles (same length as `p`).
#'
#' @examples
#' p <- seq(0.1, 0.9, by = 0.2)
#' q <- qexal(p, p0 = 0.5, mu = 0, sigma = 1, gamma = 0)
#' all.equal(p, pexal(q, p0 = 0.5, mu = 0, sigma = 1, gamma = 0), tol = 1e-4)
#'
#' @export
qexal <- function(p, p0 = 0.5, mu = 0, sigma = 1, gamma = 0, lower.tail = TRUE, log.p = FALSE) {
  p <- as.numeric(p)
  vapply(p, function(pp) qexal_cpp(pp, p0, mu, sigma, gamma, lower.tail, log.p), numeric(1))
}

#' Random Sample Generation for the exAL Distribution
#'
#' @param n Positive integer number of samples to draw (scalar).
#' @inheritParams dexal
#'
#' @return Numeric vector of length `n`.
#'
#' @examples
#' set.seed(1); length(rexal(10))
#' rexal(3, p0 = 0.5, mu = c(-1, 0, 1))
#'
#' @export
rexal <- function(n, p0 = 0.5, mu = 0, sigma = 1, gamma = 0) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) stop("n must be a positive integer")
  # check parameter lengths
  if (length(p0) != 1L && length(p0) != n) stop("p0 must have length 1 or n")
  if (length(mu) != 1L && length(mu) != n) stop("mu must have length 1 or n")
  if (length(sigma) != 1L && length(sigma) != n) stop("sigma must have length 1 or n")
  if (length(gamma) != 1L && length(gamma) != n) stop("gamma must have length 1 or n")
  # make all parameters length n
  p0 = rep_len(p0,n)
  mu = rep_len(mu,n)
  sigma = rep_len(sigma,n)
  gamma = rep_len(gamma,n)
  # apply random sampling over new vectors 
  vapply(1:n, function(nn) rexal_cpp(1, p0[nn], mu[nn], sigma[nn], gamma[nn]), numeric(1))
}

