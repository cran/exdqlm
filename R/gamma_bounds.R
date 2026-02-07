#' Bounds for the exAL shape parameter gamma
#'
#' Returns valid lower/upper bounds \code{(L, U)} for the shape parameter \code{gamma}
#' of the standardized extended Asymmetric Laplace (exAL), given \code{p0} in (0,1).
#'
#' This is a user-facing convenience wrapper around the C++ routine
#' \code{get_gamma_bounds_cpp()}, which performs the actual computation.
#'
#' @param p0 Numeric scalar in (0, 1); typically the target quantile level.
#' @return A numeric vector of length 2 named \code{c("L","U")}.
#' @examples
#' get_gamma_bounds(0.5)
#' get_gamma_bounds(0.9)
#' @export
get_gamma_bounds <- function(p0) {
  .gamma_bounds(p0)
}
