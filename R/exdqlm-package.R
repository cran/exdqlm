#' exdqlm: Extended Dynamic Quantile Linear Models
#'
#' Routines for Bayesian estimation and analysis of dynamic quantile linear models
#' using the extended asymmetric Laplace error distribution (exDQLM).
#'
#' @section Runtime options:
#' \itemize{
#'   \item `options(exdqlm.use_cpp_kf = TRUE|FALSE)` – C++ Kalman bridge (optional; default TRUE).
#'   \item `options(exdqlm.use_cpp_samplers = TRUE|FALSE)` – C++ samplers (optional; default FALSE).
#' }
#'
#' @useDynLib exdqlm, .registration = TRUE
#' @import Rcpp
#' @docType package
#' @name exdqlm-package
#' @keywords package
"_PACKAGE"
