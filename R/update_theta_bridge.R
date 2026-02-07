# R/update_theta_bridge.R
# Bridge for univariate (J = 0) only. Keeps the R-visible contract identical to update_theta().
#' Internal bridge for univariate KF used by exdqlmISVB.
#' @keywords internal
#' @noRd
update_theta_bridge <- function(ex.f, ex.q, GG, FF, y, m0, C0, df.mat,
                                use_kstep = FALSE, k = 0L,
                                debug_shapes = FALSE) {
  # Infer base sizes
  stopifnot(length(dim(GG)) == 3L)
  p  <- dim(GG)[1]
  TT <- dim(GG)[3]
  stopifnot(dim(GG)[2] == p)

  # Coerce R shapes -> C++ shapes
  GG_arr <- array(GG, dim = c(p, p, TT))

  # FF: p x TT  ->  p x 1 x TT
  if (is.matrix(FF) && nrow(FF) == p && ncol(FF) == TT) {
    FF_arr <- array(FF, dim = c(p, 1L, TT))
    FF_mat <- FF
  } else if (length(dim(FF)) == 3L && all(dim(FF) == c(p, 1L, TT))) {
    FF_arr <- FF
    FF_mat <- matrix(FF[, 1L, ], nrow = p, ncol = TT)
  } else {
    stop("FF has unexpected shape; expected p x TT or p x 1 x TT.")
  }

  # ex.f: length TT -> 1 x TT
  if (is.null(dim(ex.f))) {
    stopifnot(length(ex.f) == TT)
    ex_f <- matrix(as.numeric(ex.f), nrow = 1L)
  } else if (is.matrix(ex.f) && all(dim(ex.f) == c(1L, TT))) {
    ex_f <- ex.f
  } else stop("ex.f has unexpected shape; expected length TT or 1 x TT.")

  # ex.q: length TT -> 1 x 1 x TT
  if (is.null(dim(ex.q))) {
    stopifnot(length(ex.q) == TT)
    ex_q <- array(as.numeric(ex.q), dim = c(1L, 1L, TT))
  } else if (length(dim(ex.q)) == 3L && all(dim(ex.q) == c(1L, 1L, TT))) {
    ex_q <- ex.q
  } else stop("ex.q has unexpected shape; expected length TT or 1 x 1 x TT.")

  # y: vector length TT OR matrix TT x 1 -> 1 x TT
  if (is.null(dim(y))) {
    stopifnot(length(y) == TT)
    y_mat <- matrix(as.numeric(y), nrow = 1L)
  } else if (is.matrix(y) && all(dim(y) == c(TT, 1L))) {
    y_mat <- t(y)
  } else if (is.matrix(y) && all(dim(y) == c(1L, TT))) {
    y_mat <- y
  } else stop("y has unexpected shape; expected length TT, TT x 1, or 1 x TT.")

  # Prior bits
  stopifnot(is.numeric(m0), length(m0) == p)
  stopifnot(is.matrix(C0), all(dim(C0) == c(p, p)))
  stopifnot(is.matrix(df.mat), all(dim(df.mat) == c(p, p)))

  # helper matrices (k-step/DM unused in univariate path but kept for API parity)
  Ones      <- matrix(1, p, p)
  df.mat.k  <- df.mat
  ppx       <- 0L
  dM        <- 1L
  if (!use_kstep) k <- 0L

  if (debug_shapes) {
    msg <- paste0(
      "BRIDGE shapes - GG:", paste(dim(GG_arr), collapse = "x"),
      " FF:", paste(dim(FF_arr), collapse = "x"),
      " ex_f:", paste(dim(ex_f), collapse = "x"),
      " ex_q:", paste(dim(ex_q), collapse = "x"),
      " y:",   paste(dim(y_mat), collapse = "x")
    )
    message(msg); utils::flush.console()
  }

  # Call C++
  res <- update_theta_cpp(
    GG = GG_arr,
    m0 = as.numeric(m0),
    C0 = C0,
    ex_f = ex_f,
    ex_q = ex_q,
    FF   = FF_arr,
    y    = y_mat,
    ex_df_mat   = df.mat,
    ex_df_mat_k = df.mat.k,
    Ones = Ones,
    p = p, J = 0L, ppx = ppx, TT = TT,
    k = k, dM = dM
  )

  # Map back to R contract identical to update_theta()
  sm  <- res$sm               # p x TT
  sC  <- res$sC               # p x p x TT
  fm  <- res$fm               # p x TT
  fC  <- res$fC               # p x p x TT
  sfe <- as.numeric(res$standard_forecast_errors)  # 1 x TT -> length TT

  # exps, vars, exps2 in univariate case
  exps <- colSums(FF_mat * sm)
  vars <- vapply(seq_len(TT), function(t) {
    as.numeric(crossprod(FF_mat[, t], sC[, , t] %*% FF_mat[, t]))
  }, numeric(1))
  exps2 <- exps^2 + vars

  # \\theta-entropy: 0.5 * sum_t [ p*(1 + log 2\\pi) + log|sC_t| ]
  logdet_sC <- vapply(seq_len(TT), function(t) {
    M <- matrix(sC[, , t, drop = FALSE], nrow = p, ncol = p)
    determinant(M, logarithm = TRUE)$modulus[1]
  }, numeric(1))

  theta_entropy <- 0.5 * sum(p * (1 + log(2*pi)) + logdet_sC)

  list(
    exps = exps,
    vars = vars,
    exps2 = exps2,
    standard.forecast.errors = sfe,
    sm = sm, sC = sC, fm = fm, fC = fC,
    elbo_theta = theta_entropy
  )
}
