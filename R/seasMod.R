#' Create Fourier representation of a periodic exDQLM component
#'
#' The function creates a Fourier form periodic component for given period and harmonics.
#'
#' @param p The period.
#' @param h Vector of harmonics to be included.
#' @param m0 Prior mean of the state vector.
#' @param C0 Prior covariance of the state vector.
#'
#' @return A list of the following:
#' \itemize{
#'   \item FF - Observational vector.
#'   \item GG - Evolution matrix.
#'   \item m0 - Prior mean of the state vector.
#'   \item C0 - Prior covariance of the state vector.
#' }
#' @export
#'
#' @examples
#' # create a seasonal component with first, second and fourth harmonics of a period of 365
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#'
seasMod = function(p, h, m0, C0){
  nh = length(h)
  w = h * 2 * pi / p
  if( max(w) == pi){
    G = array(0,c(nh-1,2,2))
    for(i in 1:(nh-1)){
      G[i,1,1] =  cos(w[i])
      G[i,1,2] =  sin(w[i])
      G[i,2,1] = -sin(w[i])
      G[i,2,2] =  cos(w[i])
    }
    for(i in 1:(nh-1)){
      if(i ==1 ){ GG = G[1,,]}
      else{GG = magic::adiag(GG,G[i,,])}
    }
    GG = magic::adiag(GG,-1)
    FF = vector("numeric", 2*nh - 1)
    FF[1:(2*nh-1) %% 2 == 1] = 1
  }else{
    G = array(0,c(nh,2,2))
    for(i in 1:nh){
      G[i,1,1] =  cos(w[i])
      G[i,1,2] =  sin(w[i])
      G[i,2,1] = -sin(w[i])
      G[i,2,2] =  cos(w[i])
    }
    for(i in 1:nh){
      if(i ==1 ){ GG = G[1,,]}
      else{GG = magic::adiag(GG,G[i,,])}
    }
    FF = vector("numeric", 2*nh)
    FF[1:(2*nh) %% 2 == 1] = 1
  }
  if(methods::hasArg(m0)){
    if(length(m0) != nrow(GG)){stop("length of m0 does not match specified seasonal component(s)")}
  }else{
    m0 = rep(0,nrow(GG))
  }
  if(methods::hasArg(C0)){
    C0 = as.matrix(C0)
    if((nrow(C0) != nrow(GG)) || (ncol(C0) != nrow(GG))){stop("dimensions of C0 do not match specified seasonal component(s)")}
  }else{
    C0 = 1e3*diag(nrow(GG))
  }
  return(list(FF = FF, GG = GG, m0 = m0, C0 = C0))
}
