#' Create an n-th order polynomial exDQLM component
#'
#' The function creates an n-th order polynomial exDQLM component.
#'
#' @param order The order of the polynomial model.
#' @param m0 Prior mean of the state vector. Default is `m0 = rep(0,order)`.
#' @param C0 Prior covariance of the state vector. Default is `C0 = 1e3*diag(order)`.
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
#' # create a second order polynomial component
#' trend.comp = polytrendMod(2,rep(0,2),10*diag(2))
polytrendMod = function(order,m0,C0){
  GG = diag(rep(1,order))
  FF = rep(0,order)
  if(order > 1){
    for(i in 2:order){
      GG[i-1,i] = 1
    }
  }
  FF[1] = 1
  if(methods::hasArg(m0)){
    if(length(m0) != nrow(GG)){stop("length of m0 does not match specified polynomial component")}
  }else{
    m0 = rep(0,nrow(GG))
  }
  if(methods::hasArg(C0)){
    C0 = as.matrix(C0)
    if((nrow(C0) != nrow(GG)) || (ncol(C0) != nrow(GG))){stop("dimensions of C0 do not match specified polynomial component")}
  }else{
    C0 = 1e3*diag(nrow(GG))
  }
  return(list(FF = FF, GG = GG, m0 = m0, C0 = C0))
}
