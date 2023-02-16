#' Combines state space blocks of an exDQLM
#'
#' The function combines two models into a single state space model for an exDQLM.
#'
#' @param m1 List containing the first model to be combined.
#' @param m2 List containing the second model to be combined.
#'
#' @return List containing the new combined state space model components.
#' @export
#'
#' @examples
#' trend.comp = polytrendMod(2,rep(0,2),10*diag(2))
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = combineMods(trend.comp,seas.comp)
#' # using dlm package
#' library(dlm)
#' model = combineMods(dlmModPoly(order=2,C0=10*diag(2)),dlmModTrig(365,2,C0=10*diag(4)))
#'
combineMods = function(m1,m2){
  if(dlm::is.dlm(m1)){
    m1 = dlmMod(m1)
    message("m1 converted from a dlm object using 'dlmMod(m1)'")
  }
  if(dlm::is.dlm(m2)){
    m2 = dlmMod(m2)
    message("m2 converted from a dlm object using 'dlmMod(m2)'")
  }
  m1 = check_mod(m1)
  m2 = check_mod(m2)
  n = length(c(m1$m0,m2$m0))
  model<- NULL
  if(ncol(m1$FF)>1 | ncol(m2$FF)>1){
    if(ncol(m1$FF)>1 & ncol(m2$FF)>1 & ncol(m1$FF) != ncol(m2$FF)){
      stop("incompatible number of columns in m1$FF and m2$FF")
    }
    model$FF = matrix(0,n,max(ncol(m1$FF),ncol(m2$FF)))
    model$FF[1:nrow(m1$FF),] = m1$FF
    model$FF[(nrow(m1$FF)+1):n,] = m2$FF
  }else{
    model$FF = matrix(c(m1$FF,m2$FF),n,1)
  }
  if(!is.na(dim(m1$GG)[3]) | !is.na(dim(m2$GG)[3])){
    if(!is.na(dim(m1$GG)[3]) & !is.na(dim(m2$GG)[3]) & dim(m1$GG)[3] != dim(m2$GG)[3]){
      stop("incompatible third dimensions of m1$GG and m2$GG")
    }
    model$GG = array(0,c(n,n,max(dim(m1$GG)[3],dim(m2$GG)[3],na.rm = TRUE)))
    model$GG[1:dim(m1$GG)[1],1:dim(m1$GG)[1],] = m1$GG
    model$GG[(dim(m1$GG)[1]+1):n,(dim(m1$GG)[1]+1):n,] = m2$GG
  }else{
    model$GG = magic::adiag(m1$GG,m2$GG)
  }
  model$m0 = c(m1$m0,m2$m0)
  model$C0 = magic::adiag(m1$C0,m2$C0)
  return(model)
}
