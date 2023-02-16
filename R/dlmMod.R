#' Create state space model of exDQLM from DLM
#'
#' The function creates a state space model of an exDQLM from "`dlm`" object.
#'
#' @param m An object of class "`dlm`" representing the DLM version of the desired exDQLM state space model. Only time-invariant `dlm` objects are currently considered.
#'
#' @return List containing only the components of `m` needed for the exDQLM state space model.
#' @export
#'
#' @examples
#' library(dlm)
#' m = dlmModPoly(order=2,C0=10*diag(2)) + dlmModTrig(365,2,C0=10*diag(4))
#' model = dlmMod(m)
#'
dlmMod = function(m){
  if(!dlm::is.dlm(m)){
    stop("input must be a 'dlm' object")
  }
  if(!is.null(m$JFF) | !is.null(m$JGG) |
     !is.null(m$JV) | !is.null(m$JW)){
    stop("input must be a time-invariant 'dlm' object")
  }
  model = as.list(m[c(1,2,3,5)])
  model$FF = t(model$FF)
  return(model = model)
}
