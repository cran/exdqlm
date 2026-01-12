#' Plot exDQLM
#'
#' The function plots the MAP estimates and 95% credible intervals (CrIs) of the dynamic quantile of an exDQLM.
#'
#' @param m1 An object of class "\code{exdqlmMCMC}" or "\code{exdqlmISVB}".
#' @param add If `TRUE`, the dynamic quantile will be added to existing plot.
#' @param col Color of dynamic quantile to be plotted. Default is `purple`.
#' @param cr.percent Percentage used in the calculation of the credible intervals.
#'
#' @return A list of the following is returned:
#'  \itemize{
#'   \item `map.quant` - MAP estimate of the dynamic quantile.
#'   \item `lb.quant` - Lower bound of the 95% CrIs of the dynamic quantile.
#'   \item `ub.quant` - Upper bound of the 95% CrIs of the dynamic quantile.
#' }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                    gam.init=-3.5,sig.init=15)
#' exdqlmPlot(M0,col="blue")
#' }
#'
exdqlmPlot <- function(m1,add=FALSE,col="purple",cr.percent=0.95){

  # check inputs
  if(!is.exdqlmMCMC(m1) && !is.exdqlmISVB(m1)){
    stop("m1 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
  }
  y = m1$y
  TT = length(y)
  p = dim(m1$samp.theta)[1]
  n.samp = dim(m1$samp.theta)[3]
  if(cr.percent<=0 | cr.percent>=1){
    stop("cr.percent must be between 0 and 1")
  }
  half.alpha = (1 - cr.percent)/2

  # 95% CrIs
  big_FF = array(m1$model$FF,c(p,TT,n.samp))
  quant.samps = colSums(big_FF*m1$samp.theta)
  map.quant = rowMeans(quant.samps)
  lb.quant = matrixStats::rowQuantiles(quant.samps, probs = half.alpha)
  ub.quant = matrixStats::rowQuantiles(quant.samps, probs = cr.percent + half.alpha)

  # plot
  if(!add){
    stats::plot.ts(y,xlab="time",ylab=sprintf("quantile %s%% CrIs",100*cr.percent),ylim=range(c(y,lb.quant,ub.quant)),col="dark grey")
  }
  ts.xy = grDevices::xy.coords(y)
  graphics::lines(ts.xy$x,map.quant,col=col,lwd=1.5)
  graphics::lines(ts.xy$x,lb.quant,col=col,lwd=0.75,lty=2)
  graphics::lines(ts.xy$x,ub.quant,col=col,lwd=0.75,lty=2)

  ret = list(map.quant=map.quant,lb.quant=lb.quant,ub.quant=ub.quant)
  return(invisible(ret))
}
