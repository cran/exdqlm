#' Plot a component of an exDQLM
#'
#' The function plots the dynamic MAP estimates and 95% credible intervals (CrIs) of a specified component of an exDQLM. Alternatively, if `just.theta=TRUE` the MAP estimates and 95% credible intervals (CrIs) of a single element of the dynamic state vector are plotted.
#'
#' @inheritParams exdqlmPlot
#' @param index Index of the component or element of the state vector to be plotted.
#' @param add If `TRUE`, the dynamic component will be added to existing plot.
#' @param col Color of dynamic component to be plotted. Default is `purple`.
#' @param just.theta If `TRUE`, the function plots the dynamic distribution of the `index` element of the state vector. If `just.theta=TRUE`, `index` must have length 1.
#' @param cr.percent Percentage used in the calculation of the credible intervals.
#'
#' @return A list of the following is returned:
#'  \itemize{
#'   \item `map.comp` - MAP estimate of the dynamic component (or element of the state vector).
#'   \item `lb.comp` - Lower bound of the 95% CrIs of the dynamic component (or element of the state vector).
#'   \item `ub.comp` - Upper bound of the 95% CrIs of the dynamic component (or element of the state vector).
#' }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:365]
#' trend.comp = polytrendMod(2,rep(0,2),10*diag(2))
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = combineMods(trend.comp,seas.comp)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98,1),dim.df = c(2,6),
#'                    gam.init=-3.5,sig.init=15,tol=0.05)
#' # plot first harmonic component
#' compPlot(y,M0,index=c(3,4),col="blue")
#' }
#'
compPlot <- function(y, m1, index, add = FALSE, col="purple", just.theta = FALSE, cr.percent = 0.95){

  # check inputs
  check_ts(y)
  if(!is.exdqlm(m1)){
    stop("m1 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
  }
  TT = length(y)
  if(dim(m1$samp.theta)[2] != TT){
    stop("length of dynamic quantile does not match data")
  }
  p = length(index)
  n.samp = dim(m1$samp.theta)[3]
  if(just.theta & p != 1){
    stop("when 'just.theta = TRUE', 'index' should have length 1")
  }
  if(cr.percent<=0 | cr.percent>=1){
    stop("cr.percent must be between 0 and 1")
  }
  half.alpha = (1 - cr.percent)/2

  # 95% CrIs
  if(!just.theta){
      quant.samps = apply(array(m1$model$FF[index,],c(p,TT,n.samp))*array(m1$samp.theta[index,,],c(p,TT,n.samp)),c(2,3),sum)
  }else{
      quant.samps = matrix(m1$samp.theta[index,,],TT,n.samp)
  }
  map.quant = rowMeans(quant.samps)
  lb.quant = apply(quant.samps,1,stats::quantile,probs=half.alpha)
  ub.quant = apply(quant.samps,1,stats::quantile,probs=cr.percent + half.alpha)

  # plot
  if(!add){
    stats::plot.ts(y,xlab="time",ylab=sprintf("%s%% CrIs",100*cr.percent),ylim=range(c(lb.quant,ub.quant)),col="white")
  }
  ts.xy = grDevices::xy.coords(y)
  graphics::lines(ts.xy$x,map.quant,col=col,lwd=1.5)
  graphics::lines(ts.xy$x,lb.quant,col=col,lwd=0.75,lty=2)
  graphics::lines(ts.xy$x,ub.quant,col=col,lwd=0.75,lty=2)

  ret = list(map.comp=map.quant,lb.comp=lb.quant,ub.comp=ub.quant)
  return(invisible(ret))
}
