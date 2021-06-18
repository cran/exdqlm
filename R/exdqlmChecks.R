#' exDQLM Diagnostics
#'
#' The function computes the following for the model(s) provided: the posterior predictive loss criterion based off the check loss, the one-step-ahead distribution sequence and its KL divergence from normality. The function also plots the following: the qq-plot and ACF plot corresponding to the one-step-ahead distribution sequence, and a time series plot of the MAP standard forecast errors.
#'
#' @inheritParams exdqlmPlot
#' @param m2 An optional additional object of class "`exdqlm`" to compare with `m1`.
#' @param plot If `TRUE`, the following will be plotted for `m1` and `m2` (if provided): a qq-plot and ACF plot of the MAP one-step-ahead distribution sequence, and a time series plot of the standardized forecast errors.
#' @param cols Color(s) used to plot diagnostics.
#' @param ref Reference sample of size `length(y)` from a standard normal distribution used to compute the KL divergence.
#'
#' @return A list containing the following is returned:
#'  \itemize{
#'  \item `m1.uts` - The one-step-ahead distribution sequence of `m1`.
#'  \item `m1.KL` - The KL divergence of `m1.uts` and a standard normal.
#'  \item `m1.pplc` - The posterior predictive loss criterion of `m1` based off the check loss function.
#'  \item `m1.qq` - The ordered pairs of the qq-plot comparing `m1.uts` with a standard normal distribution.
#'  \item `m1.acf` - The autocorrelations of `m1.uts` by lag.
#'  }
#'  If `m2` is provided, analogous results for `m2` are also included in the list.
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,mean(y),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.95),dim.df = c(1),
#'                   gam.init=-3.5,sig.init=15)
#' check.out = exdqlmChecks(y,M0,plot=FALSE)
#' check.out$m1.KL
#' check.out$m1
#' }
#'
exdqlmChecks <- function(y,m1,m2=NULL,plot=TRUE,cols=c("grey","grey"),ref=NULL){

  # check inputs
  check_ts(y)
  TT = length(y)
  if(class(m1) != c("exdqlm")){
    stop("m1 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
  }
  if(dim(m1$samp.theta)[2] != TT){
    stop("length of dynamic quantile in m1 does not match data")
  }
  cols = c(matrix(cols,2,1))

  ### m1
  # m1 seq.uts
  m1.uts = stats::pnorm(m1$map.standard.forecast.errors)
  # m1 KL divergence
  TT = length(m1$map.standard.forecast.errors)
  if(is.null(ref)){
    ref = stats::rnorm(TT)
  }else{
    ref = c(ref)
    if(length(ref)!=TT){
      stop("ref should be a sample of size 'length(y)' from a standard normal distribution")
    }
  }
  m1.KL = mean(FNN::KL.divergence(ref,m1$map.standard.forecast.errors))
  # m1 pplc
  m1.loss = matrix(NA,TT,dim(m1$samp.post.pred)[2])
  for(t in 1:TT){m1.loss[t,] = CheckLossFn(m1$p0,y[t]-m1$samp.post.pred[t,])}
  m1.pplc = sum(rowMeans(m1.loss))
  # m1 qqplot
  m1.qq = stats::qqnorm(m1$map.standard.forecast.errors,plot=FALSE)
  # m1 acf
  m1.acf = stats::acf(m1.uts,plot=FALSE)
  #
  retlist = list(m1.uts=m1.uts,m1.KL=m1.KL,m1.pplc=m1.pplc,m1.qq=m1.qq,m1.acf=m1.acf)

  ### m2
  if(!is.null(m2)){
    # check inputs
    if(class(m2) != c("exdqlm")){
      stop("m2 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
    }
    if(dim(m1$samp.theta)[2] != TT){
      stop("length of dynamic quantile in m2 does not match data")
    }
    if(TT != length(m2$map.standard.forecast.errors)){
      stop("length of m1 quantile does not match length of m2 quantile")
    }
    if(m1$p0 != m2$p0){
      stop("quantiles estimated in m1 and m2 do not match")
    }
    # m2 seq.uts
    m2.uts = stats::pnorm(m2$map.standard.forecast.errors)
    retlist[["m2.uts"]] = m2.uts
    # m2 KL divergence
    retlist[["m2.KL"]] = mean(FNN::KL.divergence(ref,m2$map.standard.forecast.errors))
    # m2 pplc
    m2.loss = matrix(NA,TT,dim(m2$samp.post.pred)[2])
    for(t in 1:TT){m2.loss[t,] = CheckLossFn(m2$p0,y[t]-m2$samp.post.pred[t,])}
    retlist[["m2.pplc"]] = sum(rowMeans(m2.loss))
    # m2 qqplot
    retlist[["m2.qq"]] = stats::qqnorm(m2$map.standard.forecast.errors,plot=FALSE)
    # m2 acf
    retlist[["m2.acf"]] = stats::acf(m2.uts,plot=FALSE)
  }

  if(plot){
    if(is.null(m2)){
      qq.x.range = range(m1.qq$x)
      qq.y.range = range(m1.qq$y)
      acf.y.range = range(m1.acf$acf)
      fe.y.range = range(m1$map.standard.forecast.errors)
    }else{
      qq.x.range = range(c(m1.qq$x,retlist$m2.qq$x))
      qq.y.range = range(c(m1.qq$y,retlist$m2.qq$y))
      acf.y.range = range(c(m1.acf$acf,retlist$m2.acf$acf))
      fe.y.range = range(c(m1$map.standard.forecast.errors,m2$map.standard.forecast.errors))
    }
    # m1 qqplot
    plot(m1.qq,main="",col=cols[1],pch=20,xlab="Theoretical Quantiles",ylab="M1 Sample Quantiles",xlim=qq.x.range,ylim=qq.y.range)
    graphics::abline(a=0,b=1)
    # m1 acf
    plot(m1.acf,ylab="M1 ACF",col=cols[1],main="",ylim=acf.y.range)
    # m1 forecast errors
    ts.xy = grDevices::xy.coords(y)
    graphics::plot(ts.xy$x,m1$map.standard.forecast.errors,ylab="M1 standard forecast errors",xlab="time",col=cols[1],pch=20,type="l",ylim=fe.y.range)
    graphics::abline(h=0,lty=2)
    ### m2
    if(!is.null(m2)){
      # m2 qqplot
      plot(retlist[["m2.qq"]],main="",col=cols[2],pch=20,xlab="Theoretical Quantiles",ylab="M2 Sample Quantiles",xlim=qq.x.range,ylim=qq.y.range)
      graphics::abline(a=0,b=1)
      # m2 acf
      plot(retlist[["m2.acf"]],ylab="M2 ACF",col=cols[2],main="",ylim=acf.y.range)
      # m2 forecast errors
      graphics::plot(ts.xy$x,m2$map.standard.forecast.errors,ylab="M2 standard forecast errors", xlab="time",col=cols[2],pch=20,type="l",ylim=fe.y.range)
      graphics::abline(h=0,lty=2)
    }
  }

  # return model checks
  return(invisible(retlist))
}
