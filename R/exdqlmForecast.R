#' k-step-ahead Forecast
#'
#' The function estimates and plots the k-step-ahead forecasted quantile distribution from the filtered quantile estimates.
#'
#' @inheritParams exdqlmISVB
#' @param start.t Time index at which to start the forecast.
#' @param k Number of k-steps-ahead to forecast.
#' @param m1 An object of class "`exdqlm`".
#' @param fFF State vector for the forecast steps. `fFF` must have either 1 (non-time-varying) or k (time-varying) columns. The dimension of `fFF` must match the estimated exdqlm in `m1`.
#' @param fGG Evolution matrix for the forecast steps. `fGG` must be either a matrix (non-time-varying) or an array of depth k (time-varying). The dimensions of `fGG` must match the estimated exdqlm in `m1`.
#' @param plot If `TRUE` the forecasted quantile estimates and 95% credible intervals are plotted, along with the filtered quantile estimates and 95% credible intervals for reference. Default is `TRUE`.
#' @param add If `TRUE`, the forecasted quantile will be added to the existing plot. Default is `FALSE`.
#' @param cols Two colors used to plot filtered and forecasted quantile estimates respectively. Default is `c("purple","magenta")`.
#'
#' @return A list containing the following is returned:
#'  \itemize{
#'  \item `fa` - The forecasted state mean vectors.
#'  \item `fR` - The forecasted state covariance matrices.
#'  \item `ff` - The forecasted quantile mean estimates.
#'  \item `fQ` - The forecasted quantile variances.
#'  }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                    gam.init=-3.5,sig.init=15)
#' exdqlmForecast(y,start.t=90,k=10,M0)
#' }
#'
exdqlmForecast = function(y,start.t,k,m1,fFF=NULL,fGG=NULL,plot=TRUE,add=FALSE,cols=c("purple","magenta")){

  # check inputs
  check_ts(y)
  p = dim(m1$model$GG)[1]
  TT = dim(m1$model$GG)[3]
  if(class(m1) != c("exdqlm")){
    stop("m1 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
  }
  if(is.null(fFF)){
     if(TT-start.t < k){ stop("fFF and fGG must be provided for forecasts extending past the length of the estimated exdqlm")}
     fFF = m1$model$FF[,(start.t+1):(start.t+k)]
     fGG = m1$model$GG[,,(start.t+1):(start.t+k)]
  }else{
    fFF = as.matrix(fFF)
    if(nrow(fFF) != p){ stop("dimension of fFF must match the estimated exdqlm") }
    if(!any(ncol(fFF) == c(1,k))){ stop("fFF must have either 1 (non-time-varying) or k (time-varying) columns")}
    fGG = as.array(fGG)
    if(any(dim(fGG)[1:2] != p)){ stop("dimension of fGG must match the estimated exdqlm") }
    if(!is.na(dim(fGG)[3])){
      if(dim(fGG)[3] != k){
        stop("fGG must be either a matrix (non-time-varying) or an array of depth k (time-varying)")
        }
      }
  }
  fFF = matrix(fFF,p,k)
  fGG = array(fGG,c(p,p,TT))

  #### forecast k steps
  df.mat = make_df_mat(m1$df,m1$dim.df,p)
  fm = m1$theta.out$fm[,start.t]
  fC = m1$theta.out$fC[,,start.t]
  fa = matrix(NA,p,k)
  fR = array(NA,c(p,p,k))
  ff = rep(NA,p)
  fQ = rep(NA,p)
  for(i in 1:k){
    if(i == 1){
      fa[,1] = fGG[,,i]%*%fm
      fR[,,1] = fGG[,,i]%*%fC%*%t(fGG[,,i]) + df.mat*fC
      ff[1] = t(fFF[,i])%*%fa[,1]
      fQ[1] = t(fFF[,i])%*%fR[,,1]%*%fFF[,i]
    }else{
      fa[,i] = fGG[,,i]%*%fa[,(i-1)]
      fR[,,i] = fGG[,,i]%*%fR[,,(i-1)]%*%t(fGG[,,i]) + df.mat*fR[,,(i-1)]
      ff[i] = t(fFF[,i])%*%fa[,i]
      fQ[i] = t(fFF[,i])%*%fR[,,i]%*%fFF[,i]
    }
  }

  # plot forecast
  if(plot){
    # filtered estimate for reference
    qmap = apply(matrix(m1$model$FF[,1:start.t]*m1$theta.out$fm[,1:start.t],p,start.t),2,sum)
    qlb = qmap + sapply(1:start.t,function(t){stats::qnorm(0.025,0,sqrt(t(m1$model$FF[,t])%*%m1$theta.out$fC[,,t]%*%m1$model$FF[,t]))})
    qub = qmap + sapply(1:start.t,function(t){stats::qnorm(0.975,0,sqrt(t(m1$model$FF[,t])%*%m1$theta.out$fC[,,t]%*%m1$model$FF[,t]))})
    # forecast estimates
    fqlb = ff+stats::qnorm(0.025,0,sqrt(fQ))
    fqub = ff+stats::qnorm(0.975,0,sqrt(fQ))
    # filtered and forecasted quantiles & CrIs
    ts.xy = grDevices::xy.coords(y)
    if(!add){
      stats::plot.ts(y,xlim=c(ts.xy$x[start.t]-2*k,ts.xy$x[start.t]+k),ylim=range(c(y,qlb,qub,fqlb,fqub)),type="l",ylab="quantile forecast",col="dark grey",xlab="time")
    }
    graphics::lines(ts.xy$x[1:start.t],qlb,col=cols[1],lty=3)
    graphics::lines(ts.xy$x[1:start.t],qub,col=cols[1],lty=3)
    graphics::lines(ts.xy$x[1:start.t],qmap,col=cols[1],lwd=1.5)
    graphics::lines(ts.xy$x[start.t]:(ts.xy$x[start.t]+k),c(qmap[start.t],ff),col=cols[2])
    graphics::lines(ts.xy$x[start.t]:(ts.xy$x[start.t]+k),c(qub[start.t],fqub),col=cols[2],lty=3)
    graphics::lines(ts.xy$x[start.t]:(ts.xy$x[start.t]+k),c(qlb[start.t],fqlb),col=cols[2],lty=3)
  }

  # return forecast distributions
  return(invisible(list(fa=fa,fR=fR,ff=ff,fQ=fQ)))
}

