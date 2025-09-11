#' k-step-ahead quantile forecasts
#'
#' Computes filtered and \code{k}-step-ahead forecast quantiles from a fitted
#' dynamic quantile model and optionally adds them to an existing plot.
#'
#' @param y A univariate numeric time series (vector or \code{ts}) of the observed response.
#' @param start.t Integer index at which forecasts start (must be within the span of the fitted model in \code{m1}).
#' @param k Integer; number of steps ahead to forecast.
#' @param m1 A fitted exDQLM model object, typically returned by [exdqlmISVB()] or [exdqlmMCMC()].
#' @param fFF Optional state vector(s) for the forecast steps. A numeric matrix with
#'   \eqn{p} rows and either 1 column (non–time-varying) or \code{k} columns (time-varying).
#'   Its dimension must match the fitted model in \code{m1}.
#' @param fGG Optional evolution matrix/matrices for the forecast steps. Either a numeric
#'   \eqn{p \times p} matrix (non–time-varying) or a \eqn{p \times p \times k} array (time-varying).
#'   Its dimensions must match the fitted model in \code{m1}.
#' @param plot Logical; if \code{TRUE}, plot filtered and forecast quantiles with
#'   equal–tailed credible intervals. Default \code{TRUE}.
#' @param add Logical; if \code{TRUE}, add the forecasted quantiles to the current plot.
#'   Default \code{FALSE}.
#' @param cols Character vector of length 2 giving the colors for filtered and forecasted
#'   quantiles respectively. Default \code{c("purple","magenta")}.
#' @param cr.percent Numeric in \code{(0, 1)}; the probability mass for the credible
#'   intervals (e.g., \code{0.95}). Default \code{0.95}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{fa} Forecast state mean vectors (\eqn{p \times k} matrix).
#'   \item \code{fR} Forecast state covariance matrices (\eqn{p \times p \times k} array).
#'   \item \code{ff} Forecast quantile means (length-\code{k} numeric).
#'   \item \code{fQ} Forecast quantile variances (length-\code{k} numeric).
#' }
#'
#' @examples
#' \donttest{
#'  # Toy example; keep small and fast
#'  y <- scIVTmag[1:100]
#'  model <- polytrendMod(1, stats::quantile(y, 0.85), 10)
#'  M0 <- exdqlmISVB(y, p0 = 0.85, model, df = c(0.98), dim.df = c(1),
#'                   gam.init = -3.5, sig.init = 15)
#'  exdqlmForecast(y, start.t = 90, k = 10, m1 = M0)
#' }
#'
#' @export

exdqlmForecast = function(y,start.t,k,m1,fFF=NULL,fGG=NULL,plot=TRUE,add=FALSE,cols=c("purple","magenta"),cr.percent=0.95){

  # check inputs
  check_ts(y)
  p = dim(m1$model$GG)[1]
  TT = dim(m1$model$GG)[3]
  if(!is.exdqlm(m1)){
    stop("m1 must be an output from 'exdqlmISVB()' or 'exdqlmMCMC()'")
  }
  if(cr.percent<=0 | cr.percent>=1){
    stop("cr.percent must be between 0 and 1")
  }
  half.alpha = (1 - cr.percent)/2
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
    qlb = qmap + sapply(1:start.t,function(t){stats::qnorm(half.alpha,0,sqrt(t(m1$model$FF[,t])%*%m1$theta.out$fC[,,t]%*%m1$model$FF[,t]))})
    qub = qmap + sapply(1:start.t,function(t){stats::qnorm(cr.percent + half.alpha,0,sqrt(t(m1$model$FF[,t])%*%m1$theta.out$fC[,,t]%*%m1$model$FF[,t]))})
    # forecast estimates
    fqlb = ff+stats::qnorm(half.alpha,0,sqrt(fQ))
    fqub = ff+stats::qnorm(cr.percent + half.alpha,0,sqrt(fQ))
    # filtered and forecasted quantiles & CrIs
    ts.xy = grDevices::xy.coords(y)
    if(!add){
      stats::plot.ts(y,xlim=c(ts.xy$x[start.t]-2*k*diff(ts.xy$x)[1],ts.xy$x[start.t]+k*diff(ts.xy$x)[1]),ylim=range(c(y,qlb,qub,fqlb,fqub)),type="l",ylab="quantile forecast",col="dark grey",xlab="time")
    }
    graphics::lines(ts.xy$x[1:start.t],qlb,col=cols[1],lty=3)
    graphics::lines(ts.xy$x[1:start.t],qub,col=cols[1],lty=3)
    graphics::lines(ts.xy$x[1:start.t],qmap,col=cols[1],lwd=1.5)
    graphics::lines(seq(from = ts.xy$x[start.t], by = diff(ts.xy$x)[1], length.out = k+1),c(qmap[start.t],ff),col=cols[2])
    graphics::lines(seq(from = ts.xy$x[start.t], by = diff(ts.xy$x)[1], length.out = k+1),c(qub[start.t],fqub),col=cols[2],lty=3)
    graphics::lines(seq(from = ts.xy$x[start.t], by = diff(ts.xy$x)[1], length.out = k+1),c(qlb[start.t],fqlb),col=cols[2],lty=3)
  }

  # return forecast distributions
  return(invisible(list(fa=fa,fR=fR,ff=ff,fQ=fQ)))
}

