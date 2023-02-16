#' Transfer Function exDQLM - ISVB algorithm
#'
#' The function applies an Importance Sampling Variational Bayes (ISVB) algorithm to estimate the posterior of an exDQLM with exponential decay transfer function component.
#'
#' @inheritParams exdqlmISVB
#' @param X A univariate time-series which will be the input of the transfer function component.
#' @param lam Transfer function rate parameter lambda, a value between 0 and 1.
#' @param tf.df Discount factor(s) used for the transfer function component.
#' @param tf.m0 Prior mean of the transfer function component.
#' @param tf.C0 Prior covariance of the transfer function component.
#'
#' @return A list of the following is returned:
#' \itemize{
#'   \item `run.time` - Algorithm run time in seconds.
#'   \item `iter` - Number of iterations until convergence was reached.
#'   \item `dqlm.ind` - Logical value indicating whether gamma was fixed at `0`, reducing the exDQLM to the special case of the DQLM.
#'   \item `model` - List of the augmented state-space model including `GG`, `FF`, prior parameters `m0` and `C0`.
#'   \item `p0` - The quantile which was estimated.
#'   \item `df` - Discount factors used for each block, including transfer function component.
#'   \item `dim.df` - Dimension used for each block of discount factors, including transfer function component.
#'   \item `lam` - Transfer function rate parameter lambda.
#'   \item `sig.init` - Initial value for sigma, or value at which sigma was fixed if `fix.sigma=TRUE`.
#'   \item `seq.sigma` - Sequence of sigma estimated by the algorithm until convergence.
#'   \item `samp.theta` - Posterior sample of the state vector variational distribution.
#'   \item `samp.post.pred` - Sample of the posterior predictive distributions.
#'   \item `map.standard.forecast.errors` - MAP standardized one-step-ahead forecast errors.
#'   \item `samp.sigma` - Posterior sample of scale parameter sigma variational distribution.
#'   \item `samp.vts` - Posterior sample of latent parameters, v_t, variational distributions.
#'   \item `theta.out` - List containing the variational distribution of the state vector including filtered distribution parameters (`fm` and `fC`) and smoothed distribution parameters (`sm` and `sC`).
#'   \item `vts.out` - List containing the variational distributions of latent parameters v_t.
#'   \item `median.kt` - Median number of time steps until the effect of X_t is less than or equal to 1e-3.
#' }
#' If `dqlm.ind=FALSE`, the list also contains:
#' \itemize{
#'   \item `gam.init` - Initial value for gamma, or value at which gamma was fixed if `fix.gamma=TRUE`.
#'   \item `seq.gamma` - Sequence of gamma estimated by the algorithm until convergence.
#'   \item `samp.gamma` - Posterior sample of skewness parameter gamma variational distribution.
#'   \item `samp.sts` - Posterior sample of latent parameters, s_t, variational distributions.
#'   \item `gammasig.out` - List containing the IS estimate of the variational distribution of sigma and gamma.
#'   \item `sts.out` - List containing the variational distributions of latent parameters s_t.
#' }
#' Or if `dqlm.ind=TRUE`, the list also contains:
#' \itemize{
#'   \item `sig.out` - List containing the IS estimate of the variational distribution of sigma.
#'  }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:1095]
#' X = ELIanoms[1:1095]
#' trend.comp = polytrendMod(1,mean(y),10)
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = combineMods(trend.comp,seas.comp)
#' M1 = transfn_exdqlmISVB(y,p0=0.85,model=model,
#'                           X,df=c(1,1),dim.df = c(1,6),
#'                           gam.init=-3.5,sig.init=15,
#'                           lam=0.38,tf.df=c(0.97,0.97))
#' }
#'
transfn_exdqlmISVB<-function(y,p0,model,X,df,dim.df,lam,tf.df,fix.gamma=FALSE,gam.init=NA,fix.sigma=TRUE,sig.init=NA,dqlm.ind=FALSE,
                             exps0,tol=0.1,n.IS=500,n.samp=200,PriorSigma=NULL,PriorGamma=NULL,tf.m0=rep(0,2),tf.C0=diag(1,2),verbose=TRUE){
  # check inputs
  y = check_ts(y)
  X = check_ts(X)
  if(length(X) != length(y)){stop("y and X must be time-series of the same length")}
  model = check_mod(model)
  if(length(lam) != 1 | lam >= 1 | lam <= 0){stop("lam must be a single value between 0 and 1")}
  if(!methods::hasArg(dim.df)){
    if(length(df)!=1){
      stop("length of component discount factors does not match length of component dimensions")
    }
    dim.df = p
  }
  if(length(tf.m0)!=2){
    stop("tf.m0 should have length 2")
  }
  tf.C0 = as.matrix(tf.C0)
  if(any(dim(tf.C0)!=2)){
    stop("tf.C0 should be a 2 by 2 covariance matrix")
  }

  # initialize quantile
  if(methods::hasArg(exps0)){
    TT = length(y)
    if(length(exps0) != TT){stop("exp0 must have same length as y")}
  }else{
    TT = length(y)
    p = length(model$m0)
    if(!is.na(dim(model$GG)[3])){
      if(dim(model$GG)[3] != TT){stop("time-varying dimension of GG does not match length of y")}
    }
    GG = array(model$GG,c(p,p,TT)); model$GG = GG
    if(ncol(model$FF)>1){
      if(ncol(model$FF) != TT){stop("time-varying dimension of FF does not match length of y")}
    }
    FF = matrix(model$FF,p,TT); model$FF = FF
    init.dlm = dlm_df(y,model,df,dim.df,s.priors=list(l0=1,S0=1),just.lik=FALSE)
    exps0 = apply(FF*t(init.dlm$m),2,sum) + stats::qnorm(p0,0,sqrt(init.dlm$s[TT]))
  }

  # augment state-space model
  temp.p = length(model$m0)
  p = temp.p + 2
  FF = matrix(0,p,TT)
  FF[1:temp.p,] = model$FF
  FF[seq(temp.p+1,temp.p+2,2),] = 1
  GG = array(0,c(p,p,TT))
  GG[1:temp.p,1:temp.p,] = model$GG
  GG[(temp.p+2-1):(temp.p+2),(temp.p+2-1):(temp.p+2),] = matrix(c(lam,0,NA,1),2,2)
  GG[(temp.p+2-1),(temp.p+2),] = X

  # update model and dfs with transfer function component
  tf.model <- NULL
  tf.model$GG = GG
  tf.model$FF = FF
  tf.model$m0 = c(model$m0,tf.m0)
  tf.model$C0 = magic::adiag(model$C0,tf.C0)
  tf.model.df = c(df,matrix(tf.df,1,2))
  tf.model.dim.df = c(dim.df,rep(1,2))

  # fit transfer function exdqlm
  tf.return = exdqlmISVB(y,p0,tf.model,tf.model.df,tf.model.dim.df,fix.gamma,gam.init,fix.sigma,sig.init,dqlm.ind,
                         exps0,tol,n.IS,n.samp,PriorSigma,PriorGamma,verbose)
  tf.return$lam = lam

  k_seq = (log(1e-3)-log(abs(c(tf.model$m0[1],tf.return$theta.out$sm[(dim(tf.return$theta.out$sm)[1]-1),-TT])*c(X))))/(log(lam))
  tf.return$median.kt = stats::median(k_seq)

  # return results
  return(tf.return)
}
