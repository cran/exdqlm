#' exDQLM - ISVB algorithm
#'
#' The function applies an Importance Sampling Variational Bayes (ISVB) algorithm to estimate the posterior of an exDQLM.
#'
#' @param y A univariate time-series.
#' @param p0 The quantile of interest, a value between 0 and 1.
#' @param model List of the state-space model including `GG`, `FF`, prior parameters `m0` and `C0`.
#' @param df Discount factors for each block.
#' @param dim.df Dimension of each block of discount factors.
#' @param fix.gamma Logical value indicating whether to fix gamma at `gam.init`. Default is `FALSE`.
#' @param gam.init Initial value for gamma (skewness parameter), or value at which gamma will be fixed if `fix.gamma=TRUE`.
#' @param fix.sigma Logical value indicating whether to fix sigma at `sig.init`. Default is `TRUE`.
#' @param sig.init Initial value for sigma (scale parameter), or value at which sigma will be fixed if `fix.sigma=TRUE`.
#' @param dqlm.ind Logical value indicating whether to fix gamma at `0`, reducing the exDQLM to the special case of the DQLM. Default is `FALSE`.
#' @param exps0 Initial value for dynamic quantile. If `exps0` is not specified, it is set to the DLM estimate of the `p0` quantile.
#' @param tol Tolerance for convergence of dynamic quantile estimates. Default is `tol=0.1`.
#' @param n.IS Number of particles for the importance sampling of joint variational distribution of sigma and gamma. Default is `n.IS=500`.
#' @param n.samp Number of samples to draw from the approximated posterior distribution. Default is `n.samp=200`.
#' @param PriorSigma List of parameters for inverse gamma prior on sigma; shape `a_sig` and scale `b_sig`. Default is an inverse gamma with mean 1 (or `sig.init` if provided) and variance 10.
#' @param PriorGamma List of parameters for truncated student-t prior on gamma; center `m_gam`, scale `s_gam` and degrees of freedom `df_gam`. Default is a standard student-t with 1 degree of freedom, truncated to the support of gamma.
#' @param verbose Logical value indicating whether progress should be displayed.
#'
#' @return A object of class "\code{exdqlmISVB}" containing the following:
#' \itemize{
#'    \item `y` - Time-series data used to fit the model.
#'   \item `run.time` - Algorithm run time in seconds.
#'   \item `iter` - Number of iterations until convergence was reached.
#'   \item `dqlm.ind` - Logical value indicating whether gamma was fixed at `0`, reducing the exDQLM to the special case of the DQLM.
#'   \item `model` - List of the state-space model including `GG`, `FF`, prior parameters `m0` and `C0`.
#'   \item `p0` - The quantile which was estimated.
#'   \item `df` - Discount factors used for each block.
#'   \item `dim.df` - Dimension used for each block of discount factors.
#'   \item `sig.init` - Initial value for sigma, or value at which sigma was fixed if `fix.sigma=TRUE`.
#'   \item `seq.sigma` - Sequence of sigma estimated by the algorithm until convergence.
#'   \item `samp.theta` - Posterior sample of the state vector variational distribution.
#'   \item `samp.post.pred` - Sample of the posterior predictive distributions.
#'   \item `map.standard.forecast.errors` - MAP standardized one-step-ahead forecast errors.
#'   \item `samp.sigma` - Posterior sample of scale parameter sigma variational distribution.
#'   \item `samp.vts` - Posterior sample of latent parameters, v_t, variational distributions.
#'   \item `theta.out` - List containing the variational distribution of the state vector including filtered distribution parameters (`fm` and `fC`) and smoothed distribution parameters (`sm` and `sC`).
#'   \item `vts.out` - List containing the variational distributions of latent parameters v_t.
#'   \item `fix.sigma` Logical value indicating whether sigma was fixed at `sig.init`.
#' }
#' If `dqlm.ind=FALSE`, the object also contains:
#' \itemize{
#'   \item `gam.init` - Initial value for gamma, or value at which gamma was fixed if `fix.gamma=TRUE`.
#'   \item `seq.gamma` - Sequence of gamma estimated by the algorithm until convergence.
#'   \item `samp.gamma` - Posterior sample of skewness parameter gamma variational distribution.
#'   \item `samp.sts` - Posterior sample of latent parameters, s_t, variational distributions.
#'   \item `gammasig.out` - List containing the IS estimate of the variational distribution of sigma and gamma.
#'   \item `sts.out` - List containing the variational distributions of latent parameters s_t.
#'   \item `fix.gamma` Logical value indicating whether gamma was fixed at `gam.init`.
#' }
#' Or if `dqlm.ind=TRUE`, the object also contains:
#'  \itemize{
#'  \item `sig.out` - List containing the IS estimate of the variational distribution of sigma.
#'  }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:1095]
#' trend.comp = polytrendMod(1,mean(y),10)
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = trend.comp + seas.comp
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(1,1),dim.df = c(1,6),
#'                  gam.init=-3.5,sig.init=15,tol=0.05)
#' }
#'
exdqlmISVB<-function(y,p0,model,df,dim.df,fix.gamma=FALSE,gam.init=NA,fix.sigma=TRUE,sig.init=NA,dqlm.ind=FALSE,
                     exps0,tol=0.1,n.IS=500,n.samp=200,PriorSigma=NULL,PriorGamma=NULL,verbose=TRUE){

  # check inputs
  y = check_ts(y)
  model = check_mod(model)
  rv = check_logics(gam.init,sig.init,fix.gamma,fix.sigma,dqlm.ind)
  gam.init = rv$gam.init
  dqlm.int = rv$dqlm.ind
  fix.gamma = rv$fix.gamma

  ### Define L and U
  L = L.fn(p0); U = U.fn(p0)
  if(!is.na(gam.init)){
    if(gam.init < L | gam.init > U){
      stop(sprintf("gam.init must be between %s and %s for %s quantile",round(L,3),round(U,3),p0))
    }
  }

  ### sigma and gamma priors
  # sigma ~ IG(a_sig,b_sig)
  if(is.null(PriorSigma)){
    m_sigma = 1
    v_sigma = 10
    PriorSigma$a_sig = (m_sigma^2)/(v_sigma) + 2
    PriorSigma$b_sig = (m_sigma^3)/(v_sigma) + m_sigma
  }else{
    if(!is.list(PriorSigma) | any( is.na( match(c("a_sig", "b_sig"),names(PriorSigma)) ) )){
      stop("`PriorSigma` must be a list containing `a_sig` and `b_sig`")
      }
  }
  # gamma ~ truncated student t on L,U
  if(is.null(PriorGamma)){
    PriorGamma$m_gam = 0
    PriorGamma$s_gam = 1
    PriorGamma$df_gam = 1
   }else{
     if(!is.list(PriorGamma) | any( is.na( match(c("m_gam", "s_gam", "df_gam"),names(PriorGamma)) ) )){
       stop("`PriorGamma` must be a list containing `m_gam`,`s_gam`, and `df_gam`")
     }
   }
  PriorGammaDens<-function(gamma){ crch::dtt(gamma,location = PriorGamma$m_gam, scale = PriorGamma$s_gam, df = PriorGamma$df_gam, left = L, right = U, log = FALSE) }

  ### state-space model
  ## prior, theta ~ N(m0,C0)
  m0 = model$m0
  C0 = model$C0
  #
  TT = length(y)
  p = length(m0)
  if(!is.na(dim(model$GG)[3])){
    if(dim(model$GG)[3] != TT){stop("time-varying dimension of GG does not match length of y")}
  }
  GG = array(model$GG,c(p,p,TT)); model$GG = GG
  if(ncol(model$FF)>1){
    if(ncol(model$FF) != TT){stop("time-varying dimension of FF does not match length of y")}
  }
  FF = matrix(model$FF,p,TT); model$FF = FF
  ## discount factor blocking
  if(!methods::hasArg(dim.df)){
    if(length(df)!=1){
      stop("length of component discount factors does not match length of component dimensions")
    }
    dim.df = p
  }
  df.mat = make_df_mat(df,dim.df,p)

  ### Initialize VB
  gam0 = ifelse(!is.na(gam.init),gam.init,(L+U)/2)
  sig0 = ifelse(!is.na(sig.init),sig.init,1)
  new.gamsig.out = list(E.gam=gam0,V.gam=10,
                        E.sigma=ifelse(!is.na(sig0),sig0,m_sigma),V.sig=10,
                        E.inv.sigma=ifelse(!is.na(sig0),1/sig0,1/m_sigma),
                        E.c2.invb.absgam2.sigma = sig0*(C.fn(p0,gam0)^2)*(abs(gam0)^2)/B.fn(p0,gam0),
                        E.c.invb.absgam = C.fn(p0,gam0)*abs(gam0)/B.fn(p0,gam0),
                        E.c.a.invb.absgam = C.fn(p0,gam0)*A.fn(p0,gam0)*abs(gam0)/B.fn(p0,gam0),
                        E.a2.invb.inv.sigma = (A.fn(p0,gam0)^2)/(B.fn(p0,gam0)*sig0),
                        E.invb.inv.sigma = 1/(sig0*B.fn(p0,gam0)),
                        E.a.invb.inv.sigma = A.fn(p0,gam0)/(B.fn(p0,gam0)*sig0))
  new.sts.out = list(E.sts=rep(truncnorm::etruncnorm(a=0,b=Inf,mean=1,sd=1),TT),
                     E.sts2=rep(truncnorm::etruncnorm(a=0,b=Inf,mean=1,sd=1)^2+truncnorm::vtruncnorm(a=0,b=Inf,mean=1,sd=1),TT))
  new.uts.out = list(E.uts=rep(1/sig0,TT),
                     E.inv.uts=rep(sig0,TT))
  if(methods::hasArg(exps0)){
    if(length(exps0) != TT){ stop("exps0 must have same length as y") }
  }else{
    init.dlm = dlm_df(y,model,df,dim.df,s.priors=list(l0=1,S0=sig0),just.lik=FALSE)
    exps0 = apply(FF*t(init.dlm$m),2,sum) + stats::qnorm(p0,0,sqrt(init.dlm$s[TT]))
  }
  new.theta.out = list(exps=exps0,exps2=exps0^2)

  ### initialize convergence evaluations
  iter = 0
  conv.count = 0
  new.max = Inf
  seq.gamma = new.gamsig.out$E.gam
  seq.sigma = new.gamsig.out$E.sigma

  # function update q(st)
  update_sts<-function(exps,inv.uts,c2.invb.absgam2.sigma,c.invb.absgam,c.a.invb.absgam){
    s.sig2<-1/(1+c2.invb.absgam2.sigma*inv.uts); s.sig = sqrt(s.sig2)
    s.mu<-s.sig2*(c.invb.absgam*(y-exps)*inv.uts-c.a.invb.absgam)
    #
    E.sts = truncnorm::etruncnorm(a=rep(0,TT),b=rep(Inf,TT),mean=s.mu,sd=s.sig)
    V.sts = truncnorm::vtruncnorm(a=rep(0,TT),b=rep(Inf,TT),mean=s.mu,sd=s.sig)
    E.sts2 = s.mu^2 + s.sig2 + s.mu*s.sig*exp(stats::dnorm(-s.mu/s.sig,log = TRUE)-stats::pnorm(s.mu/s.sig,log.p = TRUE))
    return(list(sts.sig2=s.sig2,sts.mu=s.mu,
                E.sts=E.sts,E.sts2=E.sts2))
  }

  # function update q(ut)
  update_uts<-function(exps,exps2,sts,sts2,inv.sigma,a2.invb.inv.sigma,invb.inv.sigma,c.invb.absgam,c2.invb.absgam2.sigma){
    u.lambda = 0.5
    u.psi = (a2.invb.inv.sigma + 2*inv.sigma)
    u.chi = invb.inv.sigma*(y^2-2*y*exps+exps2) - 2*c.invb.absgam*sts*(y-exps) + c2.invb.absgam2.sigma*sts2
    u.chi[u.chi<=0] = 1e-3
    #
    E.uts = sapply(u.chi,function(x){sqrt(x/u.psi)*HyperbolicDist::besselRatio(sqrt(x*u.psi),u.lambda,1,Inf)})
    E.inv.uts = sapply(u.chi,function(x){sqrt(u.psi/x)*HyperbolicDist::besselRatio(sqrt(x*u.psi),u.lambda,1,Inf)-2*u.lambda/x})
    return(list(uts.lambda=u.lambda,uts.psi=u.psi,uts.chi=u.chi,E.uts=E.uts,E.inv.uts=E.inv.uts))
  }

  # function update q(theta) ffbsm
  update_theta<-function(ex.f,ex.q){
    # initialize ffbs
    m <- sm <- matrix(NA,p,TT)
    C <- sC <- array(NA,c(p,p,TT))
    standard.forecast.errors <- rep(NA,TT)
    ## forward filter
    # first iteration
    a = as.vector(GG[,,1]%*%m0)
    P = GG[,,1]%*%C0%*%t(GG[,,1])
    R = P + df.mat*P
    R = (R + t(R))/2
    f = t(FF[,1])%*%a + ex.f[1]
    q = t(FF[,1])%*%R%*%FF[,1]  + ex.q[1]
    m[,1] = a + t(R)%*%FF[,1]%*%(y[1]-f)/q[1]
    C[,,1] = R - t(R)%*%FF[,1]%*%t(FF[,1])%*%R/q[1]
    C[,,1] = (C[,,1] + t(C[,,1]))/2
    standard.forecast.errors[1] = (y[1]-f)/sqrt(q)
    # t = 2:TT
    for(t in 2:TT){
      a = as.vector(GG[,,t]%*%m[,(t-1)])
      P = GG[,,t]%*%C[,,(t-1)]%*%t(GG[,,t])
      R = P + df.mat*P
      R = (R + t(R))/2
      f = t(FF[,t])%*%a + ex.f[t]
      fB = t(FF[,t])%*%R
      q = fB%*%FF[,t] + ex.q[t]
      m[,t] = a + t(fB)%*%(y[t]-f)/q[1]
      C[,,t] = R - t(fB)%*%fB/q[1]
      C[,,t] = (C[,,t] + t(C[,,t]))/2
      standard.forecast.errors[t] = (y[t]-f)/sqrt(q)
    }
    ## backwards smoothing
    sC[,,TT] = C[,,TT]
    sm[,TT] = m[,TT]
    for(t in (TT-1):1){
      P = GG[,,(t+1)]%*%C[,,(t)]%*%t(GG[,,(t+1)])
      R = P + df.mat*P
      R = (R + t(R))/2
      svd.R = svd(R)
      inv.R = svd.R$u%*%diag(1/svd.R$d,p)%*%t(svd.R$u)
      sB = C[,,t]%*%t(GG[,,t])%*%inv.R
      sm[,t] = m[,t] + sB%*%(sm[,(t+1)]-as.vector(GG[,,(t+1)]%*%m[,(t)]))
      sC[,,t] = C[,,t] + sB%*%(sC[,,(t+1)]-R)%*%t(sB)
      sC[,,t] = (sC[,,t]+t(sC[,,t]))/2
    }
    exps =  apply(FF*sm,2,sum)
    vars = c(apply(matrix(1:TT,TT,1),1,function(x){t(FF[,x])%*%sC[,,x]%*%FF[,x]}))
    exps2 = exps^2 + vars
    return(list(exps=exps,vars=vars,exps2=exps2,standard.forecast.errors=standard.forecast.errors,sm=sm,sC=sC,fm=m,fC=C))
  }

  # function approximate q(sigma,gamma) with importance sampling
  update_gamma_sigma<-function(gamma,var.gam,sigma,var.sig,exps,exps2,sts,sts2,uts,inv.uts){
    gam.sig2 = max(var.gam,0.001)
    gam.mu = gamma
    v_sig = max(var.sig,0.001)
    m_sig = sigma
    # sampling distribution functions
    rr_gamma = function(n){
      return(crch::rtt(n, location = gam.mu, scale = sqrt(gam.sig2), df = 1, left = L+1e-3, right = U-1e-3))
    }
    dr_gamma = function(gam,log.ind=FALSE){
      return(crch::dtt(gam, location = gam.mu, scale = sqrt(gam.sig2), df = 1, left = L+1e-3, right = U-1e-3,log = log.ind))
    }
    rr_sigma = function(n){
      return(crch::rtt(n, location = m_sig, scale = sqrt(v_sig), df = 1, left = 0, right = Inf))
    }
    dr_sigma = function(sig,log.ind=FALSE){
      return(crch::dtt(sig, location = m_sig, scale = sqrt(v_sig), df = 1, left = 0, right = Inf, log = log.ind))
    }
    # variational function q(gamma,sigma) up to proportionality constant
    dq = function(sig,gam,log.ind=FALSE){
      a = A.fn(p0,gam); b = B.fn(p0,gam); c = C.fn(p0,gam);
      if(log.ind==FALSE){
        q.prior <- PriorGammaDens(gam)*(sig^(-PriorSigma$a_sig-1))*exp(-PriorSigma$b_sig/sig)
        q.lik = (sig^(-1.5*TT))*exp(-sum(uts)/sig)*
          exp( -0.5*sum( inv.uts*(y^2-2*y*exps+exps2)/sig
                         + (exps-y)*2*(inv.uts*c*abs(gam)*sts + a/sig)
                         + sig*inv.uts*(c^2)*(abs(gam)^2)*sts2
                         + 2*c*abs(gam)*sts*a
                         + (uts*a^2)/sig )/b )
        return(q.prior*q.lik)
      }else{
        log.q.prior <- log(PriorGammaDens(gam)) -(PriorSigma$a_sig+1)*log(sig)-PriorSigma$b_sig/sig
        log.q.lik = - (1.5*TT)*log(sig)-sum(uts)/sig -
          0.5*sum( inv.uts*(y^2-2*y*exps+exps2)/sig
                   + (exps-y)*2*(inv.uts*c*abs(gam)*sts + a/sig)
                   + sig*inv.uts*(c^2)*(abs(gam)^2)*sts2
                   + 2*c*abs(gam)*sts*a
                   + (uts*a^2)/sig )/b
        return(log.q.prior+log.q.lik)
      }
    }
    # importance sampling
    # sample sigma
    if(!fix.sigma){
      sigma.samples = rr_sigma(n.IS)
    }else{
      sigma.samples = rep(sig.init,n.IS)
    }
    # sample gamma
    if(!fix.gamma){
      gamma.samples = rr_gamma(n.IS)
    }else{
      gamma.samples = rep(gam.init,n.IS)
    }
    # compute weights
    if(fix.gamma && fix.sigma){
      weights = rep(1/n.IS,n.IS)
    }else{
      log.weights = apply(cbind(sigma.samples,gamma.samples),1,function(x){dq(x[1],x[2],log.ind=TRUE)}) -
        as.numeric(!fix.gamma)*dr_gamma(gamma.samples,log.ind=TRUE) - as.numeric(!fix.sigma)*dr_sigma(sigma.samples,log.ind = TRUE)
      max.log.weight = max(log.weights)
      log.sum.weights = max(log.weights) + log(sum(exp(log.weights-max.log.weight)))
      log.rescaled.weights = log.weights - log.sum.weights
      weights = exp(log.rescaled.weights)
    }
    # compute expectations
    E.gam = sum(gamma.samples*weights)
    V.gam = sum((gamma.samples^2)*weights) - E.gam^2
    E.sigma = sum(sigma.samples*weights)
    V.sigma = sum((sigma.samples^2)*weights) - E.sigma^2
    E.inv.sigma = sum((1/sigma.samples)*weights)
    E.c2.invb.absgam2.sigma = sum(((sigma.samples*(C.fn(p0,gamma.samples)^2)*(abs(gamma.samples)^2))/B.fn(p0,gamma.samples))*weights)
    E.c.invb.absgam = sum(((C.fn(p0,gamma.samples)*abs(gamma.samples))/B.fn(p0,gamma.samples))*weights)
    E.c.a.invb.absgam = sum(((C.fn(p0,gamma.samples)*A.fn(p0,gamma.samples)*abs(gamma.samples))/B.fn(p0,gamma.samples))*weights)
    E.a2.invb.inv.sigma = sum(((A.fn(p0,gamma.samples)^2)/(B.fn(p0,gamma.samples)*sigma.samples))*weights)
    E.invb.inv.sigma = sum((1/(B.fn(p0,gamma.samples)*sigma.samples))*weights)
    E.a.invb.inv.sigma = sum((A.fn(p0,gamma.samples)/(B.fn(p0,gamma.samples)*sigma.samples))*weights)
    E.log.inv.sigma = sum(log(1/sigma.samples)*weights)

    return(list(E.sigma=E.sigma,V.sigma=V.sigma,E.inv.sigma=E.inv.sigma,E.gam=E.gam,V.gam=V.gam,
                sigma.samples=sigma.samples,gamma.samples=gamma.samples,weights=weights,
                E.c2.invb.absgam2.sigma = E.c2.invb.absgam2.sigma, E.c.invb.absgam = E.c.invb.absgam,
                E.c.a.invb.absgam = E.c.a.invb.absgam, E.a2.invb.inv.sigma = E.a2.invb.inv.sigma,
                E.invb.inv.sigma = E.invb.inv.sigma, E.a.invb.inv.sigma = E.a.invb.inv.sigma,
                E.log.inv.sigma=E.log.inv.sigma))
  }

  tictoc::tic("run time")
  ### estimate posterior
  while( new.max > tol || conv.count < 5 || iter < 15){

    # counter
    iter = iter + 1
    if(verbose & iter%%5==0){
      cat(sprintf("ISVB iteration %s: %s", iter, Sys.time() ),"\n")
    }

    # update distributions
    cur.uts.out = new.uts.out
    cur.sts.out = new.sts.out
    cur.theta.out = new.theta.out
    cur.gamsig.out = new.gamsig.out

    # update q(st)
    new.sts.out <- update_sts(cur.theta.out$exps,cur.uts.out$E.inv.uts,
                              cur.gamsig.out$E.c2.invb.absgam2.sigma,cur.gamsig.out$E.c.invb.absgam,cur.gamsig.out$E.c.a.invb.absgam)

    # update q(ut)
    new.uts.out <- update_uts(cur.theta.out$exps,cur.theta.out$exps2,
                              new.sts.out$E.sts,new.sts.out$E.sts2,
                              cur.gamsig.out$E.inv.sigma,cur.gamsig.out$E.a2.invb.inv.sigma,cur.gamsig.out$E.invb.inv.sigma,
                              cur.gamsig.out$E.c.invb.absgam,cur.gamsig.out$E.c2.invb.absgam2.sigma)

    # update q(theta)
    new.theta.out <- update_theta(cur.gamsig.out$E.c.invb.absgam*new.sts.out$E.sts/cur.gamsig.out$E.invb.inv.sigma +
                                    cur.gamsig.out$E.a.invb.inv.sigma/(new.uts.out$E.inv.uts*cur.gamsig.out$E.invb.inv.sigma),
                                  (cur.gamsig.out$E.invb.inv.sigma*new.uts.out$E.inv.uts)^(-1) )

    # update q(gamma,sigma)
    new.gamsig.out<-update_gamma_sigma(cur.gamsig.out$E.gam,cur.gamsig.out$V.gam,
                                       cur.gamsig.out$E.sigma,cur.gamsig.out$V.sigma,
                                       new.theta.out$exps,new.theta.out$exps2,
                                       new.sts.out$E.sts,new.sts.out$E.sts2,
                                       new.uts.out$E.uts,new.uts.out$E.inv.uts)

    # save ISVB gamma and sigma estimates
    seq.gamma = c(seq.gamma,new.gamsig.out$E.gam)
    seq.sigma = c(seq.sigma,new.gamsig.out$E.sigma)

    # evaluate convergence
    new.max = max(abs(c(cur.theta.out$exps-new.theta.out$exps)))
    conv.count = ifelse(new.max < tol, conv.count + 1, 0)

  }
  run.time = tictoc::toc(quiet = TRUE)
  if(verbose){
    cat(sprintf("ISVB converged: %s iterations, %s seconds",iter,round(run.time$toc-run.time$tic,3)),"\n")
  }

  ### posterior samples
  # gamma and sigma
  samp.index = sample(1:n.IS,n.samp,replace=TRUE,prob=new.gamsig.out$weights)
  samp.gamma = new.gamsig.out$gamma.samples[samp.index]
  samp.sigma = new.gamsig.out$sigma.samples[samp.index]
  # uts, sts, thetas, and predicive distribution samples
  samp.uts = t(sapply(1:TT,function(t){GeneralizedHyperbolic::rgig(n.samp,chi=new.uts.out$uts.chi[t],psi=new.uts.out$uts.psi,lambda=new.uts.out$uts.lambda)}))
  samp.sts = t(sapply(1:TT,function(t){truncnorm::rtruncnorm(n.samp,a=rep(0,n.samp),b=rep(Inf,n.samp),mean=new.sts.out$sts.mu[t],sd=sqrt(new.sts.out$sts.sig2[t]))}))
  samp_theta_t = function(t){
    svd.sC = svd(new.theta.out$sC[,,t]); LL = svd.sC$u%*%diag(sqrt(svd.sC$d),p)
    new.theta.out$sm[,t] + LL%*%matrix(stats::rnorm(n.samp*p,0,1),p,n.samp)}
  samp_post_pred_t = function(t){
        rexal(n.samp,p.fn(p0,samp.gamma),colSums(matrix(FF[,t],p,n.samp)*samp.theta[,t,]) +
           samp.sigma*C.fn(p0,samp.gamma)*abs(samp.gamma)*samp.sts[t,],samp.sigma,0)
   }
  samp.theta = array(NA,c(p,TT,n.samp))
  samp.post.pred = matrix(NA,TT,n.samp)
  for(t in 1:TT){samp.theta[,t,] = samp_theta_t(t); samp.post.pred[t,] = samp_post_pred_t(t)}

  ### list results
  if(!dqlm.ind){
    retlist = list(y=y,run.time=(run.time$toc-run.time$tic),iter=iter,dqlm.ind=dqlm.ind,
                   model=model,p0=p0,df=df,dim.df=dim.df,
                   sig.init=sig.init,seq.sigma=seq.sigma,gam.init=gam.init,seq.gamma=seq.gamma,
                   samp.theta=samp.theta,samp.post.pred=samp.post.pred,
                   map.standard.forecast.errors=new.theta.out$standard.forecast.errors,
                   samp.sigma=samp.sigma,samp.gamma=samp.gamma,samp.sts=samp.sts,samp.vts=samp.uts,
                   theta.out=new.theta.out,gammasig.out=new.gamsig.out,sts.out=new.sts.out,vts.out=new.uts.out,
                   fix.sigma=fix.sigma,fix.gamma=fix.gamma)
  }else{
    retlist = list(y=y,run.time=(run.time$toc-run.time$tic),iter=iter,dqlm.ind=dqlm.ind,
                   model=model,p0=p0,df=df,dim.df=dim.df,
                   sig.init=sig.init,seq.sigma=seq.sigma,
                   samp.theta=samp.theta,samp.post.pred=samp.post.pred,
                   map.standard.forecast.errors=new.theta.out$standard.forecast.errors,
                   samp.sigma=samp.sigma,samp.vts=samp.uts,
                   theta.out=new.theta.out,sig.out=new.gamsig.out,vts.out=new.uts.out,
                   fix.sigma=fix.sigma)
  }
  # return results
  class(retlist) <- "exdqlmISVB"
  return(retlist)
}
