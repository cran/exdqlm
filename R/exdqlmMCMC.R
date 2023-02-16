#' exDQLM - MCMC algorithm
#'
#' The function applies a Markov chain Monte Carlo (MCMC) algorithm to sample the posterior of an exDQLM.
#'
#' @inheritParams exdqlmISVB
#' @param Sig.mh Covariance matrix used in the random walk MH step to jointly sample sigma and gamma.
#' @param joint.sample Logical value indicating whether or not to recompute `Sig.mh` based off the initial burn-in samples of gamma and sigma. Default is `FALSE`.
#' @param n.burn Number of MCMC iterations to burn. Default is `n.burn = 2000`.
#' @param n.mcmc Number of MCMC iterations to sample. Default is `n.mcmc = 1500`.
#' @param init.from.isvb Logical value indicating whether or not to initialize the MCMC using the ISVB algorithm. Default is `TRUE`.
#'
#' @return A list of the following is returned:
#'  \itemize{
#'   \item `run.time` - Algorithm run time in seconds.
#'   \item `model` - List of the state-space model including `GG`, `FF`, prior parameters `m0` and `C0`.
#'   \item `p0` - The quantile which was estimated.
#'   \item `df` - Discount factors used for each block.
#'   \item `dim.df` - Dimension used for each block of discount factors.
#'   \item `samp.theta` - Posterior sample of the state vector.
#'   \item `samp.post.pred` - Sample of the posterior predictive distributions.
#'   \item `map.standard.forecast.errors` - MAP standardized one-step-ahead forecast errors.
#'   \item `samp.sigma` - Posterior sample of scale parameter sigma.
#'   \item `samp.vts` - Posterior sample of latent parameters, v_t.
#'   \item `theta.out` - List containing the distributions of the state vector including filtered distribution parameters (`fm` and `fC`) and smoothed distribution parameters (`sm` and `sC`).
#' }
#' If `dqlm.ind=FALSE`, the list also contains the following:
#' \itemize{
#'   \item `samp.gamma` - Posterior sample of skewness parameter gamma.
#'   \item `samp.sts` - Posterior sample of latent parameters, s_t.
#'   \item `init.log.sigma` - Burned samples of log sigma from the random walk MH joint sampling of sigma and gamma.
#'   \item `init.logit.gamma` - Burned samples of logit gamma from the random walk MH joint sampling of sigma and gamma.
#'   \item `accept.rate` - Acceptance rate of the MH step.
#'   \item `Sig.mh` - Covariance matrix used in MH step to jointly sample sigma and gamma.
#' }
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' trend.comp = polytrendMod(1,mean(y),10)
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = combineMods(trend.comp,seas.comp)
#' M2 = exdqlmMCMC(y,p0=0.85,model,df=c(1,1),dim.df = c(1,6),
#'                 gam.init=-3.5,sig.init=15,
#'                 n.burn=100,n.mcmc=150)
#' }
#'
exdqlmMCMC <- function(y,p0,model,df,dim.df,fix.gamma=FALSE,gam.init=NA,fix.sigma=FALSE,sig.init=NA,dqlm.ind=FALSE,
                    Sig.mh,joint.sample=FALSE,n.burn=2000,n.mcmc=1500,init.from.isvb=TRUE,PriorSigma=NULL,PriorGamma=NULL,verbose=TRUE){

  # check inputs
  y = check_ts(y)
  model = check_mod(model)
  rv = check_logics(gam.init,sig.init,fix.gamma,fix.sigma,dqlm.ind)
  gam.init = rv$gam.init
  dqlm.int = rv$dqlm.ind
  fix.gamma = rv$fix.gamma

  ### MCMC iterations
  if(n.mcmc<=0){
    stop("number of mcmc samples must be positive")
    }
  if(verbose & n.burn<=0){
    warning("mcmc will be sampled without burn-in, a burn-in is recommended even if initializing using the isvb algorithm")
    n.burn=0
    }
  I = n.mcmc + n.burn

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
  PriorSigmaDens<-function(sigma){ LaplacesDemon::dinvgamma(sigma,shape=PriorSigma$a_sig,scale=PriorSigma$b_sig)  }
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

  # function to produce smoothed estimates for return value
  smoothed_theta<-function(ex.f,ex.q){
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
    return(list(standard.forecast.errors=standard.forecast.errors,sm=sm,sC=sC,fm=m,fC=C))
  }

  ### Initialize MCMC
  init.log.sigma <- init.logit.gamma <- rep(NA,n.burn)
  save.sigma <- save.gamma <- rep(NA,n.mcmc)
  save.Ut <- save.st <- matrix(NA,TT,n.mcmc)
  save.theta <- array(NA,c(p,TT,n.mcmc))
  save.post.pred <- matrix(NA,TT,n.mcmc)
  # Set initial values
  if(init.from.isvb){
    if(verbose){
      cat("running isvb algorithm to initialize mcmc","\n")
    }
    isvb.out <- exdqlmISVB(y,p0,model,df,dim.df,fix.gamma,gam.init,fix.sigma,sig.init,dqlm.ind,tol=0.5,n.IS=200,PriorSigma=PriorSigma,PriorGamma=PriorGamma,verbose=verbose)
    cursam.sigma <- ifelse(fix.sigma,sig.init,ifelse(dqlm.ind,isvb.out$sig.out$E.sigma,isvb.out$gammasig.out$E.sigma))
    cursam.Ut <- isvb.out$vts.out$E.uts
    cursam.theta <- isvb.out$theta.out$sm
  }else{
    cursam.sigma <- m_sigma
    cursam.Ut <- rep(1/m_sigma,TT)
    cursam.theta <- matrix(m0,p,TT)
  }

  ######## exDQLM
  if(!dqlm.ind){

    ### Define logit and inverse logit functions
    logit = function(x){log((x-L)/(U-x))}
    inv.logit = function(x){(U*exp(x)+L)/(exp(x)+1)}

    ### Additional initial values
    if(init.from.isvb){
      cursam.st <- isvb.out$sts.out$E.sts
      cursam.gamma <- ifelse(fix.gamma,gam.init,isvb.out$gammasig.out$E.gam)
      cursam.logit.gamma <- logit(cursam.gamma)
      cursam.log.sigma <- log(cursam.sigma)
    }else{
      cursam.st <- truncnorm::rtruncnorm(TT,a=0,b=Inf,mean=0,sd=1)
      cursam.gamma <- ifelse(!is.na(gam.init),gam.init,(L+U)/2)
      cursam.logit.gamma <- logit(cursam.gamma)
      cursam.log.sigma <- log(cursam.sigma)
    }

    ### Initialize MH
    n.accept = 0
    if(!methods::hasArg(Sig.mh)){
      if(init.from.isvb){
        Sig.mh <- stats::cov(cbind(log(isvb.out$gammasig.out$sigma.samples),logit(isvb.out$gammasig.out$gamma.samples)))
      }else{
        Sig.mh = diag(c(ifelse(fix.sigma,0,0.005),ifelse(fix.gamma,0,0.005)))
      }
    }else{
      new.Sig.mh = Sig.mh
      if(fix.sigma){
        new.Sig.mh = matrix(0,2,2)
        new.Sig.mh[2,2] = Sig.mh[2,2]
      }
      if(fix.gamma){
        new.Sig.mh = matrix(0,2,2)
        new.Sig.mh[1,1] = Sig.mh[1,1]
      }
      Sig.mh = new.Sig.mh
    }
    if(fix.gamma | fix.sigma){
      chol_Sig.mh = sqrt(Sig.mh)
    }else{
      chol_Sig.mh=t(chol(Sig.mh))
    }

    # exdqlm function sample theta ffbs
    ex_samp_theta<-function(ex.f,ex.q,gamma,sigma,sts,tau,c_tau){
      # initialize ffbs
      m <- sam.theta <- matrix(NA,p,TT)
      C <- array(NA,c(p,p,TT))
      standard.forecast.errors <- post.pred <- rep(NA,TT)
      ## forward filter
      # first iteration
      a = as.vector(GG[,,1]%*%m0)
      P = GG[,,1]%*%C0%*%t(GG[,,1])
      R = P + df.mat*P
      R = (R + t(R))/2
      f = t(FF[,1])%*%a + ex.f[1]
      q = t(FF[,1])%*%R%*%FF[,1] + ex.q[1]
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
      ## backwards sample
      svd.sC = svd(C[,,TT])
      sam.theta[,TT] = m[,TT] + svd.sC$u%*%diag(sqrt(svd.sC$d),p)%*%stats::rnorm(p,0,1)
      post.pred[TT] = brms::rasym_laplace(1,t(FF[,TT])%*%sam.theta[,TT]+c_tau*sigma*abs(gamma)*sts[TT],sigma,tau)
      for(t in (TT-1):1){
        P = GG[,,(t+1)]%*%C[,,(t)]%*%t(GG[,,(t+1)])
        R = P + df.mat*P
        R = (R + t(R))/2
        svd.R = svd(R)
        inv.R = svd.R$u%*%diag(1/svd.R$d,p)%*%t(svd.R$u)
        sB = C[,,t]%*%t(GG[,,t])%*%inv.R
        sm = m[,t] + sB%*%(sam.theta[,(t+1)]-as.vector(GG[,,(t+1)]%*%m[,(t)]))
        sC = C[,,t] - sB%*%GG[,,t]%*%C[,,t]
        svd.sC = svd((sC+t(sC))/2)
        sam.theta[,t] = sm + svd.sC$u%*%diag(sqrt(svd.sC$d),p)%*%stats::rnorm(p,0,1)
        post.pred[t] = brms::rasym_laplace(1,t(FF[,t])%*%sam.theta[,t]+c_tau*sigma*abs(gamma)*sts[t],sigma,tau)
      }
      return(list(standard.forecast.errors=standard.forecast.errors,post.pred=post.pred,sam.theta=sam.theta,fm=m,fC=C))
    }

    # exdqlm function sample uts
    ex_samp_uts<-function(reg1,gamma,sigma,sts,a_tau,b_tau,c_tau){
      apply(((y-reg1-sigma*c_tau*abs(gamma)*sts)^2)/(b_tau*sigma),1,
            function(x){GeneralizedHyperbolic::rgig(1,chi=x,psi = (a_tau^2)/(b_tau*sigma) + (2/sigma), lambda = 0.5)})
    }

    # exdqlm function sample sts
    ex_samp_sts<-function(reg1,gamma,sigma,uts,a_tau,b_tau,c_tau){
      s.sig2<-1/(1+c_tau^2*abs(gamma)^2*sigma/(b_tau*uts))
      s.mu<-s.sig2*c_tau*abs(gamma)*(y-(reg1+a_tau*uts))/(b_tau*uts)
      truncnorm::rtruncnorm(TT,rep(0,TT),rep(Inf,TT),s.mu,sqrt(s.sig2))
    }

    # exdqlm function sample sigma and gamma - MH
    logL<-function(reg1,log.sigma,logit.gamma,sts,uts){
      sigma=exp(log.sigma); gamma=inv.logit(logit.gamma)
      temp.p = p.fn(p0,gamma)
      a = (1-2*temp.p)/(temp.p*(1-temp.p))
      b = (2)/(temp.p*(1-temp.p))
      c = (as.numeric(gamma>0)-temp.p)^(-1)
      logJ<-logit.gamma-2*log(1+exp(logit.gamma))+log.sigma
      PriorGamma<- PriorGammaDens(gamma)
      PriorSigma<- PriorSigmaDens(sigma)
      sum(stats::dnorm(y,reg1+sigma*c*abs(gamma)*sts+a*uts,sqrt(sigma*b*uts),log = TRUE)) +
        sum(stats::dexp(uts,rate = 1/sigma,log=TRUE)) +
        log(PriorSigma) + log(PriorGamma) + logJ
    }
    ex_samp_lsiglgam<-function(reg1,log.sigma,logit.gamma,sts,uts,chol_Sig){
      prop<-c(log.sigma,logit.gamma)+chol_Sig%*%stats::rnorm(2)
      if(inv.logit(prop[2]) < U && inv.logit(prop[2]) > L){
        logr<-logL(reg1,prop[1],prop[2],sts,uts)-logL(reg1,log.sigma,logit.gamma,sts,uts)
        accept=(log(stats::runif(1))<logr)
      }else{
        accept=FALSE
      }
      log.sigma.new<-accept*prop[1]+(1-accept)*log.sigma
      logit.gamma.new<-accept*prop[2]+(1-accept)*logit.gamma
      return(list(log.sigma=log.sigma.new,logit.gamma=logit.gamma.new,accept=accept))
    }

    # Sample from exdqlm posterior
    tictoc::tic()
    for (i in 1:I){
      # counter
      if(verbose & i%%500==0){
        cat(sprintf("%s iteration %s, acceptance rate %s: %s", ifelse(i<=n.burn,"burn-in","MCMC"), i , round(n.accept/i,4), Sys.time()),"\n")
        }

      # exAL parameters
      tau = p.fn(p0,cursam.gamma)
      a_tau = (1-2*tau)/(tau*(1-tau))
      b_tau = (2)/(tau*(1-tau))
      c_tau = (as.numeric(cursam.gamma>0)-tau)^(-1)

      # sample theta
      ex.f = cursam.sigma*c_tau*abs(cursam.gamma)*cursam.st + cursam.Ut*a_tau
      ex.q = b_tau*cursam.Ut*cursam.sigma
      theta.out <- ex_samp_theta(ex.f,ex.q,cursam.gamma,cursam.sigma,cursam.st,tau,c_tau)
      cursam.theta = theta.out$sam.theta

      # sample uts, sts
      reg1 = apply(FF*cursam.theta,2,sum)
      cursam.Ut<-ex_samp_uts(reg1,cursam.gamma,cursam.sigma,cursam.st,a_tau,b_tau,c_tau)
      cursam.st<-ex_samp_sts(reg1,cursam.gamma,cursam.sigma,cursam.Ut,a_tau,b_tau,c_tau)

      # sample sigma and gamma - MH
      lsiglgam.out<-ex_samp_lsiglgam(reg1,cursam.log.sigma,cursam.logit.gamma,cursam.st,cursam.Ut,chol_Sig.mh)
      cursam.gamma<-inv.logit(lsiglgam.out$logit.gamma)
      cursam.logit.gamma<-lsiglgam.out$logit.gamma
      cursam.sigma<-exp(lsiglgam.out$log.sigma)
      cursam.log.sigma<-lsiglgam.out$log.sigma
      n.accept = n.accept + lsiglgam.out$accept

      # save samples after burn
      if(i <= n.burn){
        init.log.sigma[i] = cursam.log.sigma
        init.logit.gamma[i] = cursam.logit.gamma
        if(i==n.burn && joint.sample){
          Sig.mh = stats::cov(cbind(init.log.sigma[1:n.burn],init.logit.gamma[1:n.burn]))
          if(fix.gamma | fix.sigma){
            chol_Sig.mh = sqrt(Sig.mh)
          }else{
            chol_Sig.mh=t(chol(Sig.mh))
          }
          }
      }else{
        save.sigma[(i-n.burn)] = cursam.sigma
        save.gamma[(i-n.burn)] = cursam.gamma
        save.theta[,,(i-n.burn)] = cursam.theta
        save.Ut[,(i-n.burn)] = cursam.Ut
        save.st[,(i-n.burn)] = cursam.st
        save.post.pred[,(i-n.burn)] = theta.out$post.pred
      }

    }
    run.time = tictoc::toc(quiet = TRUE)
    if(verbose){
      cat(sprintf("MCMC complete: %s iterations, %s seconds",I,round(run.time$toc-run.time$tic,3)),"\n")
    }

    # exdqlm MAP standard forecast errors
    map.gam = mean(save.gamma)
    map.sig = mean(save.sigma)
    map.st = rowMeans(save.st)
    map.Ut = rowMeans(save.Ut)
    tau = p.fn(p0,map.gam)
    a_tau = (1-2*tau)/(tau*(1-tau))
    b_tau = (2)/(tau*(1-tau))
    c_tau = (as.numeric(map.gam>0)-tau)^(-1)
    theta.out <- smoothed_theta(map.sig*c_tau*abs(map.gam)*map.st+map.Ut*a_tau,b_tau*map.Ut*map.sig)
    map.standard.forecast.errors = theta.out$standard.forecast.errors

    # exdqlm results
    retlist = list(run.time=(run.time$toc-run.time$tic),model=model,p0=p0,df=df,dim.df=dim.df,
                samp.theta = coda::as.mcmc(save.theta), theta.out = theta.out,
                samp.post.pred = save.post.pred, map.standard.forecast.errors = map.standard.forecast.errors,
                samp.sigma = coda::as.mcmc(save.sigma), samp.gamma = coda::as.mcmc(save.gamma),
                init.log.sigma = coda::as.mcmc(init.log.sigma), init.logit.gamma = coda::as.mcmc(init.logit.gamma),
                samp.vts = coda::as.mcmc(save.Ut), samp.sts = coda::as.mcmc(save.st),
                accept.rate = n.accept/I, Sig.mh=Sig.mh)

  }else{
    ######## DQLM

    # fixed AL parameters
    a_tau = (1-2*p0)/(p0*(1-p0))
    b_tau = (2)/(p0*(1-p0))

    # dqlm function sample theta ffbs
    samp_theta<-function(ex.f,ex.q,sigma){
      # initialize ffbs
      m <- sam.theta <- matrix(NA,p,TT)
      C <- array(NA,c(p,p,TT))
      standard.forecast.errors <- post.pred <- rep(NA,TT)
      ## forward filter
      # first iteration
      a = as.vector(GG[,,1]%*%m0)
      P = GG[,,1]%*%C0%*%t(GG[,,1])
      R = P + df.mat*P
      R = (R + t(R))/2
      f = t(FF[,1])%*%a + ex.f[1]
      q = t(FF[,1])%*%R%*%FF[,1] + ex.q[1]
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
      ## backwards sample
      svd.sC = svd(C[,,TT])
      sam.theta[,TT] = m[,TT] + svd.sC$u%*%diag(sqrt(svd.sC$d),p)%*%stats::rnorm(p,0,1)
      post.pred[TT] = brms::rasym_laplace(1,t(FF[,TT])%*%sam.theta[,TT],sigma,p0)
      for(t in (TT-1):1){
        P = GG[,,(t+1)]%*%C[,,(t)]%*%t(GG[,,(t+1)])
        R = P + df.mat*P
        R = (R + t(R))/2
        svd.R = svd(R)
        inv.R = svd.R$u%*%diag(1/svd.R$d,p)%*%t(svd.R$u)
        sB = C[,,t]%*%t(GG[,,t])%*%inv.R
        sm = m[,t] + sB%*%(sam.theta[,(t+1)]-as.vector(GG[,,(t+1)]%*%m[,(t)]))
        sC = C[,,t] - sB%*%GG[,,t]%*%C[,,t]
        svd.sC = svd((sC+t(sC))/2)
        sam.theta[,t] = sm + svd.sC$u%*%diag(sqrt(svd.sC$d),p)%*%stats::rnorm(p,0,1)
        post.pred[t] = brms::rasym_laplace(1,t(FF[,t])%*%sam.theta[,t],sigma,p0)
      }
      return(list(standard.forecast.errors=standard.forecast.errors,post.pred=post.pred,sam.theta=sam.theta,fm=m,fC=C))
    }

    # dqlm function sample uts
    samp_uts<-function(reg1,sigma){
      apply(((y-reg1)^2)/(b_tau*sigma),1,
            function(x){GeneralizedHyperbolic::rgig(1,chi=x,psi = (a_tau^2)/(b_tau*sigma) + (2/sigma), lambda = 0.5)})
    }

    # dqlm function sample sigma
    samp_sigma<-function(reg1,uts){
      1/stats::rgamma(1, shape = PriorSigma$a_sig + 1.5*TT,
               rate = PriorSigma$b_sig + 0.5*sum( ((as.vector(y) - reg1 - a_tau*uts)^2)/(b_tau*uts) ) + sum(uts) )
    }

    # Sample from dqlm posterior
    tictoc::tic()
    for (i in 1:I){
      # counter
      if(verbose & i%%500==0){
        cat(sprintf("%s iteration %s: %s ", ifelse(i<=n.burn,"burn-in","MCMC"), i, Sys.time()), "\n")
        }

      # sample theta
      ex.f = cursam.Ut*a_tau
      ex.q = b_tau*cursam.Ut*cursam.sigma
      theta.out <- samp_theta(ex.f,ex.q,cursam.sigma)
      cursam.theta = theta.out$sam.theta

      # sample uts
      reg1 = apply(FF*cursam.theta,2,sum)
      cursam.Ut<-samp_uts(reg1,cursam.sigma)

      # sample sigma
      if(!fix.sigma){
        cursam.sigma <- samp_sigma(reg1,cursam.Ut)
      }

      # save samples after burn
      if(i > n.burn){
        save.sigma[(i-n.burn)] = cursam.sigma
        save.theta[,,(i-n.burn)] = cursam.theta
        save.Ut[,(i-n.burn)] = cursam.Ut
        save.post.pred[,(i-n.burn)] = theta.out$post.pred
      }

    }
    run.time = tictoc::toc(quiet = TRUE)
    if(verbose){
      cat(sprintf("MCMC complete: %s iterations, %s seconds",I,round(run.time$toc-run.time$tic,3)),"\n")
    }

    # dqlm MAP standard forecast errors
    map.sig = mean(save.sigma)
    map.Ut = rowMeans(save.Ut)
    theta.out <- smoothed_theta(map.Ut*a_tau,b_tau*map.Ut*map.sig)
    map.standard.forecast.errors = theta.out$standard.forecast.errors

    # dqlm results
    retlist = list(run.time=(run.time$toc-run.time$tic),model=model,p0=p0,df=df,dim.df=dim.df,
                samp.theta = coda::as.mcmc(save.theta), theta.out = theta.out,
                samp.post.pred = save.post.pred, map.standard.forecast.errors = map.standard.forecast.errors,
                samp.sigma = coda::as.mcmc(save.sigma),
                samp.vts = coda::as.mcmc(save.Ut))
  }

  # return results
  class(retlist) <- "exdqlm"
  return(retlist)
}
