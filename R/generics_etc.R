##################################
######## "exdqlm" objects ########
##################################

#' \code{exdqlm} objects
#'
#' \code{is.exdqlm} tests if its argument is a \code{exdqlm} object. 
#' 
#' @usage is.exdqlm(m)
#'
#' @param m an \strong{R} object
#'
#' @export
is.exdqlm = function(m){ return(methods::is(m,"exdqlm")) }

#' \code{exdqlm} objects
#'
#' \code{as.exdqlm} attempts to turn a list into an \code{exdqlm} object. Works for time-invariant \code{dlm} objects created using the \pkg{dlm} package. 
#' 
#' @usage as.exdqlm(m)
#'
#' @param m a list containing named elements m0, C0, FF and GG.
#'
#' @return A object of class "\code{exdqlm}" containing the state space model components:
#' \itemize{
#'   \item FF - Observational vector.
#'   \item GG - Evolution matrix.
#'   \item m0 - Prior mean of the state vector.
#'   \item C0 - Prior covariance of the state vector.
#' }
#' @export
as.exdqlm <- function(m){
  if(is.exdqlm(m)){
    return(m)
  }
  if(!is.list(m)){
    stop("Input must be a list with named elements m0, C0, FF and GG.")
  }
  if(methods::is(m,"dlm")){
    if(!is.null(m$JFF) | !is.null(m$JGG) |
       !is.null(m$JV) | !is.null(m$JW)){
      stop("'dlm' object input must be a time-invariant")
    }
    l$FF = t(m$FF)
  }
  
  # check for required components & remove extras
  refnn <- c("m0","C0","FF","GG")
  nn <- names(m)
  check <- !sapply(m, is.null)
  ind <- match(refnn,nn)
  if(anyNA(ind)){
    stop(paste("Component(s)",paste(refnn[is.na(ind)], collapse = ", "), "is (are) missing."))
  }
  final.ind = match(nn[ind][check[ind]],nn)
  model = m[final.ind]
  
  class(model) <- "exdqlm"
  model = check_mod(model)
  
  return(model)
}

#' Addition for \code{exdqlm} objects
#'
#' Combines two state space blocks into a single state space model for an exDQLM.
#' 
#' @method + exdqlm
#' @rdname plus-exdqlm
#'
#' @param m1 object of class "\code{exdqlm}" containing the first model to be combined.
#' @param m2 object of class "\code{exdqlm}" containing the second model to be combined.
#'
#' @return A object of class "\code{exdqlm}" containing the new combined state space model components:
#' \itemize{
#'   \item FF - Observational vector.
#'   \item GG - Evolution matrix.
#'   \item m0 - Prior mean of the state vector.
#'   \item C0 - Prior covariance of the state vector.
#' }
#'
#' @examples
#' trend.comp = polytrendMod(2,rep(0,2),10*diag(2))
#' seas.comp = seasMod(365,c(1,2,4),C0=10*diag(6))
#' model = trend.comp + seas.comp
#'
#' @export
"+.exdqlm" <- function(m1, m2){
  m1 = check_mod(m1)
  m2 = check_mod(m2)
  n = length(m1$m0) + length(m2$m0)
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
  model$m0 = matrix(c(m1$m0,m2$m0),n,1)
  model$C0 = magic::adiag(m1$C0,m2$C0)
  
  class(model) <- "exdqlm"
  return(model)
}

#' Print exDQLM model details
#'
#' Print the details of the exDQLM model.
#' @param x a \code{exdqlm} object.
#' @param ... further arguments (unused).
#' 
#' @export
print.exdqlm <- function(x,...){
  refnn <- c("m0","C0","FF","GG")
  descrip = c("Prior mean of the state vector:", 
              "Prior covariance of the state vector:",
              "Observational vector:",
              "Evolution matrix:")
  nn <- names(x)
  check <- !sapply(x, is.null)
  ind <- match(refnn,nn)
  ind <- ind[!is.na(ind)]
  final.ind = match(nn[ind][check[ind]],nn)
  # print
  for (i in 1:4){
    cat(descrip[i],"\n")
    print(x[final.ind[i]])
    cat("\n")
  }
}

#' Summary exDQLM model details
#'
#' Print the details of the exDQLM model.
#' @param object a \code{exdqlm} object.
#' @param ... further arguments (unused).
#' 
#' @export
summary.exdqlm <- function(object,...){
  refnn <- c("m0","C0","FF","GG")
  descrip = c("Prior mean of the state vector:", 
              "Prior covariance of the state vector:",
              "Observational vector:",
              "Evolution matrix:")
  nn <- names(object)
  check <- !sapply(object, is.null)
  ind <- match(refnn,nn)
  ind <- ind[!is.na(ind)]
  final.ind = match(nn[ind][check[ind]],nn)
  # print
  for (i in 1:4){
    cat(descrip[i],"\n")
    print(object[final.ind[i]])
    cat("\n")
  }
}



##################################
###### "exdqlmMCMC" objects ######
##################################

#' \code{exdqlmMCMC} objects
#'
#' \code{is.exdqlmMCMC} tests if its argument is a \code{exdqlmMCMC} object. 
#' 
#' @usage is.exdqlmMCMC(m)
#'
#' @param m an \strong{R} object
#'
#' @export
is.exdqlmMCMC = function(m){ return(methods::is(m,"exdqlmMCMC")) }


#' Print Method for \code{exdqlmMCMC} Objects
#'
#' @param x An \code{exdqlmMCMC} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M2 = exdqlmMCMC(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                 gam.init=-3.5,sig.init=15,
#'                 n.burn=100,n.mcmc=150)
#' print(M2)                
#' }
#'
print.exdqlmMCMC <- function(x, ...) {
  cat("Bayesian Dynamic Quantile Regression Model (exDQLM)\n")
  cat("Number of Observations:", length(x$y), "\n")
  cat("State Dimension:", length(x$model$m0), "\n")  
  cat("Discount factors ( dimensions ):", paste(x$df,"(", x$dim.df, ")",collapse = ", "),"\n \n")
  #
  cat("exDQLM fitted using MCMC\n")
  cat("Burn-in:", x$n.burn, ", MCMC samples:", x$n.mcmc , "\n")
  cat("Run-time:", x$run.time, "seconds\n")
}

#' Summary Method for \code{exdqlmMCMC} Objects
#'
#' @param object An \code{exdqlmMCMC} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M2 = exdqlmMCMC(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                 gam.init=-3.5,sig.init=15,
#'                 n.burn=100,n.mcmc=150)
#' summary(M2)                
#' }
#'
summary.exdqlmMCMC <- function(object, ...) {
  cat("Bayesian Dynamic Quantile Regression Model (exDQLM)\n")
  cat("Number of Observations:", length(object$y), "\n")
  cat("State Dimension:", length(object$model$m0), "\n")  
  cat("Discount factors ( dimensions ):", paste(object$df,"(", object$dim.df, ")",collapse = ", "),"\n \n")
  #
  cat("exDQLM fitted using MCMC\n")
  cat("Burn-in:", object$n.burn, ", MCMC samples:", object$n.mcmc , "\n")
  cat("Run-time:", object$run.time, "seconds\n")
}

#' Plot Method for \code{exdqlmMCMC} Objects
#'
#' @param x An \code{exdqlmMCMC} object.
#' @param ... Additional arguments.
#' 
#' @export
#'
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M2 = exdqlmMCMC(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                 gam.init=-3.5,sig.init=15,
#'                 n.burn=100,n.mcmc=150)
#' plot(M2)                
#' }
#'
plot.exdqlmMCMC<- function(x, ...) {
  exdqlmPlot(x,...)
}



##################################
###### "exdqlmISVB" objects ######
##################################

#' \code{exdqlmISVB} objects
#'
#' \code{is.exdqlmISVB} tests if its argument is a \code{exdqlmISVB} object. 
#' 
#' @usage is.exdqlmISVB(m)
#'
#' @param m an \strong{R} object
#'
#' @export
#' 
#' 
is.exdqlmISVB = function(m){ return(methods::is(m,"exdqlmISVB")) }

#' Print Method for \code{exdqlmISVB} Objects
#'
#' @param x An \code{exdqlmISVB} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                    gam.init=-3.5,sig.init=15)
#' print(M0)
#' }
#'
print.exdqlmISVB <- function(x, ...) {
  cat("Bayesian Dynamic Quantile Regression Model (exDQLM)\n")
  cat("Number of Observations:", length(x$y), "\n")
  cat("State Dimension:", length(x$model$m0), "\n")  
  cat("Discount factors ( dimensions ):", paste(x$df,"(", x$dim.df, ")",collapse = ", "),"\n \n")
  #
  cat("exDQLM fitted using ISVB\n")
  cat("Variational Parameters:", paste(if(!x$fix.gamma){"gamma"}, if(!x$fix.sigma){"sigma"}, if(x$fix.sigma && x$fix.gamma){"none"}, collapse=", ") , "\n")
  cat("Iterations until convergence:", x$iter, "\n")
  cat("Run-time:", x$run.time, "seconds\n")
}

#' Summary Method for \code{exdqlmISVB} Objects
#'
#' @param object An \code{exdqlmISVB} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                    gam.init=-3.5,sig.init=15)
#' summary(M0)
#' }
#'
summary.exdqlmISVB <- function(object, ...) {
  cat("Bayesian Dynamic Quantile Regression Model (exDQLM)\n")
  cat("Number of Observations:", length(object$y), "\n")
  cat("State Dimension:", length(object$model$m0), "\n")  
  cat("Discount factors ( dimensions ):", paste(object$df,"(", object$dim.df, ")",collapse = ", "),"\n \n")
  #
  cat("exDQLM fitted using ISVB\n")
  cat("Variational Parameters:", paste(if(!object$fix.gamma){"gamma"}, if(!object$fix.sigma){"sigma"}, if(object$fix.sigma && object$fix.gamma){"none"}, collapse=", ") , "\n")
  cat("Iterations until convergence:", object$iter, "\n")
  cat("Run-time:", object$run.time, "seconds\n")
}

#' Plot Method for \code{exdqlmISVB} Objects
#'
#' @param x An \code{exdqlmISVB} object.
#' @param ... Additional arguments. 
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,quantile(y,0.85),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.98),dim.df = c(1),
#'                    gam.init=-3.5,sig.init=15)
#' plot(M0)
#' }
#'
plot.exdqlmISVB <- function(x, ...) {
  exdqlmPlot(x,...)
}




##################################
### "exdqlmDiagnostic" objects ###
##################################

#' \code{exdqlmDiagnostic} objects
#'
#' \code{is.exdqlmDiagnostic} tests if its argument is a \code{exdqlmDiagnostic} object. 
#' 
#' @usage is.exdqlmDiagnostic(x)
#'
#' @param x an \strong{R} object
#'
#' @export
is.exdqlmDiagnostic = function(x){ return(methods::is(x,"exdqlmDiagnostic")) }

#' Print Method for \code{exdqlmDiagnostic} Objects
#'
#' @param x An \code{exdqlmDiagnostic} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,mean(y),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.95),dim.df = c(1),
#'                   gam.init=-3.5,sig.init=15)
#' M0.diags = exdqlmDiagnostics(M0,plot=FALSE)
#' print(M0.diags)
#' }
#'
print.exdqlmDiagnostic <- function(x, ...) {
  #
  Diagnostic <- c("KL","pplc","run-time (s)")
  M1 <- c(x$m1.KL,x$m1.pplc,as.numeric(x$m1.rt))
  #
  if(is.null(x$m2.KL)){
    print(data.frame(Diagnostic=Diagnostic,M1=M1), row.names = FALSE, digits = 3)
  }else{
    M2 <- c(x$m2.KL,x$m2.pplc,as.numeric(x$m2.rt))
    print(data.frame(Diagnostic=Diagnostic,M1=M1,M2=M2), row.names = FALSE, digits = 3)
  }
}

#' Summary Method for \code{exdqlmDiagnostic} Objects
#'
#' @param object An \code{exdqlmDiagnostic} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,mean(y),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.95),dim.df = c(1),
#'                   gam.init=-3.5,sig.init=15)
#' M0.diags = exdqlmDiagnostics(M0,plot=FALSE)
#' summary(M0.diags)
#' }
#'
summary.exdqlmDiagnostic <- function(object, ...) {
  #
  Diagnostic <- c("KL","pplc","run-time (s)")
  M1 <- c(object$m1.KL,object$m1.pplc,as.numeric(object$m1.rt))
  #
  if(is.null(object$m2.KL)){
    print(data.frame(Diagnostic=Diagnostic,M1=M1), row.names = FALSE, digits = 3)
  }else{
    M2 <- c(object$m2.KL,object$m2.pplc,as.numeric(object$m2.rt))
    print(data.frame(Diagnostic=Diagnostic,M1=M1,M2=M2), row.names = FALSE, digits = 3)
  }
}

#' Plot Method for \code{exdqlmDiagnostic} Objects
#'
#' @param x An \code{exdqlmDiagnostic} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' y = scIVTmag[1:100]
#' model = polytrendMod(1,mean(y),10)
#' M0 = exdqlmISVB(y,p0=0.85,model,df=c(0.95),dim.df = c(1),
#'                   gam.init=-3.5,sig.init=15)
#' M0.diags = exdqlmDiagnostics(M0,plot=FALSE)
#' plot(M0.diags)
#' }
#'
plot.exdqlmDiagnostic <- function(x, ...) {
  aa = list(...)
  if(is.null(aa$cols)){cols=c("red","blue")}else{cols = aa$cols}
  
  # get ranges
  if(is.null(x$m2.KL)){
    qq.x.range = range(x$m1.qq$x)
    qq.y.range = range(x$m1.qq$y)
    acf.y.range = range(x$m1.acf$acf)
    fe.y.range = range(x$m1.msfe)
  }else{
    qq.x.range = range(c(x$m1.qq$x,x$m2.qq$x))
    qq.y.range = range(c(x$m1.qq$y,x$m2.qq$y))
    acf.y.range = range(c(x$m1.acf$acf,x$m2.acf$acf))
    fe.y.range = range(c(x$m1.msfe,x$m2.msfe))
  }
  # m1 qqplot
  plot(x$m1.qq,main="",col=cols[1],pch=20,xlab="Theoretical Quantiles",ylab="M1 Sample Quantiles",xlim=qq.x.range,ylim=qq.y.range)
  graphics::abline(a=0,b=1)
  # m1 acf
  plot(x$m1.acf,ylab="M1 ACF",col=cols[1],main="",ylim=acf.y.range)
  # m1 forecast errors
  ts.xy = grDevices::xy.coords(x$y)
  graphics::plot(ts.xy$x,x$m1.msfe,ylab="M1 standard forecast errors",xlab="time",col=cols[1],pch=20,type="l",ylim=fe.y.range)
  graphics::abline(h=0,lty=2)
  ### m2
  if(!is.null(x$m2.KL)){
    # m2 qqplot
    plot(x$m2.qq,main="",col=cols[2],pch=20,xlab="Theoretical Quantiles",ylab="M2 Sample Quantiles",xlim=qq.x.range,ylim=qq.y.range)
    graphics::abline(a=0,b=1)
    # m2 acf
    plot(x$m2.acf,ylab="M2 ACF",col=cols[2],main="",ylim=acf.y.range)
    # m2 forecast errors
    graphics::plot(ts.xy$x,x$m2.msfe,ylab="M2 standard forecast errors", xlab="time",col=cols[2],pch=20,type="l",ylim=fe.y.range)
    graphics::abline(h=0,lty=2)
  }
}




##################################
#### "exdqlmForecast" objects ####
##################################

#' \code{exdqlmForecast} objects
#'
#' \code{is.exdqlmForecast} tests if its argument is a \code{exdqlmForecast} object. 
#' 
#' @usage is.exdqlmForecast(x)
#'
#' @param x an \strong{R} object
#'
#' @export
is.exdqlmForecast = function(x){ return(methods::is(x,"exdqlmForecast")) }

#' Print Method for \code{exdqlmForecast} Objects
#'
#' @param x An \code{exdqlmForecast} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#'  y <- scIVTmag[1:100]
#'  model <- polytrendMod(1, stats::quantile(y, 0.85), 10)
#'  M0 <- exdqlmISVB(y, p0 = 0.85, model, df = c(0.98), dim.df = c(1),
#'                   gam.init = -3.5, sig.init = 15)
#'  M0.forecast = exdqlmForecast(start.t = 90, k = 10, m1 = M0)
#'  print(M0.forecast)
#' }
#'
print.exdqlmForecast <- function(x, ...) {
  #
  cat("k-step-ahead Quantile Forecasts of an exDQLM\n")
  cat("Number of Observations:", length(x$m1$y), "\n")
  cat("State Dimension:", length(x$m1$model$m0), "\n")  
  cat("Forecasts start at time index", x$start.t, "and forecast k =", x$k, "steps ahead\n")
  #
}

#' Summary Method for \code{exdqlmForecast} Objects
#'
#' @param object An \code{exdqlmForecast} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#'  y <- scIVTmag[1:100]
#'  model <- polytrendMod(1, stats::quantile(y, 0.85), 10)
#'  M0 <- exdqlmISVB(y, p0 = 0.85, model, df = c(0.98), dim.df = c(1),
#'                   gam.init = -3.5, sig.init = 15)
#'  M0.forecast = exdqlmForecast(start.t = 90, k = 10, m1 = M0)
#'  summary(M0.forecast)
#' }
#'
summary.exdqlmForecast <- function(object, ...) {
  #
  cat("k-step-ahead Quantile Forecasts of an exDQLM\n")
  cat("Number of Observations:", length(object$m1$y), "\n")
  cat("State Dimension:", length(object$m1$model$m0), "\n")  
  cat("Forecasts start at time index", object$start.t, "and forecast k =", object$k, "steps ahead\n")
  #
}

#' Plot Method for \code{exdqlmForecast} Objects
#'
#' @param x An \code{exdqlmForecast} object.
#' @param ... Additional arguments (unused).
#' 
#' @export
#' 
#' @examples
#' \donttest{
#'  y <- scIVTmag[1:100]
#'  model <- polytrendMod(1, stats::quantile(y, 0.85), 10)
#'  M0 <- exdqlmISVB(y, p0 = 0.85, model, df = c(0.98), dim.df = c(1),
#'                   gam.init = -3.5, sig.init = 15)
#'  M0.forecast = exdqlmForecast(start.t = 90, k = 10, m1 = M0)
#'  plot(M0.forecast)
#' }
#'
plot.exdqlmForecast <- function(x, ...) {
  aa = list(...)
  if(is.null(aa$cols)){cols=c("purple","magenta")}else{cols = aa$cols}
  if(is.null(aa$add)){add=FALSE}else{add=aa$add}
  
  y = x$m1$y
  p = dim(x$m1$model$GG)[1]
  TT = dim(x$m1$model$GG)[3]
  half.alpha = (1 - x$cr.percent)/2
  
  # filtered estimate for reference
  FF.start.t = matrix(x$m1$model$FF[,1:x$start.t], p, x$start.t)
  fm.start.t = matrix(x$m1$theta.out$fm[,1:x$start.t], p, x$start.t)
  qmap = colSums(matrix(FF.start.t*fm.start.t,p,x$start.t))
  fC.start.t = array(x$m1$theta.out$fC[,,1:x$start.t], c(p,p,x$start.t))
  temp.var = matrix(NA,p,x$start.t)
  for(t in 1:x$start.t){ temp.var[,t] = fC.start.t[,,t] %*% FF.start.t[,t] }
  qvar = colSums(FF.start.t * temp.var)
  qsd = sqrt(qvar)
  zlb = stats::qnorm(half.alpha)
  zub = stats::qnorm(x$cr.percent + half.alpha)
  qlb = qmap + zlb * qsd
  qub = qmap + zub * qsd
  # forecast estimates
  fqlb = x$ff + zlb * sqrt(x$fQ)
  fqub = x$ff + zub * sqrt(x$fQ)
  # filtered and forecasted quantiles & CrIs
  ts.xy = grDevices::xy.coords(y)
  if(!add){
    stats::plot.ts(y,xlim=c(ts.xy$x[x$start.t]-2*x$k*diff(ts.xy$x)[1],ts.xy$x[x$start.t]+x$k*diff(ts.xy$x)[1]),ylim=range(c(y,qlb,qub,fqlb,fqub)),type="l",ylab="quantile forecast",col="dark grey",xlab="time")
  }
  graphics::lines(ts.xy$x[1:x$start.t],qlb,col=cols[1],lty=3)
  graphics::lines(ts.xy$x[1:x$start.t],qub,col=cols[1],lty=3)
  graphics::lines(ts.xy$x[1:x$start.t],qmap,col=cols[1],lwd=1.5)
  graphics::lines(seq(from = ts.xy$x[x$start.t], by = diff(ts.xy$x)[1], length.out = x$k+1),c(qmap[x$start.t],x$ff),col=cols[2])
  graphics::lines(seq(from = ts.xy$x[x$start.t], by = diff(ts.xy$x)[1], length.out = x$k+1),c(qub[x$start.t],fqub),col=cols[2],lty=3)
  graphics::lines(seq(from = ts.xy$x[x$start.t], by = diff(ts.xy$x)[1], length.out = x$k+1),c(qlb[x$start.t],fqlb),col=cols[2],lty=3)
  
}



