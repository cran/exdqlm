% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polytrendMod.R
\name{polytrendMod}
\alias{polytrendMod}
\title{Create an n-th order polynomial exDQLM component}
\usage{
polytrendMod(order, m0, C0)
}
\arguments{
\item{order}{The order of the polynomial model.}

\item{m0}{Prior mean of the state vector. Default is \code{m0 = rep(0,order)}.}

\item{C0}{Prior covariance of the state vector. Default is \code{C0 = 1e3*diag(order)}.}
}
\value{
A list of the following:
\itemize{
\item FF - Observational vector.
\item GG - Evolution matrix.
\item m0 - Prior mean of the state vector.
\item C0 - Prior covariance of the state vector.
}
}
\description{
The function creates an n-th order polynomial exDQLM component.
}
\examples{
# create a second order polynomial component
trend.comp = polytrendMod(2,rep(0,2),10*diag(2))
}
