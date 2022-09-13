\name{garchx}
\alias{garchx}

\title{Estimate a GARCH-X model}

\description{
Quasi Maximum Likelihood (ML) estimation of a GARCH(q,p,r)-X model, where q is the GARCH order, p is the ARCH order, r is the asymmetry (or leverage) order and 'X' indicates that covariates can be included. Note that the underlying estimation theory assumes the covariates are stochastic. The estimation procedure will, in general, provide consistent estimates when the standardised innovations are not normal or independent (or both), see Francq and Thieu (2018). 
}

\usage{
garchx(y, order = c(1,1), arch = NULL, garch = NULL, asym = NULL,
  xreg = NULL, vcov.type = c("ordinary", "robust"),
  initial.values = NULL, backcast.values = NULL, lower = 0,
  upper = +Inf, control = list(), hessian.control = list(),
  solve.tol = .Machine$double.eps, estimate = TRUE, c.code = TRUE,
  penalty.value = NULL, sigma2.min = .Machine$double.eps,
  objective.fun = 1, turbo = FALSE)
}

\arguments{
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{order}{integer vector of length 1, 2 or 3, for example \code{c(1,1,1)}. The first entry controls the GARCH order, the second the ARCH order and the third the ASYM (asymmetry/leverage) order}
  \item{arch}{\code{NULL} or numeric vector containing the ARCH-terms to include. Note: If not \code{NULL}, then the value of the \code{ARCH} argument overrides the value of the first entry in the \code{order} argument}
  \item{garch}{\code{NULL} or numeric vector containing the GARCH-terms to include. Note: If not \code{NULL}, then the value of the \code{GARCH} argument overrides the value of the second entry in the \code{order} argument}
  \item{asym}{\code{NULL} or numeric vector containing the ASYM-terms (asymmetry/leverage terms) to include. Note: If not \code{NULL}, then the value of the \code{ASYM} argument overrides the value of the third entry in the \code{order} argument}
  \item{xreg}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{vcov.type}{\code{character}, either \code{"ordinary"} or \code{"robust"}, see \code{\link{vcov.garchx}} }
  \item{initial.values}{\code{NULL} or a numeric vector with the initial parameter values passed on to the optimisation routine, \code{\link{nlminb}}. If \code{NULL}, the default, then the values are chosen automatically}
  \item{backcast.values}{\code{NULL} or a non-negative numeric. The backcast value is used to initiate the forward recursion of the conditional variance. If \code{NULL} (default), then the value is chosen automatically (currently the average of \code{y} squared is used). If \code{backcast.values} is a non-negative numeric, then the initial recursion values are all set to this value}
  \item{lower}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}}}
  \item{upper}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}}}
  \item{control}{a \code{\link{list}} passed on to the \code{control} argument of \code{\link{nlminb}}}
  \item{hessian.control}{a \code{\link{list}} passed on to the \code{control} argument of \code{\link{optimHess}}}
  \item{solve.tol}{numeric value passed on to the \code{tol} argument of \code{\link{solve}}, which is called whenever the coefficient variance-coariance matrix is computed. The value controls the toleranse for detecting linear dependence between columns when inverting a matrix}
  \item{estimate}{logical, if \code{TRUE} then estimation is carried out. If \code{FALSE}, then the \code{initial.values} are used}
  \item{c.code}{logical, if \code{TRUE} then compiled C code is used in the forward recursion}
  \item{penalty.value}{\code{NULL} (default) or a numeric value. If \code{NULL}, then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in \code{NA} and/or \code{+/-Inf} values. The \code{penalty.value} is the value returned by the objective function \code{\link{garchxObjective}} in the presence of \code{NA} or \code{+/-Inf} values}
  \item{sigma2.min}{numeric with default \code{.Machine$double.eps}. To avoid taking taking the log of a very small value when computing the log-likelihood, \code{sigma2.min} is used as the lower bound of the fitted conditional variances, see the code of \code{\link{garchxObjective}} }
  \item{objective.fun}{numeric, either 1 or 0}
  \item{turbo}{logical. If \code{FALSE} (default), then the coefficient variance-covariance is computed during estimation, and the fitted values and residuals are attached to the returned object. If \code{TRUE}, then these operations are skipped, and hence estimation is faster. Note, however, that if \code{turbo} is set to \code{TRUE}, then the coefficient-covariance, fitted values and residuals can still be extracted subsequent to estimation with \code{\link{vcov.garchx}}, \code{\link{fitted.garchx}} and \code{\link{residuals.garchx}}, respectively }
}

\value{A list of class 'garchx'}

\references{
Christian Francq and Le Quien Thieu (2018): 'QML inference for volatility models with covariates', Econometric Theory, doi:10.1017/S0266466617000512\cr
Christian Francq and Jean-Michel Zakoian (2019): 'GARCH Models', 2nd Edition, Wiley
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{garchxSim}}, \code{\link{nlminb}}, \code{\link{optimHess}}, \code{\link{coef.garchx}}
}

\examples{
##simulate from a garch(1,1):
set.seed(123)
y <- garchxSim(1000)

##estimate garch(1,1) model:
mymod <- garchx(y)

##print estimation results:
print(mymod)

##extract coefficients:
coef(mymod)

##extract and store conditional variances:
sigma2hat <- fitted(mymod)

##extract log-likelihood:
logLik(mymod)

##extract and store standardised residuals:
etahat <- residuals(mymod)

##extract variance-covariance matrix:
vcov(mymod)

##generate predictions:
predict(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
