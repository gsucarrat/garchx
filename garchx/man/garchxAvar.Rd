\name{garchxAvar}
\alias{garchxAvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Asymptotic Coefficient Covariance}

\description{Compute the asymptotic coefficient-covariance of a GARCH(q,p,r)-X model by simulation. Note that the principles of how to use the \code{arch}, \code{garch}, \code{asym} and \code{xreg} arguments are the same as those of \code{\link{garchx}} }

\usage{garchxAvar(pars, arch = NULL, garch = NULL, asym = NULL,
  xreg = NULL, vcov.type = c("ordinary", "robust"),
  innovations = NULL, Eeta4 = NULL, n = 1e+06, objective.fun = 1,
  seed = NULL)
}

\arguments{
  \item{pars}{\code{vector} of parameters of length 1 or more. The first component contains the coefficient-value of the intercept, the next component(s) the ARCH-coefficient(s), and so on.}
  \item{arch}{\code{NULL} or integer \code{vector} with the lags of the ARCH-terms to include. Works in the same way as the \code{arch} argument in the \code{\link{garchx}} function}
  \item{garch}{\code{NULL} or integer \code{vector} with the lags of the GARCH-terms. Works in the same way as the \code{garch} argument in the \code{\link{garchx}} function}
  \item{asym}{\code{NULL} or integer \code{vector} with the lags of the asymmetry terms to include. Works in the same way as the \code{asym} argument in the \code{\link{garchx}} function}
  \item{xreg}{\code{NULL}, or a \code{vector} or \code{matrix} with the covariates of the model. Works in the same way as the \code{xreg} argument in the \code{\link{garchx}} function}
  \item{vcov.type}{\code{character} that determines the type of coefficient-covariance}
  \item{innovations}{\code{NULL} or a \code{vector} with the standardised innovations to use. If \code{NULL}, then the innovations are standard normal}
  \item{Eeta4}{\code{numeric}, the fourth moment of the innovations. If \code{NULL}, then the value is estimated internally. Note: The value of Eeta4 is only used if \code{vcov.type = "ordinary"}, otherwise it is ignored }
  \item{n}{integer, the number of observations to use in the simulations}
  \item{objective.fun}{integer equal to 1 or 0 that determines the type of objective function to use, see the code of \code{\link{garchxObjective}} }
  \item{seed}{\code{NULL} or an integer that sets the seed (the value is passed on to \code{\link{set.seed}}. Useful for reproducibility}
}

\value{A \code{matrix}}

\references{
Christian Francq and Le Quien Thieu (2018): 'QML inference for volatility models with covariates', Econometric Theory, doi:10.1017/S0266466617000512
Christian Francq and Jean-Michel Zakoian (2019): 'GARCH Models', 2nd Edition, Wiley
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{garchx}}, \code{\link{garchxSim}}, \code{\link{vcov.garchx}}
}

\examples{

##asymptotic coefficient-covariance of a garch(1,1)
##note: the estimate is rough, since n is small
intercept <- 0.2
alpha <- 0.1
beta <- 0.8
pars <- c(intercept, alpha, beta)
seed <- 123 #for reproducibility
garchxAvar(pars, arch=1, garch=1, n=10000, seed=seed)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
