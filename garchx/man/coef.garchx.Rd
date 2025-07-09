\name{coef.garchx}
\alias{coef.garchx}
\alias{confint.garchx}
\alias{fitted.garchx}
\alias{logLik.garchx}
\alias{nobs.garchx}
\alias{predict.garchx}
\alias{print.garchx}
\alias{quantile.garchx}
\alias{residuals.garchx}
\alias{toLatex.garchx}
\alias{vcov.garchx}

\title{Extraction functions for 'garchx' objects}

\description{Extraction functions for objects of class 'garchx'}

\usage{
  \method{coef}{garchx}(object, ...)
  \method{confint}{garchx}(object, parm, level = 0.95, ...)
  \method{fitted}{garchx}(object, as.zoo = TRUE, ...)
  \method{logLik}{garchx}(object, ...)
  \method{nobs}{garchx}(object, ...)
  \method{predict}{garchx}(object, n.ahead = 10, newxreg = NULL,
    newindex = NULL, n.sim = NULL, verbose = FALSE, ...)
  \method{print}{garchx}(x, ...)
  \method{quantile}{garchx}(x, probs=0.025, names = TRUE, type = 7, as.zoo = TRUE, ...)
  \method{residuals}{garchx}(object, as.zoo = TRUE, ...)
  \method{toLatex}{garchx}(object, digits = 4, ...)
  \method{vcov}{garchx}(object, vcov.type = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'garchx'}
  \item{x}{an object of class 'garchx'}
  \item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered}
  \item{level}{\code{numeric} value between 0 and 1 (the confidence level required)}
  \item{as.zoo}{logical. If \code{TRUE}, then the returned result is of class \code{\link[zoo]{zoo}} }
  \item{n.ahead}{\code{integer} that determines how many steps ahead predictions should be generated}
  \item{newxreg}{\code{vector} or \code{matrix} with the out-of-sample regressor values}
  \item{newindex}{\code{\link[zoo]{zoo}}-index for the out-of-sample predictions. If \code{NULL} (default), then \code{1:n.ahead} is used}
  \item{n.sim}{\code{NULL} or an \code{integer}, the number of simulations}
  \item{verbose}{\code{logical}. If \code{TRUE}, then the simulations - in addition to the predictions - are returned}
  \item{probs}{\code{vector} of probabilities}
  \item{names}{\code{logical}, whether to return names or not}
  \item{type}{\code{integer} that determines the algorithm used to compute the quantile, see \code{\link{quantile}} }
  \item{digits}{\code{integer}, the number of digits in the printed LaTeX code}
  \item{vcov.type}{\code{NULL} or \code{character} that is (partially) matched to \code{"ordinary"}, \code{"robust"} or \code{"hac"}. The robust coefficient-covariance (\code{"robust"}) is that of Francq and Thieu (2019). The Heteroscedasticity and Autocorrelation Consistent (\code{"hac"}) variance-covariance estimates the variance of the score (the 'meat') based on Theorem 2.2 in De Jong and Davidson (2000). For the details of the kernel weights (Bartlett) and bandwidth, see the code in the \code{\link{vcov.garchx}} function }
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{numeric vector containing parameter estimates}
  \item{confint:}{A matrix lower and upper confidence limits for each parameter}
  \item{fitted:}{fitted conditional variance}
  \item{logLik:}{log-likelihood (normal density)}
  \item{nobs:}{the number of observations used in the estimation}
  \item{predict:}{a \code{vector} with the predictions (\code{verbose=FALSE}), or a \code{matrix} with both the predictions and the simulations (\code{verbose=TRUE}) }
  \item{print:}{print of the estimation results}
  \item{quantile:}{the fitted quantiles, i.e. the conditional standard deviation times the empirical quantile of the standardised innovations}
  \item{residuals:}{standardised residuals}
  \item{vcov:}{coefficient variance-covariance matrix}
}

\references{
Christian Francq and Le Quien Thieu (2019): 'QML inference for volatility models with covariates', Econometric Theory 35, pp. 37-72, \doi{10.1017/S0266466617000512}
Robert M. de Jong and James Davidson (2000): 'Consistency of Kernel Estimators of Heteroscedastic and Autocorrelated Covariance Matrices', Econometrica 68, pp. 407-423
}

\author{
  Genaro Sucarrat, \url{https://www.sucarrat.net/}
}

\seealso{\code{\link{garchx}}, \code{\link{garchxSim}}, \code{\link[zoo]{zoo}} }

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

##extract coefficient variance-covariance matrix:
vcov(mymod)

##generate predictions:
predict(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
