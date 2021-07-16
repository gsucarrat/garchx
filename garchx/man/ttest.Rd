\name{ttest0}
\alias{ttest0}
\alias{waldtest0}

\title{T-tests and Wald-tests under nullity}

\description{
The permissible parameter-space of GARCH-models is bounded from below by 0. This means non-standard inference is required when one or more parameters are 0 under the null hypothesis. The functions \code{ttest0} and \code{waldtest0} perform t-tests and Wald-tests when one or more parameters is 0. In the latter test, the Wald-test, the critical values are obtained by simulation, see Francq and Thieu (2018).
}

\usage{
  ttest0(x, k = NULL)
  waldtest0(x, r = 0, R = NULL, level = c(0.1,0.05,0.01),
    vcov.type = NULL, quantile.type = 7, n = 20000)
}

\arguments{
  \item{x}{an object of class 'garchx'}
  \item{k}{\code{NULL} (default) or a vector of integers with the coefficients to test. If \code{NULL}, then all coefficients apart from the intercepts are tested}
  \item{r}{\code{vector} with restrictions}
  \item{R}{\code{NULL} (default) or a full-rank \code{matrix}. If \code{NULL}, then \code{R} is specified such that a test of all coefficients - apart from the intercept - is equal to the restriction \code{r}. If \code{length(r)==1}, then it is recycled so that its dimension match that of \code{R} }
  \item{level}{\code{vector} of significance levels whose critical values should be computed}
  \item{vcov.type}{\code{NULL} or a \code{character} that determines the type of coefficient-covariance to use, see \code{\link{vcov.garchx}} }
  \item{quantile.type}{integer, the algorithm used to compute the quantile, see \code{\link{quantile}} }
  \item{n}{\code{integer}, the numer of simulations used to estimate the critical values}
}

\details{
The \code{ttest0} function performs a t-test of coefficient \code{k} with 0 as null. Under this null the parameter is on the boundary of the admissible parameter space, and so the distribution is non-standard under the null. The function \code{ttest0} returns the result(s) of these non-standard t-test(s), see Francq and Thieu (2018). If \code{k=NULL}, the default, then a test for each coefficient apart from the intercept is undertaken.

The \code{waldtest0} function performs a Wald-test of the restrictions in \code{r}, when one or more of its elements are 0, see Francq and Thieu (2018).
}

\value{
  \item{ttest0:}{a \code{matrix} with the t-tests}
  \item{waldtest0:}{a \code{list} with the test-statistic and the critical values}
}

\references{
Christian Francq and Le Quien Thieu (2018): 'QML inference for volatility models with covariates', Econometric Theory, doi:10.1017/S0266466617000512
}

\author{
  Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{garchx}}, \code{\link{quantile}}, \code{\link{vcov.garchx}}, \code{\link{rmnorm}} }

\examples{
##simulate and estimate a garch(1,1):
set.seed(123)
y <- garchxSim(1000)
mymod <- garchx(y)

##t-tests:
ttest0(mymod)

##wald-test:
waldtest0(mymod)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
