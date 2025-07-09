\name{refit}
\alias{refit}
\alias{refit.garchx}

\title{Refit a model to new data}

\description{Refit a model to new data, typically for cross-validation purposes. Re-estimation (\code{reestimate}) is optional (the default is \code{FALSE}).}

\usage{
##generic:
refit(object, ...)

##S3 method for 'garchx' objects:
\method{refit}{garchx}(object, newy = NULL, newxreg = NULL,
  backcast.value = NULL, reestimate = FALSE, ...)

}

\arguments{
\item{object}{an object of class \code{\link{garchx}} }
\item{newy}{vector, the new 'y' data, see \code{\link{garchx}} }
\item{newxreg}{the new 'xreg' data, if any, see \code{\link{garchx}} }
\item{backcast.value}{\code{NULL} or a non-negative numeric, see \code{\link{garchx}} }
\item{reestimate}{logical. If \code{FALSE} (default), then the estimates from \code{object} are used on the new data. If \code{TRUE}, then the model is re-estimated using the new data }
\item{...}{further arguments passed to or from other methods}
}

\value{A list of class 'garchx'}

\details{
\code{refit.garchx} is a convenience function to facilitate cross-validation and related analyses. 
}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\seealso{
\code{\link{garchx}}
}

\examples{
##simulate from a garch(1,1):
set.seed(123)
y <- garchxSim(1000)

##estimate garch(1,1) model:
mymod <- garchx(y)

##new data (e.g. 'out-of-sample' or 'test' data):
yy <- garchxSim(100)

##apply the estimates of 'mymod' on yy data:
refit(mymod, newy=yy)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}