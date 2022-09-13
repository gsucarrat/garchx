\name{garchxSim}
\alias{garchxSim}

\title{Simulate from a GARCH-X model}

\description{
Simulate from a GARCH(q,p,r)-X model. Optionally, if \code{verbose=TRUE}, the conditional variance and innovations are also returned.
}

\usage{
garchxSim(n, intercept = 0.2, arch = 0.1, garch = 0.8, asym = NULL, xreg = NULL,
  innovations = NULL, backcast.values = list(), verbose = FALSE, as.zoo = TRUE,
  c.code = TRUE)
}

\arguments{
  \item{n}{integer}
  \item{intercept}{numeric}
  \item{arch}{\code{NULL} or numeric vector with the values of the ARCH-coefficients}
  \item{garch}{\code{NULL} or numeric vector with the values of the GARCH-coefficients}
  \item{asym}{\code{NULL} or numeric vector with the values of the asymmetry-coefficients}
  \item{xreg}{\code{NULL} or numeric vector with the values of the X-term}
  \item{innovations}{\code{NULL} or numeric vector with the innovations. If \code{NULL}, then standard normal innovations are generated with \code{\link{rnorm}}}
  \item{backcast.values}{\code{\link{list}} with backcast values}
  \item{verbose}{\code{logical}}
  \item{as.zoo}{\code{logical}. If \code{TRUE} (default), then the returned object is of class \code{\link{zoo}}}
  \item{c.code}{\code{logical}. If \code{TRUE} (default), then compiled C code is used for the recursion (faster). Otherwise R code is used (slower)}
}

\value{
a numeric vector or matrix with the simulated values.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{garchx}}, \code{\link{zoo}}
}

\examples{
##simulate from a garch(1,1):
y <- garchxSim(1000)

##simulate from a garch(1,1) with asymmetry/leverage:
yy <- garchxSim(1000, asym=0.1)

##simulate from a garch(1,1) w/user-provided backcast values:
yyy <- garchxSim(1000, backcast.values=list(z2=1, sigma2=0.5))

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
