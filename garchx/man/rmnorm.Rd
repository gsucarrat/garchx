\name{rmnorm}
\alias{rmnorm}

\title{
Random number generation from the multivariate normal distribution
}
\description{
This function is a speed-optimised version of the \code{\link{rmnorm}} function from the mnormt package of Adelchi Azzalini (2013).
}

\usage{rmnorm(n, mean = NULL, vcov = 1)}

\arguments{
  \item{n}{integer, the number of observations to generate}
  \item{mean}{numeric vector, i.e. the mean values}
  \item{vcov}{numeric matrix, i.e. the variance-covariance matrix}
}

\value{A matrix of n rows}

\references{
Adelchi Azzalini (2013): 'mnormt: The multivariate normal and t distributions', R package version 1.4-7, \url{https://CRAN.R-project.org/package=mnormt}
}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\examples{
##generate from univariate standardised normal:
z1 <- rmnorm(100)

##generate from bivariate, independent standardised normal:
z2 <- rmnorm(100, vcov=diag(c(1,1)))

##generate from bivariate, dependent standardised normal:
z3 <- rmnorm(100, vcov=cbind(c(1,0.3),c(0.3,1)))
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
