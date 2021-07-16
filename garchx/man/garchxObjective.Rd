\name{garchxObjective}
\alias{garchxObjective}
\alias{garchxRecursion}

\title{Auxiliary functions}

\description{
Auxiliary functions used in estimation. Not intended for the average user
}

\usage{
  garchxObjective(pars, aux)
  garchxRecursion(pars, aux)
}

\arguments{
  \item{pars}{numeric vector of parameters}
  \item{aux}{\code{\link{list}} created by \code{\link{garchx}}}
}

\value{
  \item{garchxObjective:}{value of the objective function}
  \item{garchxRecursion:}{fitted conditional variance}
}

\author{
  Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
\code{\link{garchx}}, \code{\link{fitted.garchx}}, \code{\link{residuals.garchx}} 
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
