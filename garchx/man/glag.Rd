\name{glag}
\alias{glag}

\title{ Lag a vector or a matrix, with special treatment of \code{\link[zoo]{zoo}} objects}

\description{
Similar to the \code{\link{lag}} function from the \code{stats} package, but \code{glag} enables padding (e.g. \code{NA}s or 0s) of the lost entries. Contrary to the \code{\link{lag}} function in the stats package, however, the default in \code{glag} is to pad (with \code{NA}s). The \code{glag} is particularly suited for \code{\link[zoo]{zoo}} objects, since their indexing is retained
}

\usage{
glag(x, k = 1, pad = TRUE, pad.value = NA)
}

\arguments{
  \item{x}{a numeric vector or matrix}
  \item{k}{integer equal to the lag (the default is 1)}
  \item{pad}{logical. If TRUE (default), then the lost entries are padded with pad.value. If FALSE, then no padding is undertaken}
  \item{pad.value}{the pad-value}
}

\value{A vector or matrix with the lagged values}

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{lag}}, \code{\link{gdiff}}, \code{\link{diff}}
}

\examples{
##lag series with NA for the missing entries:
x <- rnorm(5)
glag(x)

##lag series with no padding:
x <- rnorm(5)
glag(x, pad=FALSE)

##lag series and retain the original zoo-index ordering:
x <- as.zoo(rnorm(5))
glag(x)

##lag two periods:
glag(x, k=2)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
