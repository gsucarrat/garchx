\name{gdiff}
\alias{gdiff}

\title{ Difference a vector or a matrix, with special treatment of zoo objects}

\description{
Similar to the \code{\link{diff}} function from the base package, but \code{gdiff} enables padding (e.g. \code{NA}s or 0s) of the lost entries. Contrary to the \code{\link{diff}} function in the base package, however, the default in \code{gdiff} is to pad (with \code{NA}s). The \code{gdiff} function is particularly suited for \code{\link{zoo}} objects, since their indexing is retained
}

\usage{
gdiff(x, lag = 1, pad = TRUE, pad.value = NA)
}

\arguments{
  \item{x}{a numeric vector or matrix}
  \item{lag}{integer equal to the difference-length (the default is 1)}
  \item{pad}{logical. If \code{TRUE} (default), then the lost entries are padded with pad.value. If \code{FALSE}, then no padding is undertaken}
  \item{pad.value}{numeric, the pad-value}
}

\value{A vector or matrix with the differenced values}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\note{Empty}

\seealso{
\code{\link{diff}}, \code{\link{glag}}, \code{\link{lag}}
}

\examples{
##1st difference of a series:
x <- rnorm(5)
gdiff(x)

##1st difference with no padding:
gdiff(x, pad=FALSE)

##1st difference retaining the original zoo-index ordering:
gdiff(as.zoo(x))

##1st difference of a matrix:
y <- matrix(rnorm(8),4,2)
gdiff(y)

##2nd difference of the same matrix:
gdiff(y, lag=2)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
