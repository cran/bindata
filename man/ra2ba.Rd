\name{ra2ba}
\title{Convert Real Valued Array to Binary Array}
\usage{
ra2ba(x)
}
\alias{ra2ba}
\arguments{
 \item{x}{Array of arbitrary dimension}
}
\description{
  Converts all values of the real valued array \code{x} to binary values
  by thresholding at 0.
}
\author{Friedrich Leisch}
\examples{
x <- array(rnorm(10), dim=c(2,5))
ra2ba(x)
}
\keyword{array}
