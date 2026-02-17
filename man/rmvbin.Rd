\name{rmvbin}
\alias{rmvbin}
\title{Multivariate Binary Random Variates}
\description{
  Creates correlated multivariate binary random variables by
  thresholding a normal distribution.  The correlations of the
  components can be specified either as common probabilities,
  correlation matrix of the binary distribution, or covariance matrix of
  the normal distribution.
}
\usage{
rmvbin(n, margprob, commonprob=diag(margprob),
       bincorr=diag(length(margprob)),
       sigma=diag(length(margprob)),
       colnames=NULL, simulvals=NULL)
}
\arguments{
  \item{n}{number of observations.}
  \item{margprob}{margin probabilities that the components are 1.}
  \item{commonprob}{matrix of probabilities that components \code{i} and
    \code{j} are simultaneously 1.}
  \item{bincorr}{matrix of binary correlations.}
  \item{sigma}{covariance matrix for the normal distribution.}
  \item{colnames}{vector of column names for the resulting observation
    matrix.}
  \item{simulvals}{result from \code{simul.commonprob}, a default data
    array is automatically loaded if this argument is omitted.}
}
\details{
  Only one of the arguments \code{commonprob}, \code{bincorr} and
  \code{sigma} may be specified.  Default are uncorrelated components.
  
  \code{n} samples from a multivariate normal distribution with mean and
  variance chosen in order to get the desired margin and common
  probabilities are sampled.  Negative values are converted to 0,
  positive values to 1.
}
\references{
  \bibshow{bindata::Leisch+Weingessel+Hornik:1998}  
}
\author{Friedrich Leisch}
\seealso{
  \code{\link{commonprob2sigma}},
  \code{\link{check.commonprob}},
  \code{\link{simul.commonprob}}
}
\examples{
## uncorrelated columns:
rmvbin(10, margprob=c(0.3,0.9))

## correlated columns
m <- cbind(c(1/2,1/5,1/6),c(1/5,1/2,1/6),c(1/6,1/6,1/2))
rmvbin(10, commonprob=m)

## same as the second example, but faster if the same probabilities are
## used repeatedly (commonprob2sigma rather slow)
sigma <- commonprob2sigma(m)
rmvbin(10, margprob = diag(m), sigma = sigma)

## The default 'simulvals' may not work for very small marginal
## probabilities.  E.g., for 
p1 <- 0.0206
p2 <- 0.0318
m <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
## rmvbin(1, margprob = c(p1,p2), bincorr = m) fails with 'Extrapolation
## occurred ... margprob and commonprob not compatible?'
## (reported by <Ping.Mahling@pfizer.com> in 2025-12).
## Using e.g.
s <- simul.commonprob(margprob = c(p1, p2), corr = seq(-1, 1, by = 0.05))
## makes things work:
rmvbin(10, margprob = c(p1,p2), bincorr = m, simulvals = s)
}
\keyword{multivariate}
\keyword{distribution}
