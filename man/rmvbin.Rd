\name{rmvbin}
\title{Multivariate Binary Random Variates}
\usage{
rmvbin(n, margprob, commonprob=diag(margprob),
       bincorr=diag(length(margprob)),
       sigma=diag(length(margprob)),
       colnames=NULL, simulvals=NULL)
}
\alias{rmvbin}
\arguments{
 \item{n}{Number of observations.}
 \item{margprob}{Margin probabilities that the components are 1.}
 \item{commonprob}{Matrix of probabilities that components \code{i} and
   \code{j} are simultaneously 1.}
 \item{bincorr}{Matrix of binary correlations.}
 \item{sigma}{Covariance matrix for the normal distribution.}
 \item{colnames}{Vector of column names for the resulting observation
   matrix.}
 \item{simulvals}{Result from \code{simul.commonprob(), a default data
     array is automatically loaded if this argument is omitted.}
 }
\description{
  Creates correlated multivariate binary random variables by
  thresholding a normal distribution. The correlations of the components
  can be specified either as common probabilities, correlation matrix  of the
  binary distribution, or covariance matrix of the normal
  distribution. Hence, only one of the arguments \code{commonprob},
  \code{bincorr} and \code{sigma} may be specified. Default are
  uncorrelated components.
  
  \code{n} samples from a multivariate normal distribution with mean and
  variance chosen in order to get the desired margin and common
  probabilities are sampled. Negative values are converted to 0,
  positive values to 1.
}
\references{Friedrich Leisch, Andreas Weingessel and Kurt Hornik
  (1998). On the generation of correlated artificial binary
  data. Working Paper Series, SFB ``Adaptive Information Systems and
  Modelling in Economics and Management Science'', Vienna University of
  Economics, \url{http://www.wu-wien.ac.at/am}
}
\author{Friedrich Leisch}
\seealso{\code{\link{commonprob2sigma}},\code{\link{check.commonprob}},
  \code{\link{simul.commonprob}}}
\examples{
# uncorrelated columns:
rmvbin(10, margprob=c(0.3,0.9))

# correlated columns
m <- cbind(c(1/2,1/5,1/6),c(1/5,1/2,1/6),c(1/6,1/6,1/2))
rmvbin(10,commonprob=m)

# same as the second example, but faster if the same probabilities are
# used repeatedly (coomonprob2sigma rather slow)
sigma <- commonprob2sigma(m)
rmvbin(10,margprob=diag(m),sigma=sigma)
}
\keyword{binary}
\keyword{multivariate}
\keyword{distribution}




