\name{commonprob2sigma}
\title{Calculate a Covariance Matrix for the Normal Distribution out of
  a Matrix of Joint Probabilites.}
\usage{
commonprob2sigma(commonprob, simulvals)
}
\alias{commonprob2sigma}
\arguments{
  \item{commonprob}{Matrix of Pairwise Probabilities.}
  \item{simulvals}{Array received by \code{simul.commonprob}.}
}
\description{
  Computes a covariance matrix for a normal distribution which
  corresponds to a binary distribution with marginal probabilites given
  by \code{diag(commonprob)} and pairwise probabilities given by
  \code{commonprob}. 

  For the simulations the values of \code{simulvals} are used.

  If a non-valid covariance matrix is the result, the program stops with
  an error in the case of NA arguments and yields are warning message if
  the matrix is not positive definite.
}
\value{
A covariance matrix is returned with the same dimensions as
\code{commonprob}. 
}
\references{Friedrich Leisch, Andreas Weingessel and Kurt Hornik
  (1998). On the generation of correlated artificial binary
  data. Working Paper Series, SFB ``Adaptive Information Systems and
  Modelling in Economics and Management Science'', Vienna University of
  Economics, \url{http://www.wu-wien.ac.at/am}
}
\author{Friedrich Leisch}
\seealso{\code{\link{simul.commonprob}}}
\examples{
m <- cbind(c(1/2,1/5,1/6),c(1/5,1/2,1/6),c(1/6,1/6,1/2))
sigma <- commonprob2sigma(m)
}
\keyword{binary}
\keyword{multivariate}

