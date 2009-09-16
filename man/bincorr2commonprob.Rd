\name{bincorr2commonprob}
\alias{bincorr2commonprob}
\title{Convert Binary Correlation Matrix to Matrix of Joint
  Probabilities}
\description{
  Compute a matrix of common probabilities for a binary random vector
  from given marginal probabilities and correlations.
}
\usage{
bincorr2commonprob(margprob, bincorr)
}
\arguments{
  \item{margprob}{vector of marginal probabilities.}
  \item{bincorr}{matrix of binary correlations.}
}
\value{
  The matrix of common probabilities.  This has the probabilities that
  variable \eqn{i} equals 1 in element \eqn{(i,i)}, and the joint
  probability that variables \eqn{i} and \eqn{j} both equal 1 in element
  \eqn{(i,j)} (if \eqn{i \ne j}{i != j}).
}
\references{Friedrich Leisch, Andreas Weingessel and Kurt Hornik
  (1998). On the generation of correlated artificial binary
  data. Working Paper Series, SFB ``Adaptive Information Systems and
  Modelling in Economics and Management Science'', Vienna University of
  Economics, \url{http://www.wu-wien.ac.at/am}
}
\author{Friedrich Leisch}
\seealso{
  \code{\link{commonprob2sigma}},
  \code{\link{simul.commonprob}}.
}
\keyword{distribution}
\keyword{multivariate}
