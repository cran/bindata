\name{check.commonprob}
\title{Check Joint Binary Probabilities}
\alias{check.commonprob}
\usage{
check.commonprob(commonprob)
}
\arguments{
  \item{commonprob}{Matrix of pairwise probabilities.}
}
\description{
  The main diagonal elements \code{commonprob[i,i]} are interpreted as
  probabilities \eqn{p_{A_i}}{p(A_i)} that a binary variable \eqn{A_i}
  equals 1. The 
  off-diagonal elements \code{commonprob[i,j]} are the probabilites
  \eqn{p_{A_iA_j}}{p(A_iA_j)} that both \eqn{A_i} and \eqn{A_j} are 1.

  This programs checks some necessary conditions on these probabilities
  which must be fulfilled in order that a joint distribution of the
  \eqn{A_i} with the given probabilities can exist.

  The conditions checked are
  \deqn{0 \leq p_{A_i} \leq 1}{0 <= p(A_i) <= 1}
  
  \deqn{\max(0, p_{A_i} + p_{A_j} - 1) \leq p_{A_iA_j} \leq
    \min(p_{A_i}, p_{A_j}), i \neq j}{
    max(0, p(A_i)+p(A_j)-1) <= p(A_iA_j) <= min(p(A_i), p(A_j)), i != j}

   \deqn{p_{A_i} + p_{A_j} + p_{A_k} - p_{A_iA_j} -p_{A_iA_k} - p_{A_jA_k}
     \leq 1, i \neq j, i \neq k, j \neq k}{
     p(A_i)+p(A_j)+p(A_k)-p(A_iA_j)-p(A_iA_k)-p(A_jA_k) <= 1,
     i != j, i != k, j != k}
  
}
\value{
  \code{check.commonprob} returns \code{TRUE}, if all conditions are
  fulfilled. The attribute \code{"message"} of the return value contains
  some information on the errors that were found.
}
\references{Friedrich Leisch, Andreas Weingessel and Kurt Hornik
  (1998). On the generation of correlated artificial binary
  data. Working Paper Series, SFB ``Adaptive Information Systems and
  Modelling in Economics and Management Science'', Vienna University of
  Economics, \url{http://www.wu-wien.ac.at/am}
}
\author{Andreas Weingessel}
\seealso{
  \code{\link{simul.commonprob}},
  \code{\link{commonprob2sigma}}}
\examples{
check.commonprob(cbind(c(0.5, 0.4), c(0.4, 0.8)))

check.commonprob(cbind(c(0.5, 0.25), c(0.25, 0.8)))

check.commonprob(cbind(c(0.5, 0, 0), c(0, 0.5, 0), c(0, 0, 0.5)))
}
\keyword{distribution}
\keyword{multivariate}
