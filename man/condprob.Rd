\name{condprob}
\title{Conditional Probabilities of Binary Data}
\usage{
condprob(x)
}
\alias{condprob}
\arguments{
 \item{x}{Matrix of binary data with rows corresponding to cases and
   columns corresponding to variables.}
}
\description{
  Returns a matrix containing the conditional probabilities
  \eqn{P(x_i=1|x_j=1)} where \eqn{x_i} corresponds to the \code{i}-th
  column of \code{x}. 
}
\author{Friedrich Leisch}
\keyword{binary}
\keyword{multivariate}
