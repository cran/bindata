\name{condprob}
\alias{condprob}
\title{Conditional Probabilities of Binary Data}
\usage{
condprob(x)
}
\arguments{
  \item{x}{matrix of binary data with rows corresponding to cases and
    columns corresponding to variables.}
}
\description{
  Returns a matrix containing the conditional probabilities
  \eqn{P(x_i=1|x_j=1)} where \eqn{x_i} corresponds to the \code{i}-th
  column of \code{x}. 
}
\author{Friedrich Leisch}
\keyword{distribution}
