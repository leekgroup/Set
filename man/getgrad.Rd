\name{getgrad}
\alias{getgrad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
a third helper function
}
\description{
something
}
\usage{
getgrad(t, rho, k, n, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(t, rho, k, n, q)
  {
    s <- sum(n[t==1])
    sN <- sum(n)-s

    sq <- sum(q[t==1])
    sQ <- sum(q)-sq

    g1 <- (1-k)/s^2 * (sq*n - s*q)
    g2 <- k/sN^2 * (sQ*n - sN*q)
    g1+g2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
