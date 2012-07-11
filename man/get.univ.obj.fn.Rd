\name{get.univ.obj.fn}
\alias{get.univ.obj.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
another helper function
}
\description{
more stuff
}
\usage{
get.univ.obj.fn(t, p, k, n, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function(t, p, k ,n ,q)
  {
    a1 <- p\%*\%(n-q)
    b1 <- t\%*\%(n-q)
    c1 <- p\%*\%n
    d1 <- t\%*\%n

    a2 <- -p\%*\%q
    b2 <- sum(q) - t\%*\%q
    c2 <- -p\%*\%n
    d2 <- sum(n) - t\%*\%n

    function(sigma)
      {
        res <- (1-k) * (sigma*a1 + b1) / (sigma*c1 + d1)
        res <- res + k * (sigma*a2 + b2) / (sigma*c2 + d2)
        res
      }

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
