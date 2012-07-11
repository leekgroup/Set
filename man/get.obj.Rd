\name{get.obj}
\alias{get.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
helper function
}
\description{
stuff
}
\usage{
get.obj(t, k, n, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
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
function(t, k, n, q)
  {
    res <- (1-k) * sum((n-q)[t==1])
    res1 <- res / sum(n[t==1])
    res <- k * sum(q[t==0])
    res2 <- res / sum(n[t==0])
    res1 + res2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
