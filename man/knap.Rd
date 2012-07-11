\name{knap}
\alias{knap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
yet another delightful little helper
}
\description{
yet more
}
\usage{
knap(rho, k, n, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function(rho, k, n, q)
  {
    N <- sum(n)
    c1 <- (1-k)/rho * (n-q)
    c2 <- k/(N-rho) * q

    obj <- c1-c2
    o <- order(obj/n)

    keep1 <- n[o] <= rho
    ii <- which(keep1)

    s <- cumsum(n[o][ii])
    keep <- which(s <= rho)

    if (s[max(keep)] < rho) {
      new.rho <- s[max(keep)+1]
      return(knap(new.rho, k, n, q))
    }

    keep <- ii[keep]
    ids <- o[keep]

    if (length(ids) == length(n)) {
      ids <- ids[-length(ids)]
    }

    res <- rep(0, length(n))
    res[ids] <- 1
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
