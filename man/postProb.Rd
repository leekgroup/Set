\name{postProb}
\alias{postProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates the posterior probabilities that genes are differentially expressed.
}
\description{
  This function inputs gene-level statistics and returns a vector of estimated
  posterior probabilities of differential expression for each gene.
}
\usage{
postProb(obs.stat, null, B = NULL, pos = FALSE, p0 = 1,
K = length(obs.stat),dim.basis = 10,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs.stat}{
    Vector of observed gene-level statistics.
  }
  \item{null}{
    Vector of null gene-level statistics or a character (such as
    "rnorm") indicating the random generating mechanism for the null statistics.
  }
  \item{B}{
    Number of simulations under the null.
  }
  \item{pos}{
    If set to \code{TRUE}, applies \code{log} transformation to the
    gene-level statistics.
  }
  \item{p0}{
    Value (between \code{0} and \code{1}) for the prior probability that
    a gene is not differentially expressed.
  }
  \item{K}{
    Number of intervals used in probability estimation. A good value in
    practice seems to be the length of \code{obs.stat}.
  }
  \item{dim.basis}{
    Basis dimension (\code{k}) in \code{gam}.
  }
  \item{...}{
    If \code{null} is a character, then pass the parameters of the
    relevant random generating mechanism.
  }
}
\value{
  Vector of the estimated gene-level probabilities of differential expression.
}
\references{
  Storey, J., J. Akey, and L. Kruglyak (2005): Multiple locus linkage
  analysis of genomewide expression in yeast, \emph{PLoS Biology}, 3.
  
  Boca S.M., H. Corrada Bravo, B. Caffo, J.T. Leek, and G. Parmigiani (2010):
  A decision-theory approach to interpretable set analysis for
  high-dimensional
  data, \emph{JHU Biostat Working Paper 211},
  \url{http://www.bepress.com/jhubiostat/paper211/}
}
\author{
  John D. Storey, Jeffrey T. Leek, Simina M. Boca
}

\examples{
set.seed(10831)
stat <- rnorm(1000,sd=sqrt(1/10),mean=c(rep(1,50),rep(0,950)))
set.seed(23134)
stat0 <- matrix(rnorm(1000*20,sd=sqrt(1/10)),nrow=1000)

post.prob.genes.1 <-
    postProb(stat,
             null = stat0,
             K = 1000)

set.seed(23134)
post.prob.genes.2 <-
    postProb(stat,
             null = "rnorm",
             B = 20,
             K = 1000,
             mean = 0, sd = sqrt(1/10))

##The results are the same, since in both case we are considering 20
##simulations of 1000 null statistics from a N(0, 1/10) distribution,
##using the same seed!
identical(post.prob.genes.1, post.prob.genes.2)
}
