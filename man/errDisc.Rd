\name{errDisc}
\alias{errDisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtains estimates of the following quantities:
  expected number of false discoveries (EFD), the expected
  number of missed discoveries (EMD), the expected ratio of false
  discoveries (EFDR), and the expected ratio of missed discoveries
  (EMDR) for given atoms.
}
\description{
  This function inputs a list of non-overlapping sets (atoms) and a vector of
  estimated gene-level posterior probabilities of differential
  expression, and outputs a matrix with estimates of the following quantities:
  expected number of false discoveries (EFD), the expected
  number of missed discoveries (EFD), the expected fraction of false
  discoveries (EFDR), and the expected fraction of missed discoveries
  (EMDR) for the given atoms.
}
\usage{
errDisc(atoms, post.prob.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{atoms}{
    A list of non-overlapping sets (atoms).
  }
  \item{post.prob.genes}{
    Vector of the estimated gene-level probabilities of differential
    expression (output from \code{postProb}.)
  }
}
\value{
  Matrix with the estimated EFD, EMD, EFDR, and EMDR for each atom in \code{atoms} as columns.
}
\references{
  Boca S.M., H. Corrada Bravo, B. Caffo, J.T. Leek, and G. Parmigiani (2010):
  A decision-theory approach to interpretable set analysis for
  high-dimensional
  data, \emph{JHU Biostat Working Paper 211},
  \url{http://www.bepress.com/jhubiostat/paper211/}
}
\author{
  Simina M. Boca
}
\seealso{postProb}
\examples{
##first obtain estimates of the gene-level posterior probabilities
data(GenderStats)
data(ChromSets)

gene.stat <- (gene.stat - mean(gene.stat))/sd(gene.stat)
set.seed(70790707)
post.prob.genes <- postProb(gene.stat, null = rnorm(20 * length(gene.stat),
+ 0, 1), K = length(gene.stat))

##now run errDisc function
EFD.EMD.res <- errDisc(sets, post.prob.genes)
head(EFD.EMD.res)
}

