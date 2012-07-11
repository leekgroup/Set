\name{bayesEst}
\alias{bayesEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtains the Bayes estimator for each atom considered.
}
\description{
  This function inputs the estimated quantities obtained from
  \code{errDisc} (EFD, EMD, EFDR, EMDR) and outputs
  the non-overlapping sets (atoms) in the Bayes estimator.
}
\usage{
bayesEst(EFD.EMD, w, type = "count", pen.for.feats =
0, pen.for.atoms = 0, atom.sizes = EFD.EMD[,"EFD"]/EFD.EMD[,"EFDR"],
ilim = 10, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EFD.EMD}{
    Matrix with estimated EFD, EMD, EFDR, EMDR for each of the
    non-overlapping sets (atoms) in the Bayes estimator
    (output of the \code{errDisc} function.)
  }
  \item{w}{
    Value of threshold used to determine whether or not an atom is in
    the Bayes estimator.
  }
  \item{type}{
    Either \code{"count"} or \code{"ratio"}, depending on whether the
    loss function is based on the EFD and EMD or on the EFDR and EMDR.
  }
  \item{pen.for.feats}{
    Penalty for the number of features in the Bayes estimator.
  }
  \item{pen.for.atoms}{
    Penalty for the number of atoms in the Bayes estimator.
  }
  \item{atom.sizes}{
    Vector of the number of genes in each atom.
  }
  \item{ilim}{
    Number of iterations performs to obtain Bayes estimator in ratio case.
  }
  \item{verbose}{
    If \code{TRUE}, displays intermediate messages for ratio case.
  }
}
\value{
  Vector of the atoms in the Bayes estimator.
  }
\references{
  Boca S.M., H. Corrada Bravo, B. Caffo, J.T. Leek, and G. Parmigiani (2010):
  A decision-theory approach to interpretable set analysis for
  high-dimensional
  data, \emph{JHU Biostat Working Paper 211},
  \url{http://www.bepress.com/jhubiostat/paper211/}
}
\author{
  Hector Corrada Bravo, Simina M. Boca, Jeffrey T. Leek
}
\seealso{postProb, errDisc}  
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

##now get Bayes estimator in various cases
bayesEst(EFD.EMD.res, w = 0.75, type = "count", atom.sizes =
sapply(sets, length))
bayesEst(EFD.EMD.res, w = 0.97, type = "count", atom.sizes =
sapply(sets, length))

bayesEst(EFD.EMD.res, w = 0.75, type = "count", atom.sizes =
sapply(sets, length), pen.for.feats = 0.2)
bayesEst(EFD.EMD.res, w = 0.97, type = "count", atom.sizes =
sapply(sets, length), pen.for.feats = 0.2)

bayesEst(EFD.EMD.res, w = 0.75, type = "count", atom.sizes =
sapply(sets, length), pen.for.atoms = 5)
bayesEst(EFD.EMD.res, w = 0.97, type = "count", atom.sizes =
sapply(sets, length), pen.for.atoms = 5)

bayesEst(EFD.EMD.res, w = 0.75, type = "ratio", atom.sizes =
sapply(sets, length))
}
