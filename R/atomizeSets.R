#' A function to make a mosaic-style plot for sets and their atoms.
#'
#' This function allows you to input a list of sets of features and visualize their atoms via a mosaic-style plot.
#' @param sets List representing sets of features.
#' @param values
#' @param col
#' @param bgCol
#' @param orderBy
#' @param setLab
#' @param atomAvLab
#' @param ylab 
#' @export
#' @examples
#' sets <- list(A=1:50, B=31:130, C=131:160, D=c(40:50, 161:260))
#' atomizeSets(sets)

atomizeSets <- function(sets, values=NULL, col="red", bgCol="ivory3", orderBy="none", setLab=NULL, atomAvLab=NULL, ylab="Atom size")
{
  ##get atoms
  atoms <- enumAtoms(sets)
  nrAtoms <- length(atoms)
  ##get atom sizes
  sizesAtoms <- sapply(atoms, length)
  
  ##get average values in atoms
  if(!is.null(values))
  {
    ##make values between 0 and 1
    if(sum(values < 0)>=1)
    {
      warning("Some values are less than 0. Mean values per atom will not be plotted.")
      values <- NULL
    }
    if(sum(values > 1)>=1)
    {
      values <- values/max(values)
    }
    avAtoms <- sapply(atoms, function(a,v){mean(v[a])}, values)
    avAtoms[avAtoms < 0.01] <- 0.01
    if(orderBy=="mean")
    {
      atoms <- atoms[order(avAtoms)]
      sizesAtoms <- sizesAtoms[order(avAtoms)]
      avAtoms <- sort(avAtoms)
    }
  }
  
  if(orderBy=="size")
  {
    atoms <- atoms[order(sizesAtoms)]
    if(!is.null(values)) avAtoms <- avAtoms[order(sizesAtoms)]
    sizesAtoms <- sort(sizesAtoms)
  }
  
  sizeTotal <- sum(sizesAtoms)
  cumSizesAtoms <- cumsum(sizesAtoms)
  cumSizesAtoms <- c(0, cumSizesAtoms)
  
  ##get set names
  setNames <- names(sets)
  if(is.null(setNames)) 
  {
    names(sets) <- setNames <- as.character(1:length(sets))
  }
  nrSets <- length(setNames)
  
  plot(c(0,nrSets*0.5+(nrSets-1)*0.005+0.005), c(1,sizeTotal+nrAtoms), type="n", xlab="", ylab=ylab, axes="FALSE")
  if(is.null(setLab)) setLab <- setNames
  axis(1, at=0.51*(0:(nrSets-1))+0.25, labels=setLab, lty="blank")
  if(is.null(atomAvLab)) atomAvLab=rep("Average in\n atom", nrSets)
  if(!is.null(values)) axis(3, at=0.51*(0:(nrSets-1))+0.25, labels=atomAvLab, lty="blank")
  
  for(s in 1:nrSets)
  {
    ##print(s)
    set <- setNames[s]
    ##get atoms that are in this set
    atomsInSet <- sapply(atoms, function(a,s){all(is.element(a,s))}, sets[[set]])
    
    ##make all rectangles bgcol, except for those of the atoms belonging to the set
    polyCol <- rep(bgCol, nrAtoms)
    polyCol[atomsInSet] <- col
    
    for(l in 1:length(atoms))
    {
      polyAtoml <- cbind(c(0.51*(s-1),0.51*(s-1),0.51*(s-1)+0.5,0.51*(s-1)+0.5),  
                         c(cumSizesAtoms[l],cumSizesAtoms[l+1],
                           cumSizesAtoms[l+1],cumSizesAtoms[l])+l)
      if(is.null(values))
      {
        polygon(polyAtoml, col=polyCol[l], border=NA)
      } else {
        polyAtomlVal <- cbind(c(0.51*(s-1),0.51*(s-1),0.51*(s-1)+0.5*avAtoms[l],0.51*(s-1)+0.5*avAtoms[l]),  
                              c(cumSizesAtoms[l],cumSizesAtoms[l+1],
                                cumSizesAtoms[l+1],cumSizesAtoms[l])+l)
        polygon(polyAtoml, col=bgCol, border=NA)
        polygon(polyAtomlVal, col=polyCol[l], border=NA)		
      }
    }
  }
}

