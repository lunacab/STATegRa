\name{modelSelection}
\alias{modelSelection}
\alias{modelSelection,list,numeric,character-method}
\title{
Number of optimal common and distinctive components in object-wise data
}
\description{
This function combine the function \code{\link{selectCommonComps}} and \code{\link{PCA.selection}} to stimate the optimal number of common and distinctive components in object-wise data.  That is, it calculates the optimal number of common and individual components depending on the maximal number of common components and the individual components selection criteria provided by the user. The result is a list with the optimal number of common and distinctive components. This optimal number can be use for component analysis with function \code{\link{omicsCompAnalysis}}
}
\usage{
modelSelection(Input,Rmax,fac.sel,varthreshold=NULL,nvar=NULL,PCnum=NULL)
}
\arguments{
  \item{Input}{List of \code{ExpressionSet} objects, one for each block of data}
  \item{Rmax}{Maximum number of common components}
  \item{fac.sel}{Criterium for selecting number of components. The posible option are: "\%accum", "single\%", "rel.abs" and "fixed.num"}
  \item{varthreshold}{Threshold for the selection of components in "\%accum", "single\%" criterions}
  \item{nvar}{Threshold applied when the option "rel.abs" is selected}
  \item{PCnum}{Fixed number of components to select when the option "fixed.num" is selected}
}
\value{
The function returns a list with the following components:
\item{common}{Number of optimal common components}
\item{dist}{Number of optimal distictive components for each block}
}
\author{
Patricia Sebastian-Leon
}
\seealso{
\code{\link{selectCommonComps}},\code{\link{PCA.selection}},\code{\link{omicsCompAnalysis}}
}
\examples{
data(STATegRa_S3)

## Create ExpressionSets
# Block1 - Expression data
B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,pDataDescr=c("classname"))
# Block2 - miRNA expression data
B2 <- createOmicsExpressionSet(Data=Block2.PCA,pData=ed.PCA,pDataDescr=c("classname"))

## Model Selection
ms <- modelSelection(Input=list(B1,B2),Rmax=4,fac.sel="single\%",
                      varthreshold=0.03)
ms
}
\keyword{multivariate}

