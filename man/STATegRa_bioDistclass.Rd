\name{bioDistclass-class}
\docType{class}
\alias{class:bioDistclass}
\alias{bioDistclass-class}
\title{
bioDistclass
}
\description{
Class to manage and store surrogates distances between features.
}
\section{Constructor}{
  bioDistclass
}
\author{
David Gomez-Cabrero
}
\seealso{
\code{\link{bioDistclass-constructor}}
}
\examples{
# Definition of a bioDistclass Object
data(STATegRa_S1)
## Create ExpressionSets
mRNA.ds <- createOmicsExpressionSet(Data=Block1,pData=ed,pDataDescr=c("classname"))

require(Biobase)

## Create bioDistclass Object
bioDistmRNA<-bioDistclass(name = "mRNAbymRNA",
                 distance = cor(t(exprs(mRNA.ds)),method="spearman"),
                 map.name = "id",
                 map.metadata = list(),
                 params = list())
}
