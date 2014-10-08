\name{bioDistclass}
\alias{bioDistclass-constructor}
\alias{bioDistclass}
\title{
bioDistclass
}
\description{
Class to manage mappings between genomic features.
}
\usage{
bioDistclass(name,distance, map.name,map.metadata,params)
}
\arguments{
  \item{name}{name assigned to the object}
  \item{distance}{A matrix object that stores the distance between features}.                          
  \item{map.name}{Character that stores the name of the bioMap object used to compute the distance.}
  \item{map.metadata}{ a list of parameters used to generate the mapping}
  \item{params}{ a list of parameters used to generate the distance}
}
\author{
David Gomez-Cabrero
}
\value{
An object of class \code{bioDistclass}
}
\seealso{
\code{\link{bioMap}}
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
