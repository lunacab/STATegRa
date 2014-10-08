\name{bioDistFeaturePlot}
\alias{bioDistFeaturePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bioDistFeaturePlot
}
\description{
Function that pltos the results from a bioDistFeature analysis.
}
\usage{
bioDistFeaturePlot(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix, the output of a bioDistFeature analysis}
}
\value{
Generates a plot with the computed summary.
}
\author{
David Gomez-Cabrero
}
\seealso{
\code{\link{bioDistFeature}},\code{\link{bioDistW}}
}
\examples{
# Definition of a bioDistclass Object
data(STATegRa_S1)
data(STATegRa_S2)

# Truncate data for brevity
Block1 <- Block1[1:100,]
Block2 <- Block2[1:100,]

## Create ExpressionSets
mRNA.ds <- createOmicsExpressionSet(Data=Block1,pData=ed,pDataDescr=c("classname"))
miRNA.ds <- createOmicsExpressionSet(Data=Block2,pData=ed,pDataDescr=c("classname"))
  
## Create the bioMap  
map.gene.miRNA<-bioMap(name = "Symbol-miRNA",
                metadata =  list(type_v1="Gene",type_v2="miRNA",
                                 source_database="targetscan.Hs.eg.db",
                                 data_extraction="July2014"),
                map=mapdata)  

require(Biobase)

# Create Gene-gene distance computed through miRNA data
bioDistmiRNA<-bioDist(referenceFeatures = rownames(Block1),     
             reference = "Var1",
             mapping = map.gene.miRNA,
             surrogateData = miRNA.ds,  ### miRNA data
             referenceData = mRNA.ds,  ### mRNA data
             maxitems=2,
             selectionRule="sd",
             expfac=NULL,
             aggregation = "sum",
             distance = "spearman",
             noMappingDist = 0,
             filtering = NULL,
             name = "mRNAbymiRNA")

# Create Gene-gene distance through mRNA data
bioDistmRNA<-new("bioDistclass",
                 name = "mRNAbymRNA",
                 distance = cor(t(exprs(mRNA.ds)),method="spearman"),
                 map.name = "id",
                 map.metadata = list(),
                 params = list())

###### Generation of the list of Surrogated distances.

bioDistList<-list(bioDistmRNA,bioDistmiRNA)
sample.weights<-matrix(0,4,2)
sample.weights[,1]<-c(0,0.33,0.67,1)
sample.weights[,2]<-c(1,0.67,0.33,0)

###### Generation of the list of bioDistWclass objects.

bioDistWList<-bioDistW(referenceFeatures = rownames(Block1),
                       bioDistList = bioDistList,
                       weights=sample.weights)
                       
###### Computing the matrix of features/distances associated.

fm<-bioDistFeature(Feature = rownames(Block1)[1] ,
                       listDistW = bioDistWList,
                       threshold.cor=0.7)
bioDistFeaturePlot(data=fm)
                       
}
