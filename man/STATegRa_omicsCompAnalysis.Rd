\name{omicsCompAnalysis}
\alias{omicsCompAnalysis}
\alias{omicsCompAnalysis,list,character,character,numeric,numeric-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Components analysis for the analysis of object wise data
}
\description{
This function performs a components analysis of object wise omics data to understand the mechanisms that underlay all the data blocks under study (common mechanisms) and the mechanisms underlying each of the data block independently (distinctive mechanisms). This analysis include both, the preprocessing of data and the components analysis by using three different methodologies.
}
\usage{
omicsCompAnalysis(Input,Names,method,Rcommon,Rspecific,
                  convThres=10^(-10),maxIter=600,center=FALSE,
                  scale=FALSE,weight=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Input}{List of \code{ExpressionSet} objects, one for each block of data}
  \item{Names}{Vector indicating a distintive name for each block of data}
  \item{method}{Method used for components analysis. The methods included are DISCO_SCA, JIVE and O2PLS}
  \item{Rcommon}{Number of common components between all blocks}
  \item{Rspecific}{Vector indicating the individual components of each block of data}
  \item{convThres}{Maximum value for convergence (stop criterium. Default value: 0.0000001)}
  \item{maxIter}{Maximum number of iterations (stop criterium. Default value: 600)}
  \item{center}{Character indicating which type of centering must be applied to data before the analysis. Data can be centered considering all blocks together (PERBLOCKS) or each block separately (PERBLOCKS). If the parameter has a NULL (default) value, no center is applied}
  \item{scale}{Character indicating which type of scaling must be applied to data before the analysis. Data can be scaled considering all blocks together (PERBLOCKS) or each block separately (PERBLOCKS). If the parameter has a NULL (default) value, no scale is applied}
  \item{weight}{Logical indicating if the blocks has to be weighted before the analysis}
}
\value{
A object of class \code{\link{caClass-class}}
}
\author{
Patricia Sebastian Leon
}
\seealso{
\code{\link{caClass-class}}
}
\examples{
# Loading data
data("STATegRa_S3")
B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,
                               pDataDescr=c("classname"))
B2 <- createOmicsExpressionSet(Data=Block2.PCA,
                               pData=ed.PCA,pDataDescr=c("classname"))
# Omics components analysis
discoRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
                              method="DISCOSCA",Rcommon=2,Rspecific=c(2,2),
                              center=TRUE,scale=TRUE,weight=TRUE)
jiveRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
                              method="JIVE",Rcommon=2,Rspecific=c(2,2),
                              center=TRUE,scale=TRUE,weight=TRUE)
o2plsRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
                              method="O2PLS",Rcommon=2,Rspecific=c(2,2),
                              center=TRUE,scale=TRUE,weight=TRUE)
}
\keyword{multivariate}