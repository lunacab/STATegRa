\name{caClass-class}
\docType{class}
\alias{class:caClass}
\alias{caClass-class}
\title{
caClass
}
\description{
Stores the results of any of the omicsPCA analyses.
}
\author{
Patricia Sebastian Leon
}
\seealso{
\code{\link{omicsCompAnalysis}}
Accessor methods: \code{\link{getInitialData}}, \code{\link{getMethodInfo}}, \code{\link{getScores}}, \code{\link{getLoadings}}, \code{\link{getVAF}}, \code{\link{getPreprocessing}}
}
\examples{
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
# all caClass instances
class(discoRes)
class(jiveRes)
class(o2plsRes)

}
