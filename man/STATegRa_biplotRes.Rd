\name{biplotRes}
\alias{biplotRes}
\alias{biplotRes,caClass,character,numeric,character-method}
\title{
Biplot results of Components Analysis
}
\description{
Function to display a biplot of an \code{caClass} object coming from a Component Analysis.
}
\usage{
biplotRes(object,type,comps,block,title=NULL,colorCol=NULL,sizeValues=c(2,4),
          shapeValues=c(17,0),background=TRUE,pointSize=4,labelSize=NULL,
          axisSize=NULL,titleSize=NULL)
}
\arguments{
  \item{object}{\code{caClass} object containing results of Components Analysis}
  \item{type}{Indicate which components should be represented, "common", "individual" or "both" in the case of combined plots}
  \item{comps}{Components to plot. If combined=FALSE, it indicates the x and y components of the type and block chosen. If \code{combined=TRUE}, it indicates the component to plot for the first block of information and the component for the second block of information to plot together. By default the components are set to c(1,2) if \code{combined=FALSE} and to c(1,1) if \code{combined=TRUE}.}
  \item{block}{Indicates which block has to be represented. It can be specify by a numeric value (1 or 2) or a character (name of the block in the input data provided to the \code{omicCompAnalysis()} analysis)}
  \item{title}{Title of the biplot}
  \item{colorCol}{Character indicating with column of the pData in the initial ExpressionSet has to be used to set color in the graph. By default is the first column of the pData}
  \item{sizeValues}{Vector indicating which sizes has to be use for scores and loadings?}
  \item{shapeValues}{Vector indicating which shapes has to be used for scores and loadings}
  \item{background}{Logical indicating if the plot have to be represented with grey background or not}
  \item{pointSize}{Size of the points represented in the plot}
  \item{labelSize}{Size of the labels if label is not NULL}
  \item{axisSize}{Size of the text in the axis of the plot}
  \item{titleSize}{Size of main title}
}
\value{
Plot of class \code{\link{ggplot}} 
}
\author{
Patricia Sebastian-Leon
}
\seealso{
\code{\link{ggplot}},\code{\link{omicsCompAnalysis}}
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

# Biplot common part. DISCO-SCA

biplotRes(object=discoRes,type="common",comps=c(1,2),block="",
          title=NULL,colorCol="classname",sizeValues=c(2,4),
          shapeValues=c(17,0),background=TRUE,pointSize=4,
          labelSize=NULL,axisSize=NULL,titleSize=NULL)

# Biplot common part. O2PLS

p1 <- biplotRes(object=o2plsRes,type="common",comps=c(1,2),
          block="expr",title=NULL,colorCol="classname",
          sizeValues=c(2,4),shapeValues=c(17,0),
          background=TRUE,pointSize=4,labelSize=NULL,
          axisSize=NULL,titleSize=NULL)
p2 <- biplotRes(object=o2plsRes,type="common",comps=c(1,2),
          block="mirna",title=NULL,colorCol="classname",
          sizeValues=c(2,4),shapeValues=c(17,0),
          background=TRUE,pointSize=4,labelSize=NULL,
          axisSize=NULL,titleSize=NULL)

# Biplot distinctive part. O2PLS

p1 <- biplotRes(object=discoRes,type="individual",comps=c(1,2),
          block="expr",title=NULL,colorCol="classname",
          sizeValues=c(2,4),shapeValues=c(17,0),
          background=TRUE,pointSize=4,labelSize=NULL,
          axisSize=NULL,titleSize=NULL)
p2 <- biplotRes(object=discoRes,type="individual",comps=c(1,2),
          block="mirna",title=NULL,colorCol="classname",
          sizeValues=c(2,4),shapeValues=c(17,0),
          background=TRUE,pointSize=4,labelSize=NULL,
          axisSize=NULL,titleSize=NULL)
}
