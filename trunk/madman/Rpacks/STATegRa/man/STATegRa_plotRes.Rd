\name{plotRes}
\alias{plotRes}
\alias{plotRes,caClass,numeric,character,character,logical,character-method}
\title{
Plot results of Components Analysis
}
\description{
Function that allows to plot scatterplots for scores or loadings, for common and distinctive part as well as combined plots for plotting both parts together. Different combinations of the parameters produce different plots.
}
\usage{
plotRes(object, comps=c(1,2), what, type, combined, block, color=NULL,
        shape=NULL, labels=NULL, background=TRUE, palette=NULL, pointSize=4,
        labelSize=NULL, axisSize=NULL, titleSize=NULL)
}
\arguments{
  \item{object}{\code{caClass} object containing results of Components Analysis}
  \item{comps}{Components to plot. If combined=FALSE, it indicates the x and y components of the type and block chosen. If \code{combined=TRUE}, it indicates the component to plot for the first block of information and the component for the second block of information to plot together. By default the components are set to c(1,2) if \code{combined=FALSE} and to c(1,1) if \code{combined=TRUE}.}
  \item{what}{what=c("scores","loadings"). Are scores or loadings had to be represented?}
  \item{type}{Indicate which components should be represented, "common", "individual" or "both" in the case of combined plots}
  \item{combined}{Logical indicating if the plot is a simple plot representing two components from the same block of information, or combined representing in the same plot two component each of them from a different block of information.}
  \item{block}{Indicates which block has to be represented. It can be specify by a numeric value (1 or 2) or a character (name of the block in the input data provided to the \code{omicCompAnalysis()} analysis)}
  \item{color}{Character indicating with column of the pData in the initial ExpressionSet has to be used to set color in the graph. By default is the first column of the pData}
  \item{shape}{Character indicating with column of the pData in the initial ExpressionSet has to be used to set shape in the graph.}
  \item{labels}{Character indicating with column of the pData in the initial ExpressionSet has to be used to set labels in the graph.}
  \item{background}{Logical indicating if the plot have to be represented with grey background or not}
  \item{palette}{Vector indicating the pallete of colors used for the plot. See \code{\link{ggplot}}}
  \item{pointSize}{Size of the points represented in the plot}
  \item{labelSize}{Size of the labels if label is not NULL}
  \item{axisSize}{Size of the text in the axis of the plot}
  \item{titleSize}{Size of main title}
}
\value{
Plot of class \code{\link{ggplot}} 
}
\author{
Patricia Sebastian Leon
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

# Scatterplot of scores variables associated to common components

# DISCO-SCA
plotRes(object=discoRes,comps=c(1,2),what="scores",type="common",
        combined=FALSE,block="",color="classname",shape=NULL,labels=NULL,
        background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
        axisSize=NULL,titleSize=NULL)
# JIVE
plotRes(object=jiveRes,comps=c(1,2),what="scores",type="common",
        combined=FALSE,block="",color="classname",shape=NULL,labels=NULL,
        background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
        axisSize=NULL,titleSize=NULL)
        
# O2PLS
# Scatterplot of scores variables associated to common components
# Associated to first block
p1 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
        combined=FALSE,block="expr",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)
# Associated to second block
p2 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
        combined=FALSE,block="mirna",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)
        
# Combined plot of scores variables assocaited to common components
plotRes(object=o2plsRes,comps=c(1,1),what="scores",type="common",
        combined=TRUE,block="",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)

# Loadings plot for individual components
# Separately for each block
p1 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
        combined=FALSE,block="expr",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
        combined=FALSE,block="mirna",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)

# Combined plot
plotRes(object=discoRes,comps=c(1,1),what="loadings",type="individual",
        combined=TRUE,block="",color="classname",shape=NULL,
        labels=NULL,background=TRUE,palette=NULL,pointSize=4,
        labelSize=NULL,axisSize=NULL,titleSize=NULL)
}