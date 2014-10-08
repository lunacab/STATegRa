\name{getScores}
\alias{getScores}
\alias{getScores,caClass-method}
\title{
Retrieve information about the scores obtained using \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\description{
Generic function to retrieve scores (common and distinctive) obtained in the \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\usage{
getScores(x, part=NULL, block=NULL)
}
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified scores for common and distinctive part are diplayed}
  \item{part}{Character indicating if scores for "common" or "distinctive" part has to be displayed}
  \item{block}{Character indicating the block of data whose scores has to be returned. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both scores are displayed}
}
\value{
By default, a list with scores associated to both parts (common and distinctive). If any parameter is specified, the selected information is displayed.
}
\author{
Patricia Sebastian-Leon
}
\seealso{
\code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
}
\examples{
data("STATegRa_S3")
B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
                               pDataDescr=c("classname"))
B2 <- createOmicsExpressionSet(Data=Block2.PCA,
                               pData=ed.PCA, pDataDescr=c("classname"))
# Omics components analysis
res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
                         method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
                         center=TRUE, scale=TRUE, weight=TRUE)
getScores(res)
getScores(res, part="common")
getScores(res, part="distinctive", block="expr")
}
\keyword{multivariate}
