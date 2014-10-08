\name{getLoadings}
\alias{getLoadings}
\alias{getLoadings,caClass-method}
\title{
Retrieve information about the loadings obtained using \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\description{
Generic function to retrieve loadings (common and distinctive) obtained in the \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\usage{
getLoadings(x, part=NULL, block=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified loadings for common and distinctive part are diplayed}
  \item{part}{Character indicating if loadings for "common" or "distinctive" part has to be displayed}
  \item{block}{Character indicating the block of data whose scores has to be returned. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both loadings are displayed}
}
\value{
By default, a list with loadings associated to both parts (common and distinctive). If any parameter is specified, the selected information is displayed.
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
getLoadings(res)
getLoadings(res, part="common", block="expr")
getLoadings(res, part="distinctive", block="expr")
}
\keyword{multivariate}
