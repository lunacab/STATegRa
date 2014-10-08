\name{getVAF}
\alias{getVAF}
\alias{getVAF,caClass-method}
\title{
Retrieve information about the VAF obtained using \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\description{
Generic function to retrieve loadings (common and distinctive) obtained in the \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\usage{
getVAF(x, part=NULL, block=NULL)
}
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified loadings for common and distinctive part are diplayed}
  \item{part}{Character indicating if VAF for "common" or "distinctive" part has to be displayed. If the method is DISCO-SCA, another option, "cross", indicating the associated error due to the rotation is allowed}
  \item{block}{Character indicating the block of data whose VAF has to be returned. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both VAFs are displayed}
}
\details{
It is important to note that in the case of O2PLS method, it is not possible to calculate the VAF because the components are not orthogonal
}
\value{
By default, a list with VAF associated to both parts (common and distinctive). If any parameter is specified, the selected information is displayed.
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
getVAF(res)
getVAF(res, part="common")
getVAF(res, part="distinctive", block="expr")
}
\keyword{multivariate}
