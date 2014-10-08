\name{getMethodInfo}
\alias{getMethodInfo}
\alias{getMethodInfo,caClass-method}
\title{
Retrieve information about the components analysis done using \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\description{
Generic function to retrieve information about the method and the number of components (common and individual) employed in the \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\usage{
getMethodInfo(x, method=FALSE, comps=NULL, block=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified all information about the method and number of components is diplayed}
  \item{method}{Logical indicating if method employed has to be returned (FALSE by default)}
  \item{comps}{Character indicating which components number has to be returned (The options are: "common" for common components, "distinctive" for distinctive components}
  \item{block}{Character indicating the block of data whose number of components has to be returned. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both number of distinctive components are displayed}
}
\value{
It returns a list with the information indicated by the parameters. By default a list with the following components:
\item{method}{The component analysis method used for the analysis}
\item{comps}{Number of components (common and distinctive) used as input for the method}
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
getMethodInfo(res)
getMethodInfo(res, method=TRUE)
getMethodInfo(res, comps="all", block="expr")
}
\keyword{multivariate}
