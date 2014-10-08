\name{getInitialData}
\alias{getInitialData}
\alias{getInitialData,caClass-method}
\title{
Retrieve Initial data from caClass objects
}
\description{
Generic function to retrieve the initial data used in a components analysis \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}} object
}
\usage{
getInitialData(x, block=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified all information about the preprocessing of the data, as well as the data after preprocessing is applied, is returned}
  \item{block}{Character indicating the block of data that has to be returned. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed}
}
\value{
It returns a list with the information indicated by the parameters. By default a list with the following components:
\item{process}{methods applied to preprocess the data}
\item{preproData}{Preprocessed data}
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
getInitialData(res)
getInitialData(res, block="expr")
}
\keyword{multivariate}