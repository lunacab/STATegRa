\name{getPreprocessing}
\alias{getPreprocessing}
\alias{getPreprocessing,caClass-method}
\title{
Retrieve information about the preprocessing done in \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\description{
Generic function to retrieve information about preprocessing done in \code{\link{omicsCompAnalysis}} from an \code{\link{caClass-class}}
}
\usage{
getPreprocessing(x, process=FALSE, preproData=FALSE, block=NULL)
}
\arguments{
  \item{x}{\code{\link{caClass-class}} object. If only this parameter is specified loadings for common and distinctive part are diplayed}
  \item{process}{Logical indicating whether to return information about the processing done. Defaults to FALSE.}
  \item{preproData}{Logical indicating whether to return the pre-processed data matrices. Defaults to FALSE.}
  \item{block}{Character indicating the block of data whose pre-processed state has to be returned, if \code{preproData=TRUE}. It can be specified by using the position of the block ("1" or "2") or the name assigned to this block in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed}
}
\value{
If \code{process} and \code{preproData} are both requested, a list containing:
\item{process}{methods applied to preprocess the data}
\item{preproData}{Preprocessed data (either a single block or a list of blocks)}
Otherwise, either a character vector for \code{process} or either a matrix or list of matrices for \code{preproData}, depending on \code{block}.
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
getPreprocessing(res, process=TRUE)
getPreprocessing(res, preproData=TRUE, block="1")
}
\keyword{multivariate}
