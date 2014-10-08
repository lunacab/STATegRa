\name{plotVAF}
\alias{plotVAF}
\alias{plotVAF,caClass-method}
\title{
Function to plot the VAF (Variance Explained For) obtained from Component Analysis.
}
\description{
This function allows to the user to visualize the VAF results coming from a component analysis. The input is an object resulting from \code{omicsCompAnalysis} function of class \code{caClass}. The variance explained for (VAF) each component is given in VAF slot, that is
a list containing the VAF for common and distinctive components. In the case of
O2PLS, VAF cannot be calculated, because components are not orthogonal. Structure of
plots associated to DISCO-SCA and JIVE are different. In the case of DISCO-SCA,
components of individual blocks have an associated error due to the rotation is not
perfect. This is the reason because the DISCO-SCA distinctive components have VAF
in the other block. This VAF not associated to the corresponding block could be
interpreted as the error for not having a perfect rotation.
}
\usage{
plotVAF(object,mainTitle)
}
\arguments{
  \item{object}{\code{caClass} object containing results of Components Analysis}
  \item{mainTitle}{Title of the plot}
}
\value{
Plot representing the VAF of common and distintive components of each block.
}
\author{
Patricia Sebastian -Leon
}
\seealso{
\code{\link{omicsCompAnalysis}}
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
require(ggplot2)

# DISCO-SCA plotVAF
plotVAF(discoRes)

# JIVE plotVAF
plotVAF(jiveRes)
}
\keyword{hplot}

