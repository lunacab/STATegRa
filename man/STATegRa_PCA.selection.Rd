\name{PCA.selection}
\alias{PCA.selection}
\alias{PCA.selection,matrix,character-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of optimal number of components using PCA
}
\description{
Function to select the number of optimal components using Principal Components Analysis. There are four different criterions to select the number of components: accumulated variance explained, the individual explained variance of each component, the absolute value of its 
variability or just a fixed number of components.
}
\usage{
PCA.selection(Data,fac.sel,varthreshold=NULL,nvar=NULL,PCnum=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{Matrix with data. Samples has to be in columns and variables in rows}
  \item{fac.sel}{Criterium for selecting number of components. The posible option are: "\%accum", "single\%", "rel.abs" and "fixed.num"}
  \item{varthreshold}{Threshold for the selection of components in "\%accum", "single\%" criterions}
  \item{nvar}{Threshold applied when the option "rel.abs" is selected}
  \item{PCnum}{Fixed number of components to select when the option"fixed.num" is selected}
}
\value{
The function returns a list with the following components:
\item{PCAres}{A list containing the results of the PCA decomposition}
\item{numComps}{Number of selected components applying the selected criterion}

The components of \code{PCAres} list are the following:
\item{eigen}{Eigen values of the decomposition}
\item{var.exp}{Explained variance of each component}
\item{scores}{Matrix of scores}
\item{loadings}{Matrix of loadings}

}
\author{
Patricia Sebastian Leon
}
\seealso{
\code{\link{modelSelection}}
}
\examples{
data(STATegRa_S3)
ps <- PCA.selection(Data=Block2.PCA,fac.sel="single\%",varthreshold=0.03)
ps$numComps
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

