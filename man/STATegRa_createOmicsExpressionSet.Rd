\name{createOmicsExpressionSet}
\alias{createOmicsExpressionSet}
\alias{createOmicsExpressionSet,matrix-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an ExpressionSet object since omics dataset.
}
\description{
This function allow to the user to create a ExpressionSet object from a matrix representing an omics dataset. It allows to include the experimental design and annotation in the ExpressionSet object.
}
\usage{
createOmicsExpressionSet(Data,pData=NULL,pDataDescr=NULL,feaData=NULL,feaDataDescr=NULL)
}
\arguments{
  \item{Data}{Matrix with omics dataset. }
  \item{pData}{Matrix indicating the phenotype of the samples (columns). rownames(pData)==colnames(exprs)}
  \item{pDataDescr}{Description of the phenotypic data}
  \item{feaData}{Matrix indicating the variables/features annotation}
  \item{feaDataDescr}{Description of the feature annotation}
  }
\details{
  In Data matrix, samples has to be in columns and variables has to be in rows
}
\value{
A object of class \code{\link{ExpressionSet}}  
}
\author{
Patricia Sebastian-Leon
}
\examples{
data(STATegRa_S3)
B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,
                               pDataDescr=c("classname"))
B2 <- createOmicsExpressionSet(Data=Block2.PCA,pData=ed.PCA,
                               pDataDescr=c("classname"))
}
\keyword{datagen}

