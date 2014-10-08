\name{selectCommonComps}
\alias{selectCommonComps}
\alias{selectCommonComps,matrix,matrix,numeric-method}
\title{
Select common components of two blocks of object wise omics data.
}
\description{
This function applies a Simultaneous Component Analysis (SCA). The idea is that the scores for both blocks
should have a similar behavior if the components are in the common mode. For evaluate if a component is common or not, the ratios be-
tween the explained variances (SSQ) of each block and its estimationare used. The highest component having its ratios between 0.8
and 1.5 is selected as the optical number of common components. 
}
\usage{
selectCommonComps(X,Y,Rmax)
}
\arguments{
  \item{X}{First block of omics data}
  \item{Y}{Second block of omics data}
  \item{Rmax}{Maximum number of common components}
}
\value{
The function returns a list with the following components:
\item{common}{Number of optimal common components}
\item{ssqs}{Matrix with SSQ for each block and its estimation}
\item{pssq}{\code{\link{ggplot}} object representing SSQ for each block and its estimation}
\item{pratios}{\code{\link{ggplot}} object representing the SSQ ratios between each block and its estimation}
}
\author{
Patricia Sebastian-Leon
}
\seealso{
\code{\link{modelSelection}},\code{\link{ggplot}}
}
\examples{
data(STATegRa_S3)

cc <- selectCommonComps(X=Block1.PCA,Y=Block2.PCA,Rmax=3)
cc$common
# Graphical output
cc$pssq
cc$pratios
}
\keyword{multivariate}

