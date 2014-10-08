\docType{data}
\name{STATegRa_initial_data}
\alias{STATegRa_initial_data}
\alias{Block1}
\alias{Block2}
\alias{ed}
\alias{mapdata}
\alias{Block1.PCA}
\alias{Block2.PCA}
\alias{ed.PCA}

\title{STATegRa_initial_data}
\format{Two matrices with mRNA and miRNA expression data, a design matrix that describes both and a mapping between miRNA and genes.}
\description{
mRNA data (\code{Block1}), miRNA data (\code{Block2}) and the design matrix (\code{ed}), from \code{STATegRa_S1}, provides selected data downloaded from \url{https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/}. The mapping between miRNA and mRNA (\code{mapdata}, available in \code{STATegRa_S2}) contains, as a processed matrix, selected information available from TargetScan; we selected the set of miRNA target predictions for humans for those miRNA-mRNA pairs where both miRNA and mRNA were in \code{Block1} and \code{Block2} respectively. 

The PCA version of the data (\code{Block1.PCA}, \code{Block2.PCA}, \code{ed.PCA}; available in \code{STATegRa_S3}), provides a similar data-set to \code{Block1}, \code{Block2} and \code{ed} data; howere in thiscase the data has been processed in order to provide a pedagogic example of OmicsPCA. Results obtained from OmicsPCA (\code{\link{omicsCompAnalysis}}) with the existing data should not be taken as clinically valid. 
}
\source{
     (a)  See \url{https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/}.
     (b)  Gabor Csardi, \code{targetscan.Hs.eg.db}: TargetScan miRNA target predictions for human. R package version 0.6.1.
}
\author{
      David Gomez-Cabrero, Patricia Sebastian-Leon, Gordon Ball 2014-10-03
}
\usage{STATegRa_S1}
\keyword{datasets}