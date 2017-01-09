# This file contains no code, it just serves to hold roxygen documentation comments
# for the package and data files so that we can do away with all Rd files

#' @title STATegRa
#' @name STATegRa
#' @aliases STATegRa-package
#' @docType package
#' @description
#' STATegRa is a package for the integrative analysis of multi-omic data-sets.
#'
#' For full information, see the user's guide.
#' @seealso \code{\link{STATegRaUsersGuide}}
NULL

#' @title STATegRa data
#' @name STATegRa_data
#' @aliases Block1 Block2 ed mapdata Block1.PCA Block2.PCA ed.PCA
#' @docType data
#' @format Two matrices with mRNA and miRNA expression data, a design matrix that describes both and a mapping between miRNA and genes.
#' @description
#' mRNA data (\code{Block1}), miRNA data (\code{Block2}) and the design matrix (\code{ed}), from \code{STATegRa_S1}, provides selected data downloaded from \url{https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/}. The mapping between miRNA and mRNA (\code{mapdata}, available in \code{STATegRa_S2}) contains, as a processed matrix, selected information available from TargetScan; we selected the set of miRNA target predictions for humans for those miRNA-mRNA pairs where both miRNA and mRNA were in \code{Block1} and \code{Block2} respectively. 
#'
#' The PCA version of the data (\code{Block1.PCA}, \code{Block2.PCA}, \code{ed.PCA}; available in \code{STATegRa_S3}), provides a similar data-set to \code{Block1}, \code{Block2} and \code{ed} data; however in this case the data has been processed in order to provide a pedagogic example of OmicsPCA. Results obtained from OmicsPCA (\code{\link{omicsCompAnalysis}}) with the existing data should not be taken as clinically valid.
#'
#' @source
#' (a)  See \url{https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/}.
#' (b)  Gabor Csardi, \code{targetscan.Hs.eg.db}: TargetScan miRNA target predictions for human. R package version 0.6.1
#'
#' @author David Gomez-Cabrero, Patricia Sebastian-Leon, Gordon Ball
#' @examples
#' data(STATegRa_S1)
#' data(STATegRa_S2)
#' data(STATegRa_S3)
NULL

#' @title STATegRa data
#' @name STATegRa_data_TCGA_BRCA
#' @aliases TCGA_BRCA_Data
#' @docType data
#' @format One list, which contains three ExpressionSet objects.
#' @description
#'
#' Data were downloaded from TCGA data portal, \url{https://tcga-data.nci.nih.gov/tcga/}.
#' We downloaded sixteen tumour samples and the sixteen matching normal, for Breast invasive carcinoma, BRCA, batch 93.
#' Herein, three types of data modalities are included, RNAseq (\code{TCGA_BRCA_Data$RNAseq}), RNAseqV2 (\code{TCGA_BRCA_Data$RNAseqV2})
#' and Expression-Genes (\code{TCGA_BRCA_Data$Microarray}). The Data Level was set to Level 3.
#' For each data type, we pooled all data to one matrix, where rows corresponded to genes and columns to samples.
#' Only the first 100 genes are included.
#'
#' @source
#' See \url{https://tcga-data.nci.nih.gov/tcga/}.
#'
#' @author Nestoras Karathanasis, Vincenzo Lagani
#' @examples
#' # load data
#' data(TCGA_BRCA_Batch_93)
NULL

#' @name STATegRa-defunct
#' @title Defunct functions in STATegRa
#' @description
#' These functions have are defunct and no longer available
#' @details
#' \itemize{
#' \item{holistOmics: replaced by \code{\link{omicsNPC}}}
#' }
NULL