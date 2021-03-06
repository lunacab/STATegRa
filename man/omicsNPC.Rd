% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STATegRa_omicsNPC.R
\name{omicsNPC}
\alias{omicsNPC}
\alias{omicsNPC,list,data.frame-method}
\alias{omicsNPC,list,missing-method}
\title{omicsNPC, applying the Non-Parametric Combination (NPC) on omics datasets}
\usage{
omicsNPC(dataInput, dataMapping, dataTypes = rep('continuous', length(dataInput)), 
               combMethods = c("Fisher", "Liptak", "Tippett"), numPerms = 1000, 
               numCores = 1, verbose = FALSE, functionGeneratingIndex = NULL, 
               outcomeName = NULL, allCombinations = FALSE, 
               dataWeights = rep(1, length(dataInput))/length(dataInput), 
               returnPermPvalues = FALSE, ...)
}
\arguments{
\item{dataInput}{List of ExpressionSet objects, one for each data modality.}

\item{dataMapping}{A data frame describing how to map measurements across datasets. See details for more information.}

\item{dataTypes}{Character vector with possible values: 'continuous', 'count'. Alternatively, a list of functions for assessing deregulation / association with an outcome}

\item{combMethods}{Character vector with possible values: 'Fisher', 'Liptak', 'Tippett'. If more than one is specified, all will be used.}

\item{numPerms}{Number of permutations}

\item{numCores}{Number of CPU cores to use}

\item{verbose}{Logical, if set to TRUE omicsNPC prints out the step that it performs}

\item{functionGeneratingIndex}{Function generating the indices for randomly permuting the samples}

\item{outcomeName}{Name of the outcome of interest / experimental factor, as reported in the design matrices. If NULL, the last column of the design matrices is assumed to be the outcome of interest.}

\item{allCombinations}{Logical, if TRUE all combinations of omics datasets are considered}

\item{dataWeights}{A vector specifying the weigth to give to each dataset. Note that sum(dataWeights) should be 1.}

\item{returnPermPvalues}{Logical, should the p-values computed at each permutation being returned?}

\item{...}{Additional arguments to be passed to the user-defined functions}
}
\value{
A list containing:
stats0 Partial deregulation / association statistics 
pvalues0 The partial p-values computed on each dataset
pvaluesNPC The p-values computed through NPC.
permPvalues The p-values computed at each permutation
}
\description{
This function applies the NonParametric Combination methodology on the integrative analysis 
of different omics data modalities.
It retrieves genes associated to a given outcome, taking into account all omics data.
First, each datatype is analyzed independently using the appropriate method.
omicsNPC analyses continuous data (microarray) using limma, while count data (e.g., RNAseq) 
are first preprocessed with using the "voom" function. The user can also specify their own 
function for computing deregulation / association 
The p-values from the single dataset analysis are then combined employing Fisher, 
Liptak and Tippett combining functions.
The Tippett function returns findings which are supported by at least one omics modality.
The Liptak function returns findings which are supportd by most modalities.
The Fisher function has an intermediate behavior between those of Tippett and Liptak.
}
\examples{
# Load the data
data("TCGA_BRCA_Batch_93")
# Setting dataTypes, the first two ExpressionSets include RNAseq data,
# the third ExpressionSet includes Microarray data.
dataTypes <- c("count", "count", "continuous")
# Setting methods to combine pvalues
combMethods = c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms = 1000
# Setting number of cores
numCores = 1
# Setting omicsNPC to print out the steps that it performs.
verbose = TRUE
# Run omicsNPC analysis.
# The output contains a data.frame of p-values, where each row corresponds to a gene, 
# and each column corresponds to a method used in the analysis.

\dontrun{out <- omicsNPC(dataInput = TCGA_BRCA_Data, dataTypes = dataTypes,
                            combMethods = combMethods, numPerms = numPerms,
                            numCores = numCores, verbose = verbose)}
}
\references{
Pesarin, Fortunato, and Luigi Salmaso. Permutation tests for complex data: theory, applications and software.
John Wiley & Sons, 2010.
Nestoras Karathanasis, Ioannis Tsamardinos and Vincenzo Lagani. omicsNPC: applying the Non-Parametric Combination 
methodology to the integrative analysis of heterogeneous omics data. PlosONE 11(11): e0165545. doi:10.1371/journal.pone.0165545
}
\author{
Nestoras Karathanasis, Vincenzo Lagani
}
