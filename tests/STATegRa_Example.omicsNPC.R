rm(list = ls())
require("STATegRa")
# Load the data
data("TCGA_BRCA_Batch_93")
# Setting dataTypes
dataTypes <- c("count", "count", "continuous")
# Setting methods to combine pvalues
combMethods = c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms = 1000
# Setting number of cores
numCores = 1
# Setting holistOmics to print out the steps that it performs.
verbose = TRUE
# Run holistOmics analysis.
output <- omicsNPC(dataInput = TCGA_BRCA_Data, dataTypes = dataTypes, combMethods = combMethods, numPerms = numPerms, numCores = numCores, verbose = verbose)
