
#set up
rm(list = ls())
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsNPC.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')
set.seed(12345)
library(foreach)
library(data.table)

#control panel
combMethods <- c('Fisher', 'Liptak', 'Tippett');
numPerms <- 10;
numCores <- 1;

#creating the data matrices
dataset1 <- matrix(rnorm(100*20), 100, 20)
rownames(dataset1) <- paste('g', 1:100, sep = '_')
colnames(dataset1) <- paste('sample', 1:20, sep = '_')
dataset2 <- matrix(rnorm(100*20), 100, 20)
rownames(dataset2) <- paste('p', 1:100, sep = '_')
colnames(dataset2) <- paste('sample', 11:30, sep = '_')

#creating the phenotype data
phenotypeData <- matrix(rnorm(30*3), 30, 3)
rownames(phenotypeData) <- paste('sample', 1:30, sep = '_')
colnames(phenotypeData) <- c('cov1', 'cov2', 'outcome')

#data mapping
dataMapping <- expand.grid(rownames(dataset1), rownames(dataset2))
dataMapping <- dataMapping[sample(1:nrow(dataMapping), 100), ]
colnames(dataMapping) <- c('dataset1', 'dataset2')

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                    phenotypeData = phenotypeData, numPerms = numPerms, numCores = 1)

print(results)

# dataInput = list(dataset1, dataset2)
# dataTypes = rep('continuous', length(dataInput))
# outcomeName = NULL
# combMethods = c("Fisher", "Liptak", "Tippett")
# allCombinations = FALSE
# verbose = FALSE
# functionGeneratingIndex = NULL
# dataWeights = rep(1, length(dataInput))/length(dataInput)
# returnPermPvalues = FALSE
