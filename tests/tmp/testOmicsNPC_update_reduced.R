
#remember: different results between NPC and PC depending how much overlapping are the samples

#### Script for testing the omicsNPC function ####

#### Set up ####

#memory, library and sourcing 
rm(list = ls())
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsNPC.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')
library(foreach)
library(data.table)

#control panel
set.seed(12345)
combMethods <- c('Fisher', 'Liptak', 'Tippett');
numPerms <- 10000;
numCores <- 3;

#### Binary outcome, two datasets ####

#creating the data matrices
dataset1 <- matrix(rnorm(5*20), 5, 20)
rownames(dataset1) <- paste('g', 1:5, sep = '_')
colnames(dataset1) <- paste('sample', 1:20, sep = '_')
dataset2 <- matrix(rnorm(4*20), 4, 20)
rownames(dataset2) <- paste('p', 1:4, sep = '_')
colnames(dataset2) <- paste('sample', 21:40, sep = '_')

#creating the phenotype data
phenotypeData <- matrix(rnorm(40*3), 40, 3)
rownames(phenotypeData) <- paste('sample', 1:40, sep = '_')
colnames(phenotypeData) <- c('cov1', 'cov2', 'outcome')
phenotypeData[, 'outcome'] <- c(rep(1, 10), rep(0, 10),
                                rep(1, 10), rep(0, 10))

#one data matrix is actually differentially expressed
dataset1[ , 11:20] <- dataset1[ , 11:20] + 2 * matrix(runif(5 * 5), 5, 10)

#data mapping
dataMapping <- expand.grid(rownames(dataset1), rownames(dataset2))
dataMapping <- dataMapping[sample(1:nrow(dataMapping), 2), ]
colnames(dataMapping) <- c('dataset1', 'dataset2')

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                    phenotypeData = phenotypeData, numPerms = numPerms, numCores = numCores)
print(results$pvaluesNPC)

(x2 <- -2 * sum(log(c(results$dataset1$P.Value[dataMapping$dataset1[1]], 
                      results$dataset2$P.Value[dataMapping$dataset2[1]]))))
print(pchisq(x2, 4, lower.tail = FALSE))

(x2 <- -2 * sum(log(c(results$dataset1$P.Value[dataMapping$dataset1[2]], 
                      results$dataset2$P.Value[dataMapping$dataset2[2]]))))
print(pchisq(x2, 4, lower.tail = FALSE))

# 
# dataInput = list(dataset1 = dataset1, dataset2 = dataset2)
# dataTypes = rep('continuous', length(dataInput))
# outcomeName = NULL
# combMethods = c("Fisher", "Liptak", "Tippett")
# allCombinations = FALSE
# verbose = FALSE
# functionGeneratingIndex = NULL
# dataWeights = rep(1, length(dataInput))/length(dataInput)
# returnPermPvalues = FALSE

