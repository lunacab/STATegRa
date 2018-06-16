
#### Script for testing the omicsNPC function ####

#### Set up ####

#memory, library and sourcing 
rm(list = ls())
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsNPC.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')
source('../../R/STATegRa_omicsPC.R')
source('../../R/STATegRa_omicsPC_internal.R')
library(foreach)
library(data.table)

#control panel
set.seed(12345)
combMethods <- c('Fisher', 'Liptak', 'Tippett');
numPerms <- 10;
numCores <- 3;

#### Binary outcome, two datasets ####

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
phenotypeData[, 'outcome'] <- c(rep(1, 15), rep(0, 15))

#one data matrix is actually differentially expressed
dataset1[ , 16:20] <- dataset1[ , 16:20] + 3 * matrix(runif(100 * 5), 100, 5)

#data mapping
dataMapping <- expand.grid(rownames(dataset1), rownames(dataset2))
dataMapping <- dataMapping[sample(1:nrow(dataMapping), 100), ]
colnames(dataMapping) <- c('dataset1', 'dataset2')

#omicsNPC
resultsNPC <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                  phenotypeData = phenotypeData, numPerms = numPerms, numCores = numCores)
print(resultsNPC)

#omicsPC
resultsPC <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), 
                   dataMapping = dataMapping, phenotypeData = phenotypeData)
print(resultsPC)

#computing fisher p-values from scratch
p1 <- resultsPC$dataset1[dataMapping$dataset1, 'P.Value']
p2 <- resultsPC$dataset2[dataMapping$dataset2, 'P.Value']
x2 <- -2 * (log(p1) + log(p2))
fisherPvalue <- pchisq(x2, 4, lower.tail = FALSE)

#comparing with PC and NPC
cor(resultsPC$pvaluesPC[, 'Fisher'], fisherPvalue)
cor(resultsNPC$pvaluesNPC[, 'Fisher'], fisherPvalue)

#computing FDR from scratch
fisherQvalue <- p.adjust(fisherPvalue, method = 'fdr')

#comparing with PC and NPC
cor(resultsPC$qvaluesPC[, 'Fisher'], fisherQvalue)
cor(resultsNPC$qvaluesNPC[, 'Fisher'], fisherQvalue)
plot(resultsNPC$qvaluesNPC[, 'Fisher'], fisherQvalue)

#### Multi-class outcome, two datasets ####

#changing phenotype
phenotypeData <- as.data.frame(phenotypeData)
phenotypeData[, 'outcome'] <- rep(letters[1:3], 10)

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                      phenotypeData = phenotypeData, numPerms = numPerms, numCores = numCores)
print(results)

#### Survival outcome, two datasets ####

#changing phenotype
phenotypeData[, 'outcome'] <- Surv(time = runif(30), 
                                   event = sample(0:1, size = 30, replace = TRUE))

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                    phenotypeData = phenotypeData, numPerms = numPerms, numCores = numCores)
print(results)

#### Binary outcome, three datasets ####

#changing phenotype
phenotypeData[, 'outcome'] <- c(rep(1, 15), rep(0, 15))

#additional dataset
dataset3 <- matrix(rnorm(100*30), 100, 30)
rownames(dataset3) <- paste('m', 1:100, sep = '_')
colnames(dataset3) <- paste('sample', 1:30, sep = '_')

#data mapping
dataMapping <- expand.grid(rownames(dataset1), rownames(dataset2), rownames(dataset3))
dataMapping <- dataMapping[sample(1:nrow(dataMapping), 100), ]
colnames(dataMapping) <- c('dataset1', 'dataset2', 'dataset3')

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2, dataset3=dataset3), 
                      dataMapping = dataMapping, phenotypeData = phenotypeData,  
                      allCombinations = FALSE, numPerms = numPerms, numCores = numCores)
print(results)

#### Binary outcome, three datasets, all combinations ####

#omicsNPC
results <- omicsNPC(dataInput = list(dataset1=dataset1, dataset2=dataset2, dataset3=dataset3), 
                    dataMapping = dataMapping, phenotypeData = phenotypeData,  
                    allCombinations = TRUE, numPerms = numPerms, numCores = numCores)
print(results)


