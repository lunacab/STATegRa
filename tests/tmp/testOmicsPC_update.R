
#### Script for testing the omicsPC function ####

#### Set up ####

#memory, library and sourcing 
rm(list = ls())
source('../../R/STATegRa_omicsPC_internal.R')
source('../../R/STATegRa_omicsPC.R')
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')

#control panel
set.seed(12345)

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

#data mapping
dataMapping <- expand.grid(rownames(dataset1), rownames(dataset2))
dataMapping <- dataMapping[sample(1:nrow(dataMapping), 100), ]
colnames(dataMapping) <- c('dataset1', 'dataset2')

#omicsPC
results <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                    phenotypeData = phenotypeData)
print(results)

#### Multi-class outcome, two datasets ####

#changing phenotype
phenotypeData <- as.data.frame(phenotypeData)
phenotypeData[, 'outcome'] <- rep(letters[1:3], 10)

#omicsNPC
results <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                      phenotypeData = phenotypeData)
print(results)

#### Survival outcome, two datasets ####

#changing phenotype
phenotypeData[, 'outcome'] <- Surv(time = runif(30), 
                                   event = sample(0:1, size = 30, replace = TRUE))

#omicsPC
results <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2), dataMapping = dataMapping, 
                    phenotypeData = phenotypeData)
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

#omicsPC
results <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2, dataset3=dataset3), 
                      dataMapping = dataMapping, phenotypeData = phenotypeData)
print(results)

#### Binary outcome, three datasets, all combinations ####

#omicsPC
results <- omicsPC(dataInput = list(dataset1=dataset1, dataset2=dataset2, dataset3=dataset3), 
                    dataMapping = dataMapping, phenotypeData = phenotypeData,  
                    allCombinations = TRUE)
print(results)
