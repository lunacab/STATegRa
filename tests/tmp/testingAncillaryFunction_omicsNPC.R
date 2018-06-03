
rm(list = ls());
set.seed(12345)
library(survival)
library(limma)
library(edgeR)
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')

#### normal ####

dataMatrix <- matrix(rnorm(1000), 100, 10);
design <- data.frame(matrix(rnorm(20), 10, 2))
outcomeName <- 'X1'
a <- computeAssocContinuousData(dataMatrix, design, outcomeName)
class(a)
a

dataMatrixCount <- matrix(rnbinom(n = 1000, size = 1000, prob = 0.2), 100, 10);
a <- computeAssocCountData(dataMatrixCount, design, outcomeName)
a

design <- cbind(design, target = factor(rep(c(1,2), 5)))
outcomeName <- 'target'
a <- computeAssocContinuousData(dataMatrix, design, outcomeName)
a
a <- computeAssocCountData(dataMatrixCount, design, outcomeName)
a

levels(design$target) <- 1:3;
design$target[7:10] <- 3
a <- computeAssocContinuousData(dataMatrix, design, outcomeName)
a
a <- computeAssocCountData(dataMatrixCount, design, outcomeName)
a

design$target <- Surv(time = runif(10), event = rbinom(n = 10, size = 1, prob = 0.65))
a <- computeAssocContinuousData(dataMatrix, design, outcomeName)
a
a <- computeAssocCountData(dataMatrixCount, design, outcomeName)
a
