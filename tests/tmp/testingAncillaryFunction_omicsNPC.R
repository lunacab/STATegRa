
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

#### projecting ####

set.seed(1234)

dataMatrix <- matrix(rnorm(30), 3, 10);
design <- data.frame(matrix(rnorm(20), 10, 2))
rownames(dataMatrix) <- c('a','b', 'd');
outcomeName <- 'X1'
mapping <- list(a = paste('V', 1:33, sep = ''), 
                b = paste('V', 34:70, sep = ''), 
                d = paste('V', 71:100, sep = ''))

a <- computeProjectedAssocContinuousData(dataMatrix, design, outcomeName, mapping = mapping)
class(a)
a

dataMatrixCount <- matrix(rnbinom(n = 30, size = 1000, prob = 0.2), 3, 10);
rownames(dataMatrixCount) <- c('a','b', 'd');
a <- computeProjectedAssocCountData(dataMatrixCount, design, outcomeName, mapping = mapping)
a

design <- cbind(design, target = factor(rep(c(1,2), 5)))
outcomeName <- 'target'
a <- computeProjectedAssocContinuousData(dataMatrix, design, outcomeName, mapping = mapping)
a
a <- computeProjectedAssocCountData(dataMatrixCount, design, outcomeName, mapping = mapping)
a

levels(design$target) <- 1:3;
design$target[7:10] <- 3
a <-computeProjectedAssocContinuousData(dataMatrix, design, outcomeName, mapping = mapping)
a
a <- computeProjectedAssocCountData(dataMatrixCount, design, outcomeName, mapping = mapping)
a

design$target <- Surv(time = runif(10), event = rbinom(n = 10, size = 1, prob = 0.65))
a <- computeProjectedAssocContinuousData(dataMatrix, design, outcomeName, mapping = mapping)
a
a <- computeProjectedAssocCountData(dataMatrixCount, design, outcomeName, mapping = mapping)
a
