
#set up
rm(list = ls())
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsPC_internal.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')
set.seed(12345)

#control panel
numMeasurements <- 5;
numDatasets <- 2;
combMethods <- c('Fisher', 'Liptak', 'Tippett');
numPerms <- 10000;
numCores <- 4;

#simulating data
numCombMethods <- length(combMethods)
measurements <- paste('m', 1:numMeasurements, sep = '');
pvaluesPerm <- array(runif(prod(c(numMeasurements, numDatasets, numPerms + 1))), 
                     dim = c(numMeasurements, numDatasets, numPerms + 1))
pvalues0 <- as.matrix(pvaluesPerm[ , , 1])
dataWeights <- rep(1, numDatasets)/numDatasets

#initializing
statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                  dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]))

#combining the pvalues
for(i in 1:numCombMethods){
  statsNPC[ , i, ] <-  apply(pvaluesPerm, 3, combiningPvalues, 
                             method = combMethods[i], dataWeights = dataWeights)
}

#ensuring no statistics are < 0
statsNPC <- statsNPC - min(statsNPC) + 1 

#computing the NPC p-values
pvaluesNPC <- statisticsToPvalues(statsNPC);
pvaluesNPC <- aperm(pvaluesNPC, c(2, 3, 1)) # magic numbers!
pvaluesNPC <- as.matrix(pvaluesNPC[, , 1])

#PC pvalues
pvaluesPC <- matrix(NA, numMeasurements, numCombMethods)
for(i in 1:numCombMethods){
  pvaluesPC[, i] <- combiningPvaluesParametric(pvalues0, combMethods[i])
}

x2 <- -2 * (log(pvalues0[, 1]) + log(pvalues0[, 2]))
fisherPvalues <- pchisq(x2, 4, lower.tail = FALSE)

pvalues0
pvaluesNPC
pvaluesPC
fisherPvalues
