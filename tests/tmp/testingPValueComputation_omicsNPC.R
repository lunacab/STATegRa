
#set up
rm(list = ls())
source('../../R/STATegRa_omicsNPC_internal.R')
source('../../R/STATegRa_omicsPC_internal.R')
source('../../R/STATegRa_omicsNPC_ancillaryFunctions.R')
library(doSNOW)
library(foreach)

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
dataWeights <- rep(1, numDatasets)/numDatasets

#initializing
statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                      dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]))
statsNPC_par <- statsNPC

#combining the pvalues
etime <- system.time(
  for(i in 1:numCombMethods){
    statsNPC[ , i, ] <-  apply(pvaluesPerm, 3, combiningPvalues, 
                               method = combMethods[i], dataWeights = dataWeights)
  }
)
print(etime)

#combining the pvalues parallel
etime_par <- system.time({
  
  #creating the cluster
  cl <- parallel::makeCluster(numCores)
  
  #registering the cluster
  doSNOW::registerDoSNOW(cl)   
  
  for(i in 1:numCombMethods){
    
    statsNPC_par[ , i, ] <- foreach(j = 1:dim(pvaluesPerm)[3], .combine = 'cbind') %dopar% {
      combiningPvalues(as.matrix(pvaluesPerm[ , , j]), method = combMethods[i], dataWeights = dataWeights);
    }
    
  }
  
  #closing the cluster
  stopCluster(cl)
  
})
print(etime_par)

#ensuring all stats are positive
statsNPC <- statsNPC - min(statsNPC) + 1
statsNPC_par <- statsNPC_par - min(statsNPC_par) + 1

#computing the NPC p-values
etime2 <- system.time(pvaluesNPC <- statisticsToPvalues(statsNPC));
print(etime2)
etime2 <- system.time(pvaluesNPC_par <- statisticsToPvalues(statsNPC_par));
print(etime2)

