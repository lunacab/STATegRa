# omicsNPC_internal function definition
# This function is not meant to be invoked directly by the user.
# See the omicsNPC function instead.

#' @importFrom utils combn
omicsNPC_internal <- function(dataMatrices, 
                              designs, 
                              dataMapping,
                              combMethods, 
                              functionsAnalyzingData, 
                              functionGeneratingIndex, 
                              outcomeName = NULL,
                              numPerms = 1000, numCores = 1,
                              allCombinations = FALSE,
                              dataWeights = rep(1, length(dataMatrices))/length(dataMatrices),
                              verbose = FALSE,
                              returnPermPvalues = FALSE,
                              ...){
  
  #### STEP 1 Compute initial statistics ####
  
  #message
  if(verbose){
      message("Compute initial statistics on data")
  }
  
  #information
  datasetsNames <- names(dataMatrices);
  numDataMatrices <- length(dataMatrices)
  numCombMethods <- length(combMethods);
  
  #can we run in parallel?
  parallelAvailable <- FALSE;
  
  # computing the statistics on the original data
  stats0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, dataMapping = dataMapping, functionsAnalyzingData = functionsAnalyzingData, outcomeName = outcomeName, ...)
  #system.time(stats0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, dataMapping = dataMapping, functionsAnalyzingData = functionsAnalyzingData, outcomeName = outcomeName))
  
  #information
  measurements <- dimnames(stats0)[[1]];
  numMeasurements <- length(measurements);
  
  #### STEP 2 Building NULL distributions by permuting data ####

  # message
  if(verbose){
      message("Building NULL distributions by permuting data")
  }

  # set function for permutation
  executePermutation <- function(...){
      
    #save(designs, outcomeName,  file = 'E:\\Dropbox (Personal)\\inprogress\\TwoPagesBionformatics\\3_newAnalysis_M\\safe1.RData')
    
    # permute designs
    permDesigns <- permutingData(designs = designs, functionGeneratingIndex = functionGeneratingIndex, outcomeName = outcomeName, ...)

    #save(designs, outcomeName, permDesigns, dataMapping, dataMatrices,  file = 'E:\\Dropbox (Personal)\\inprogress\\TwoPagesBionformatics\\3_newAnalysis_M\\safe2.RData')

    # recompute pvalues
    tempRes <- computeAssociation(dataMatrices, designs = permDesigns, dataMapping = dataMapping, functionsAnalyzingData = functionsAnalyzingData, outcomeName = outcomeName, ...)
    
    #return
    return(tempRes)
    
  }

  # try to run in parallel
  if (numCores != 1){
   
    # check if required packages are available
    if(requireNamespace("doSNOW", quietly = TRUE)){
    
        #creating the cluster
        cl <- parallel::makeCluster(numCores)
    
        #registering the cluster
        doSNOW::registerDoSNOW(cl)
    
        # compute DE expression
        statsPerm <- foreach(i = 1:numPerms, .packages = c('limma', 'survival', 'mdgsa')) %dopar% {
          source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_internal.R');
          source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_ancillaryFunctions.R');
          tryCatch({executePermutation(...)}, 
                   error = function(e){
                     tempRes <- matrix(NA, dim(dataMapping)[1], dim(dataMapping)[2] - 1);
                     rownames(tempRes) <- rownames(dataMapping);
                     colnames(tempRes) <- colnames(dataMapping)[1:(dim(dataMapping)[2] - 1)];
                     return(tempRes)
                     }
                   )
         
          }
        #system.time(statsPerm <- foreach(i = 1:numPerms, .packages = c('foreach', 'limma', 'survival', 'mdgsa')) %dopar% {executePermutation()})
        
        #stopping the cluster
        parallel::stopCluster(cl)
        
        #parallel is available
        parallelAvailable <- TRUE;
        
    }else{
      
      #not possible to run in parallel
      message("Install doSNOW to run omicsNPC in parallel...")
      message("omicsNPC will run using ONE core")
      statsPerm <- foreach(i = 1:numPerms) %do% {
      tryCatch({executePermutation(...)}, 
                 error = function(e){
                   tempRes <- matrix(NA, dim(dataMapping)[1], dim(dataMapping)[2] - 1);
                   rownames(tempRes) <- rownames(dataMapping);
                   colnames(tempRes) <- colnames(dataMapping)[1:(dim(dataMapping)[2] - 1)];
                   return(tempRes)
                 }
        )
      }

    }
    
  }else{
    
    #running serially
    statsPerm <- foreach(i = 1:numPerms) %do% {
      tryCatch({executePermutation(...)}, 
               error = function(e){
                 tempRes <- matrix(NA, dim(dataMapping)[1], dim(dataMapping)[2] - 1);
                 rownames(tempRes) <- rownames(dataMapping);
                 colnames(tempRes) <- colnames(dataMapping)[1:(dim(dataMapping)[2] - 1)];
                 return(tempRes)
               }
      )
    }

  }
   
  # combining original and permuted statistics in a single list
  statsPerm <- c(list(stats0), statsPerm)

  # transforming the list in a multi-dimensional array
  statsPerm <- array(unlist(statsPerm), 
                     dim = c(nrow(statsPerm[[1]]), ncol(statsPerm[[1]]), length(statsPerm)),
                     dimnames = list(rownames(statsPerm[[1]]),
                                     colnames(statsPerm[[1]]),
                                     1:length(statsPerm)))
  
  #### STEP 3 Compute p-values based on NULL distributions ####
  
  #message
  if(verbose){
      message("Compute pseudo p-values based on NULL distributions...")
  }
  
  # compute permutation pvalues based on stats
  pvaluesPerm <- statisticsToPvalues(statsPerm)
  
  #rearraging so that the order of dimension is measurements, datasets, permutations
  pvaluesPerm <- aperm(pvaluesPerm, c(2, 3, 1)) # magic numbers!
  
  #taking the partial p-values
  pvalues0 <- matrix(pvaluesPerm[ , , 1],
                       nrow = numMeasurements, 
                       ncol = numDataMatrices, 
                       dimnames = dimnames(pvaluesPerm)[1:2]);
  
  #### STEP 4 Compute NPC p-values ####
  
  #message
  if(verbose){
      message("NPC p-values calculation...")
  }
  
  #all combination of datasets?
  if(allCombinations){
    
    #which are the combinations of all data matrices?
    combinations <- list();
    for(i in 2:numDataMatrices){
      combinations <- c(combinations, combn(1:numDataMatrices, i, simplify = FALSE));
    }
    numCombinations <- length(combinations);
    
    #initializing 
    dataMappings <- vector('list', numCombinations);
    pvaluesNPC <- vector('list', numCombinations);
    
    #for each combination, let's compute the NPC stats and p-values
    for(j in 1:numCombinations){
      
      #initializing
      statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                        dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]));
      
      #is parallel available?
      if(parallelAvailable){
        
        #creating the cluster
        cl <- parallel::makeCluster(numCores)
        
        #registering the cluster
        doSNOW::registerDoSNOW(cl)   
        
        #combining the pvalues. Here we choose which data matrice to combine
        for(i in 1:numCombMethods){
          statsNPC[ , i, ] <- foreach(k = 1:dim(pvaluesPerm)[3], .combine = 'cbind') %dopar% {
            source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_internal.R');
            source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_ancillaryFunctions.R');
            combiningPvalues(as.matrix(pvaluesPerm[ , combinations[[j]], k]), method = combMethods[i], dataWeights = dataWeights);
          }
        }
        
        #closing the cluster
        stopCluster(cl)
        
      }else{
        
        #combining the pvalues. Here we choose which data matrice to combine
        for(i in 1:numCombMethods){
          statsNPC[ , i, ] <-  apply(pvaluesPerm[ , combinations[[j]], ], 3, combiningPvalues, method = combMethods[i], dataWeights = dataWeights);
        }
        
      }
      
      #let's ensure that all statistics are larger than 0
      #this does not affect the computation of the pvalues
      finiteStats <- !is.infinite(statsNPC);
      minStats <- min(statsNPC[finiteStats], na.rm = TRUE);
      maxStats <- max(statsNPC[finiteStats], na.rm = TRUE);
      statsNPC[!finiteStats & sign(statsNPC) < 0] <- minStats - 1;
      statsNPC[!finiteStats & sign(statsNPC) > 0] <- maxStats + 1;
      minStats <- min(statsNPC, na.rm = TRUE);
      statsNPC <- statsNPC - min(statsNPC) + 1 
      
      #computing the NPC p-values
      pvaluesNPCTemp <- statisticsToPvalues(statsNPC);
      pvaluesNPCTemp <- aperm(pvaluesNPCTemp, c(2, 3, 1)) # magic numbers!
      
      #creating the data mapping for the combination at hand
      dataMappingTmp <- dataMapping[, combinations[[j]]];
      toKeep <- !duplicated(dataMappingTmp) & !apply(dataMappingTmp, 1, function(x){all(is.na(x))})
      dataMappingTmp <- cbind(dataMappingTmp[toKeep, ], dataMapping[toKeep, ncol(dataMapping)]) #adding the reference
      rownames(dataMappingTmp) <- apply(dataMappingTmp, 1, function(x){paste(x, collapse = ':')});
      
      #keeping only the relevant pvalues
      if(!returnPermPvalues){
        pvaluesNPC[[j]] <- matrix(pvaluesNPCTemp[toKeep, , 1], 
                             nrow = sum(toKeep), 
                             ncol = numCombMethods, 
                             dimnames = list(dimnames(pvaluesNPCTemp)[[1]][toKeep], dimnames(pvaluesNPCTemp)[[2]]))
        rownames(pvaluesNPC[[j]]) <- rownames(dataMappingTmp);
      }else{
        pvaluesNPC[[j]] <- pvaluesNPCTemp[toKeep, , ];
        dimnames(pvaluesNPC[[j]])[[1]] <- rownames(dataMappingTmp);
      }
      
      #storing and naming the results
      dataMappings[[j]] <- dataMappingTmp;
      names(pvaluesNPC)[j] <- paste(datasetsNames[combinations[[j]]], collapse = '_');
      names(dataMappings)[j] <- names(pvaluesNPC)[j];
      
    }
      #not all combinations
  }else{
    
    #initializing
    statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                      dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]));
    
    #is parallel available ?
    if(parallelAvailable){
      
      #creating the cluster
      cl <- parallel::makeCluster(numCores)
      
      #registering the cluster
      doSNOW::registerDoSNOW(cl)   
      
      for(i in 1:numCombMethods){
        statsNPC[ , i, ] <- foreach(k = 1:dim(pvaluesPerm)[3], .combine = 'cbind') %dopar% {
          source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_internal.R');
          source('../../omicsNPC/code/stategra/R/STATegRa_omicsNPC_ancillaryFunctions.R');
          combiningPvalues(as.matrix(pvaluesPerm[ , , k]), method = combMethods[i], dataWeights = dataWeights);
        }
      }
      
      #closing the cluster
      stopCluster(cl)
      
    }else{
      #combining the pvalues
      for(i in 1:numCombMethods){
        statsNPC[ , i, ] <-  apply(pvaluesPerm, 3, combiningPvalues, 
                                   method = combMethods[i], dataWeights = dataWeights)
      }
      
    }
    
    #let's ensure that all statistics are larger than 0
    #this does not affect the computation of the pvalues
    finiteStats <- !is.infinite(statsNPC);
    minStats <- min(statsNPC[finiteStats], na.rm = TRUE);
    maxStats <- max(statsNPC[finiteStats], na.rm = TRUE);
    statsNPC[!finiteStats & sign(statsNPC) < 0] <- minStats - 1;
    statsNPC[!finiteStats & sign(statsNPC) > 0] <- maxStats + 1;
    minStats <- min(statsNPC, na.rm = TRUE);
    statsNPC <- statsNPC - min(statsNPC) + 1 
    
    #computing the NPC p-values
    pvaluesNPC <- statisticsToPvalues(statsNPC);
    pvaluesNPC <- aperm(pvaluesNPC, c(2, 3, 1)) # magic numbers!
    
    #keeping only the relevant pvalues, if not specified otherwise
    if(!returnPermPvalues){
      pvaluesNPC <- matrix(pvaluesNPC[ , , 1], 
                            nrow = numMeasurements, 
                            ncol = numCombMethods, 
                            dimnames = dimnames(pvaluesNPC)[1:2])
    }

  }
  
  #creating the object to return
  #Note: for the moment is only a list. It will be a proper object in a next release
  #pvaluesNPC: either a matrix (allCombinations = FALSE) or a list (allCombinations = TRUE)
  toReturn <- list(stats0 = stats0, pvalues0 = pvalues0, pvaluesNPC = pvaluesNPC);
  if(returnPermPvalues){
    toReturn$pvaluesPerm = pvaluesPerm;
  }
  if(allCombinations){
    toReturn$dataMappings <- dataMappings;
  }
  
  return(toReturn)

}

#### Functions to run omicsNPC. These should not be changed. ####

# This method applies the functions for computing the association between measurements and outcome
computeAssociation <- function(dataMatrices, designs, dataMapping, functionsAnalyzingData, outcomeName, ...){
  
  #applying to each data matrix the corresponding function
  results <- vector('list', length(dataMatrices));
  for(j in 1:length(dataMatrices)){
    results[[j]] <- functionsAnalyzingData[[j]](dataMatrices[[j]], designs[[j]], outcomeName, ...);
  }
  
  #Return the results as a matrix: we do now assume the same number of results for each dataset; only a coherent naming
  
  #matrix to return
  toReturn <- matrix(NA, nrow = dim(dataMapping)[1], ncol = length(dataMatrices));
  colnames(toReturn) <- names(dataMatrices);
  
  #filling the matrix
  for(j in 1:length(results)){
    toReturn[ , j] <- results[[j]][ as.character(dataMapping[[j]]) ];
  }
  rownames(toReturn) <- rownames(dataMapping);
  
  #return
  return(toReturn)
  
}

# This method permutes the samples likewise across datasets
permutingData <- function(designs, functionGeneratingIndex, outcomeName, ...){
  
  #if the outcomeName is null, we assume it to be the last column of the design matrices
  if(is.null(outcomeName)){
    outcomeName <- names(designs[[1]])[ncol(designs[[1]])]; #this is a bit of a hack...
  }
  
  #Which columns are common in all designs?
  commonElements <- colnames(designs[[1]]);
  for(j in 2:length(designs)){
    commonElements <- intersect(commonElements, colnames(designs[[j]]));
  }
  
  #Extracting the common elements. 
  commonDesign <- designs[[1]][commonElements];
  rowNamesCommonDesign <- rownames(designs[[1]])
  for(j in 2:length(designs)){
    commonDesign <- rbind(commonDesign, designs[[j]][commonElements], make.row.names = FALSE);
    rowNamesCommonDesign <- c(rowNamesCommonDesign, rownames(designs[[j]]));
  }
  toEliminate <- which(duplicated(rowNamesCommonDesign));
  if(length(toEliminate) > 0){
    commonDesign <- as.data.frame(commonDesign[-toEliminate, ]); 
    colnames(commonDesign) <- commonElements;
    rowNamesCommonDesign <- rowNamesCommonDesign[-toEliminate];
  }
  rownames(commonDesign) <- rowNamesCommonDesign;
  
  # computing the permutation index
  index <- functionGeneratingIndex(commonDesign)
  permSampleNames <- rownames(commonDesign)[index];
  
  #getting the indexes for each matrix. Those indexes are actually sample names
  indexes <- vector('list', length(designs));
  for(j in 1:length(designs)){
    currentSamples <- rownames(designs[[j]]);
    toEliminate <- which(! permSampleNames %in% currentSamples);
    if(length(toEliminate)>0){
      tempIndex <- permSampleNames[-toEliminate]; 
    }else{
      tempIndex <- permSampleNames
    }
    indexes[[j]] <- tempIndex;
  }
  
  #permuting the outcome / experimental factor of each design matrix
  permDesigns <- designs;
  for(j in 1:length(designs)){
    tempDesign <- designs[[j]];
    outcomeId <- which(colnames(tempDesign) == outcomeName);
    tempDesign[ , outcomeId] <- tempDesign[indexes[[j]], outcomeId];
    permDesigns[[j]] <- tempDesign;
  }
  
  #returning the permuted designs
  return(permDesigns)
  
}

#function for computing the p-value for each element of a vector of statistics
computePvaluesVect <- function(statsVect){        
  
  #taking the absolute value of the statistics
  statsVect <- abs(statsVect);
  
  #how many values are present in the vector of statistics?
  numStatistics <- sum(!is.na(statsVect))
  
  #the rank function with ties.method = 'min' indicates how many statistics are 
  #larger or equal than each value in the vector
  #Note: we assume that the greater the statistics, the more significant the finding
  #This is the case with, for example, chi-square statistics
  #This is not the case for on-tailed t-test where more negative statistics are actually more extreme
  #The functions computing the statistics should be programmed accordingly; 
  #for example, by changing the sign of negative statistics and putting to zero the positive ones, 
  #or by returning the -log(desiredOneTailedPvalue)
  #THIS IS NOT CHANGED: we take the absolute value of the statistics! See above
  numLargerValues <- numStatistics - rank(statsVect, na.last = "keep", ties.method = "min");
  
  #correction for avoiding zero p-values
  pvaluesVect <- (numLargerValues + 1)/(numStatistics);
  
  #returning the p-values
  return(pvaluesVect)
  
}  

# function for transforming an multi-dimensional array of statistics to pvalues
statisticsToPvalues <- function(stats){
  
  #stats is supposed to be an array where the permutations are stored in the last dimension
  
  #in case no dimension are defined
  if (is.null(dim(stats))) {
    stats <- array(stats, dim = c(1, length(stats)))
  }
  
  #how many dimensions?
  numDims <- length(dim(stats));
  
  #computing the p-values
  pvalues <- apply(stats, 1:(numDims-1), computePvaluesVect)
  
  #returning the values
  return(pvalues)
}

# function for combining p-values
combiningPvalues <- function(pvalues, method, dataWeights){
  
  #if no valid p-values, then return NA
  
  #depending by the method, a different combination function is applied
  if(method == "Fisher"){
    stats <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                -2*sum( dataWeights[idx] * log(x[idx]) ) })
  }else if(method == "Liptak"){
    stats <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));                          
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                sum( dataWeights[idx] * qnorm(1 - x[idx]) ) })
  }else if(method == "Tippett"){
    stats <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                max( dataWeights[idx] * (1 - x[idx]) ) })
  }
  
  #returning the values
  return(stats)
  
}

#function for merging p-values from different permutation runs
combinePermPvalues <- function(p1, p2, numPerms1, numPerms2){
  
  #this function assumes that both p1 and p2 were produced as p = (numLargerStat + 1) / (numPerms + 1)
  
  #num of original stats larger than stat0
  nS1 <- p1 * (numPerms1 + 1) - 1;
  nS2 <- p2 * (numPerms2 + 1) - 1;
  nS <- nS1 + nS2;
  
  #merged p-values
  p <- (nS + 1) / (numPerms1 + numPerms2 + 1)
  
  #return output
  return(p)
  
}