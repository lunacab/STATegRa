# omicsPC_internal function definition
# This function is not meant to be invoked directly by the user.
# See the omicsPC function instead.

#' @importFrom utils combn
omicsPC_internal <- function(dataMatrices, 
                              designs, 
                              dataMapping,
                              combMethods, 
                              functionsAnalyzingData, 
                              outcomeName = NULL,
                              allCombinations = FALSE,
                              verbose = FALSE,
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
  
  # computing the statistics on the original data
  stats0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, 
                               dataMapping = dataMapping, 
                               functionsAnalyzingData = functionsAnalyzingData, 
                               outcomeName = outcomeName,
                               returnPValues = FALSE, ...)
  # system.time(stats0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, 
  #                                          dataMapping = dataMapping, 
  #                                          functionsAnalyzingData = functionsAnalyzingData, 
  #                                          outcomeName = outcomeName,
  #                                          returnPValues = FALSE))
  
  #information
  measurements <- dimnames(stats0)[[1]];
  numMeasurements <- length(measurements);
  
  # computing the p-values on the original data
  pvalues0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, 
                                 dataMapping = dataMapping, 
                                 functionsAnalyzingData = functionsAnalyzingData, 
                                 outcomeName = outcomeName,
                                 returnPValues = TRUE, ...)
  
  #### STEP 2 Compute PC p-values ####
  
  #message
  if(verbose){
      message("PC p-values calculation...")
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
    pvaluesPC <- vector('list', numCombinations);
    
    #for each combination, let's compute the NPC stats and p-values
    for(j in 1:numCombinations){
      
      #combining the pvalues. Here we choose which data matrices to combine
      pvaluesPCTemp <- matrix(NA, nrow = dim(dataMapping)[1], ncol = numCombMethods)
      colnames(pvaluesPCTemp) <- combMethods;
      for(i in 1:numCombMethods){
        pvaluesPCTemp[ , i] <-  combiningPvaluesParametric(pvalues = pvalues0[ , combinations[[j]]], method = combMethods[i]);
      }
      
      #creating the data mapping for the combination at hand
      dataMappingTmp <- dataMapping[, combinations[[j]]];
      toKeep <- which(!duplicated(dataMappingTmp) & !apply(dataMappingTmp, 1, function(x){all(is.na(x))}))
      dataMappingTmp <- cbind(dataMappingTmp[toKeep, ], dataMapping[toKeep, ncol(dataMapping)]) #adding the reference
      colnames(dataMappingTmp)[ncol(dataMappingTmp)] <- colnames(dataMapping)[ncol(dataMapping)]
      pvaluesPCTemp <- pvaluesPCTemp[toKeep, ]
      rownames(dataMappingTmp) <- apply(dataMappingTmp, 1, function(x){paste(x, collapse = ':')});
      rownames(pvaluesPCTemp) <- rownames(dataMappingTmp);
      
      #storing and naming the results
      pvaluesPC[[j]] <- pvaluesPCTemp;
      dataMappings[[j]] <- dataMappingTmp;
      names(pvaluesPC)[j] <- paste(datasetsNames[combinations[[j]]], collapse = '_');
      names(dataMappings)[j] <- names(pvaluesPC)[j];
      
    }
      
  }else{
    
    #combining the pvalues
    pvaluesPC <- matrix(NA, nrow = dim(dataMapping)[1], ncol = numCombMethods)
    rownames(pvaluesPC) <- rownames(dataMapping);
    colnames(pvaluesPC) <- combMethods;
    for(i in 1:numCombMethods){
      pvaluesPC[ , i] <-  combiningPvaluesParametric(pvalues = pvalues0, method = combMethods[i])
    }

  }
  
  #creating the object to return
  #Note: for the moment is only a list. It will be a proper object in a next release
  #pvaluesPC: either a matrix (allCombinations = FALSE) or a list (allCombinations = TRUE)
  if(allCombinations){
    toReturn <- list(stats0 = stats0, pvalues0 = pvalues0, 
                     pvaluesPC = pvaluesPC, dataMappings = dataMappings);
  }else{
    toReturn <- list(stats0 = stats0, pvalues0 = pvalues0, pvaluesPC = pvaluesPC);
  }
  
  return(toReturn)

}


# function for combining p-values
combiningPvaluesParametric <- function(pvalues, method){
  
  #depending by the method, a different combination function is applied
  if(method == "Fisher"){
    combPvalues <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                stat <- -2*sum(log(x[idx]));
                                1 - pchisq(stat, 2 * length(idx)) })
  }else if(method == "Liptak"){
    combPvalues <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));                          
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                stat <- sum(qnorm(1 - x[idx]));
                                1 - pnorm(stat/sqrt(length(idx))) })
  }else if(method == "Tippett"){
    combPvalues <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                min(1, length(idx) * min(x[idx])) })
  }else if(method == "Simes" || method == "Benjamini"){
    combPvalues <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                x <- sort(x[idx], decreasing = FALSE);
                                min(length(x) * x / (1:length(x))) })
  }else if(method == "Sidak"){
    combPvalues <- apply(pvalues,1,function(x){
                                idx <- which(!is.na(x));
                                if(length(idx) == 0){ #in case all p-values are NA
                                  return(NA);
                                }
                                1 - (1 - min(x[idx])) ^ length(idx) })
  }
  
  #returning the values
  return(combPvalues)
  
}