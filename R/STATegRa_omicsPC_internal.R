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
  
  # computing the p-values on the original data
  tmp <- computeAssociation(dataMatrices = dataMatrices, designs = designs, 
                               dataMapping = dataMapping, 
                               functionsAnalyzingData = functionsAnalyzingData, 
                               outcomeName = outcomeName,
                               returnPValues = TRUE, ...)
  # tmp <- computeAssociation(dataMatrices = dataMatrices, designs = designs,
  #                           dataMapping = dataMapping,
  #                           functionsAnalyzingData = functionsAnalyzingData,
  #                           outcomeName = outcomeName,
  #                           returnPValues = TRUE)
  pvalues0 <- tmp$stats
  results0 <- tmp$results
  
  #information
  measurements <- dimnames(pvalues0)[[1]];
  numMeasurements <- length(measurements);
  
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
    adjPvaluesPC <- vector('list', numCombinations);
    
    #for each combination, let's compute the PC stats and p-values
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
      dataMappingTmp <- dataMappingTmp[toKeep, ]
      pvaluesPCTemp <- pvaluesPCTemp[toKeep, ]
      rownames(dataMappingTmp) <- NULL
      rownames(pvaluesPCTemp) <- NULL
      
      #storing and naming the results
      pvaluesPC[[j]] <- pvaluesPCTemp;
      dataMappings[[j]] <- dataMappingTmp;
      names(pvaluesPC)[j] <- paste(datasetsNames[combinations[[j]]], collapse = '_');
      names(dataMappings)[j] <- names(pvaluesPC)[j];
      
      #adjusted pvalues
      adjPvaluesPC[[j]] <- apply(pvaluesPC[[j]], 2, function(x){p.adjust(x, method = 'fdr')})
      names(adjPvaluesPC)[j] <- paste(datasetsNames[combinations[[j]]], collapse = '_');
      
      #adding data mapping to pvaluesPC
      pvaluesPC[[j]] <- cbind(dataMappings[[j]], pvaluesPC[[j]])
      adjPvaluesPC[[j]] <- cbind(dataMappings[[j]], adjPvaluesPC[[j]])
      rownames(pvaluesPC[[j]]) <- NULL;
      rownames(adjPvaluesPC[[j]]) <- NULL;
      
    }
      
  }else{
    
    #combining the pvalues
    pvaluesPC <- matrix(NA, nrow = dim(dataMapping)[1], ncol = numCombMethods)
    rownames(pvaluesPC) <- rownames(dataMapping);
    colnames(pvaluesPC) <- combMethods;
    for(i in 1:numCombMethods){
      pvaluesPC[ , i] <-  combiningPvaluesParametric(pvalues = pvalues0, method = combMethods[i])
    }
    
    #adjusting pvalues
    adjPvaluesPC <- apply(pvaluesPC, 2, function(x){p.adjust(x, method = 'fdr')})
    
    #adding data mapping to pvaluesPC
    pvaluesPC <- cbind(dataMapping[, 1:length(dataMatrices)], pvaluesPC)
    adjPvaluesPC <- cbind(dataMapping[, 1:length(dataMatrices)], adjPvaluesPC)

  }
  
  #creating the object to return
  toReturn <- results0;
  toReturn$pvaluesPC <- pvaluesPC;
  toReturn$qvaluesPC <- adjPvaluesPC;
  
  #returning
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