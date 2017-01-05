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
  
  # computing the statistics on the original data
  stats0 <- computeAssociation(dataMatrices = dataMatrices, designs = designs, dataMapping = dataMapping, functionsAnalyzingData = functionsAnalyzingData, outcomeName = outcomeName, ...)
  
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
      
    # permute designs
    permDesigns <- permutingData(designs = designs, functionGeneratingIndex = functionGeneratingIndex, outcomeName = outcomeName, ...)


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
        statsPerm <- foreach(i = 1:numPerms, .packages = 'limma') %dopar% {executePermutation(...)}

        #stopping the cluster
        parallel::stopCluster(cl)
        
    }else{
      
      #not possible to run in parallel
      message("Install doSNOW to run omicsNPC in parallel...")
      message("omicsNPC will run using ONE core")
      statsPerm <- foreach(i = 1:numPerms) %do% executePermutation(...)

    }
    
  }else{
    
    #running serially
    statsPerm <- foreach(i = 1:numPerms) %do% executePermutation(...)

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
  #keeping only the relevant pvalues
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
    pvaluesNPC <- vector('list', numCombinations);
    
    #for each combination, let's compute the NPC stats and p-values
    for(j in 1:numCombinations){
      
      #initializing
      statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                        dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]));
      
      #combining the pvalues. Here we choose which data matrice to combine
      for(i in 1:numCombMethods){
        statsNPC[ , i, ] <-  apply(pvaluesPerm[ , combinations[[j]], ], 3, combiningPvalues, method = combMethods[i], dataWeights = dataWeights);
      }
      
      #computing the NPC p-values
      pvaluesNPCTemp <- statisticsToPvalues(statsNPC);
      pvaluesNPCTemp <- aperm(pvaluesNPCTemp, c(2, 3, 1)) # magic numbers!
      
      #keeping only the relevant pvalues, if not specified otherwise
      if(!returnPermPvalues){
        pvaluesNPC[[j]] <- matrix(pvaluesNPCTemp[ , , 1], 
                             nrow = numMeasurements, 
                             ncol = numCombMethods, 
                             dimnames = dimnames(pvaluesNPCTemp)[1:2])
      }
      
      #naming the list component
      names(pvaluesNPC)[j] <- paste(datasetsNames[combinations[[j]]], collapse = '_');
      
    }
      
  }else{
    
    #initializing
    statsNPC <- array(dim = c(numMeasurements, numCombMethods, numPerms + 1), 
                      dimnames = list(measurements, combMethods, dimnames(pvaluesPerm)[[3]]));
    
    #combining the pvalues
    for(i in 1:numCombMethods){
      statsNPC[ , i, ] <-  apply(pvaluesPerm, 3, combiningPvalues, 
                                 method = combMethods[i], dataWeights = dataWeights)
    }
    
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
  if(returnPermPvalues){
    toReturn <- list(stats0 = stats0, pvalues0 = pvalues0, pvaluesNPC = pvaluesNPC, pvaluesPerm = pvaluesPerm);
  }else{
    toReturn <- list(stats0 = stats0, pvalues0 = pvalues0, pvaluesNPC = pvaluesNPC);
  }
  
  return(toReturn)

}

#### Ancillary functions. The user may provide different ones ####

#this fucntion computes association for continuous data
computeAssocContinuousData <- function(dataMatrix, design, outcomeName, ...){

  #ensuring the outcome is the last column of the design
  if(!is.null(outcomeName)){
    
    #check
    if(!any(colnames(design) == outcomeName)){
      stop(paste('Ensure', outcomeName, 'is present in all design matrices'));
    }
    
    #moving the outcome to last position
    outcomeTemp <- design[[outcomeName]];
    design[[outcomeName]] <- NULL;
    design <- cbind(design, outcomeTemp);
    names(design)[ncol(design)] <- outcomeName;
    
  }
  
  #we assume that the outcome is the last column
  outcomeId <- ncol(design);

  #what type of outcome?
  multiGroups <- FALSE;
  if((is.factor(design[[outcomeId]]) & length(levels(design[[outcomeId]])) >= 3) ||
     (is.character(design[[outcomeId]]) & length(unique(design[[outcomeId]])) >= 3)){
    multiGroups <- TRUE;
  }
  
  #creating the model matrix 
  numGroups <- length(unique(design[[outcomeId]]));
  modelFormula <- as.formula(paste('~', paste(colnames(design), collapse = '+')))
  design <- model.matrix(modelFormula, design)
  numCols <- ncol(design)
  
  #if there are n groups, the last n-1 columns represent them
  numColGroups <- (ncol(design) - (numGroups - 1) + 1):ncol(design)
  
  #fitting the models
  fit <- lmFit(dataMatrix, design = design)
  fit2 <- eBayes(fit)
  
  #if not multiple groups
  if(!multiGroups){
    results <- topTable(fit2, numCols, number = Inf, sort.by = 'none')
  }else{
    results <- topTable(fit2, numColGroups, number = Inf, sort.by = 'none')
  }
  
  #returning statistics (log transformed p-values)
  logPvalues <- (-1) * log10(results$P.Value);
  names(logPvalues) <- rownames(results);
  return(logPvalues)

}

#this fucntion computes association for count data
computeAssocCountData <- function(dataMatrix, design, outcomeName, ...){

  
  #ensuring the outcome is the last column of the design
  if(!is.null(outcomeName)){
    
    #check
    if(!any(colnames(design) == outcomeName)){
      stop(paste('Ensure', outcomeName, 'is present in all design matrices'));
    }
    
    #moving the outcome to last position
    outcomeTemp <- design[[outcomeName]];
    design[[outcomeName]] <- NULL;
    design <- cbind(design, outcomeTemp);
    names(design)[ncol(design)] <- outcomeName;
    
  }
  
  #we assume that the outcome is the last column
  outcomeId <- ncol(design);
  
  #what type of outcome?
  multiGroups <- FALSE;
  if((is.factor(design[[outcomeId]]) & length(levels(design[[outcomeId]])) >= 3) ||
     (is.character(design[[outcomeId]]) & length(unique(design[[outcomeId]])) >= 3)){
    multiGroups <- TRUE;
  }
  
  #creating the model matrix 
  numGroups <- length(unique(design[[outcomeId]]));
  modelFormula <- as.formula(paste('~', paste(colnames(design), collapse = '+')))
  design <- model.matrix(modelFormula, design)
  numCols <- ncol(design)
  
  #if there are n groups, the last n-1 columns represent them
  numColGroups <- (ncol(design) - (numGroups - 1) + 1):ncol(design)
  
  #applying voom
  nf = calcNormFactors(dataMatrix, method = "TMM")
  voom.data = voom(dataMatrix, design = design,
                   lib.size = colSums(dataMatrix) * nf)
  voom.data$genes = rownames(dataMatrix)  
  
  #fitting the models
  fit <- lmFit(voom.data, design = design)
  fit2 <- eBayes(fit)
  
  #if not multiple groups
  if(!multiGroups){
    results <- topTable(fit2, numCols, number = Inf, sort.by = 'none')
  }else{
    results <- topTable(fit2, numColGroups, number = Inf, sort.by = 'none')
  }
  
  #returning statistics (log transformed p-values)
  logPvalues <- (-1) * log10(results$P.Value);
  names(logPvalues) <- rownames(results);
  return(logPvalues)
  
}

#function for generating random, iid permutation indices
generate_iid_data_index <- function(design){
  numSamples <- dim(design)[1];
  index <- sample(x = 1:numSamples, size = numSamples, replace = FALSE)
  return(index)
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
  
  #how many values are present in the vector of statistics?
  numStatistics <- sum(!is.na(statsVect))
  
  #the rank function with ties.method = 'min' indicates how many statistics are 
  #larger or equal than each value in the vector
  #Note: we assume that the greater the statistics, the more significant the finding
  #This is the case with, for example, chi-square statistics
  #This is not the case for on-tailed t-test where more negative statistics are actually more extreme
  #The functions computing the statistics should be programmed accordingly; 
  #for example, by changing the sign of negative statistics and putting to zero positive once, 
  #or by returning the -log(desiredOneTailedPvalue)
  numLargerValues <- numStatistics - rank(statsVect, na.last = "keep", ties.method = "min");
  
  #correction for avoiding zero p-values
  pvaluesVect <- (numLargerValues + 1)/(numStatistics + 1);
  
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