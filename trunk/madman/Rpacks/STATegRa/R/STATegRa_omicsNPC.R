#' @export
#' @import affy
#' @import edgeR
#' @import limma
#' @import foreach
#' @title omicsNPC, applying the Non-Parametric Combination (NPC) on omics datasets
#' @aliases omicsNPC,list,data.frame-method omicsNPC,list,missing-method
#' @description
#' This function applies the NonParametric Combination methodology on the integrative analysis 
#' of different omics data modalities.
#' It retrieves genes associated to a given outcome, taking into account all omics data.
#' First, each datatype is analyzed independently using the appropriate method.
#' omicsNPC analyses continuous data (microarray) using limma, while count data (e.g., RNAseq) 
#' are first preprocessed with using the "voom" function. The user can also specify their own 
#' function for computing deregulation / association 
#' The p-values from the single dataset analysis are then combined employing Fisher, 
#' Liptak and Tippett combining functions.
#' The Tippett function returns findings which are supported by at least one omics modality.
#' The Liptak function returns findings which are supportd by most modalities.
#' The Fisher function has an intermediate behavior between those of Tippett and Liptak.
#' @usage omicsNPC(dataInput, dataMapping, dataTypes = rep('continuous', length(dataInput)), 
#'                combMethods = c("Fisher", "Liptak", "Tippett"), numPerms = 1000, 
#'                numCores = 1, verbose = FALSE, functionGeneratingIndex = NULL, 
#'                outcomeName = NULL, allCombinations = FALSE, 
#'                dataWeights = rep(1, length(dataInput))/length(dataInput), 
#'                returnPermPvalues = FALSE, ...)
#'
#' @param dataInput List of ExpressionSet objects, one for each data modality. 
#' @param dataMapping A data frame describing how to map measurements across datasets. See details for more information.
#' @param dataTypes Character vector with possible values: 'continuous', 'count'. Alternatively, a list of functions for assessing deregulation / association with an outcome
#' @param combMethods Character vector with possible values: 'Fisher', 'Liptak', 'Tippett'. If more than one is specified, all will be used.
#' @param numPerms Number of permutations
#' @param numCores Number of CPU cores to use
#' @param verbose Logical, if set to TRUE omicsNPC prints out the step that it performs
#' @param functionGeneratingIndex Function generating the indices for randomly permuting the samples
#' @param outcomeName Name of the outcome of interest / experimental factor, as reported in the design matrices. If NULL, the last column of the design matrices is assumed to be the outcome of interest.
#' @param allCombinations Logical, if TRUE all combinations of omics datasets are considered
#' @param dataWeights A vector specifying the weigth to give to each dataset. Note that sum(dataWeights) should be 1.
#' @param returnPermPvalues Logical, should the p-values computed at each permutation being returned?
#' @param ... Additional arguments to be passed to the user-defined functions
#'
#' @return A list containing:
#' stats0 Partial deregulation / association statistics 
#' pvalues0 The partial p-values computed on each dataset
#' pvaluesNPC The p-values computed through NPC.
#' permPvalues The p-values computed at each permutation
#'
#' @author Nestoras Karathanasis, Vincenzo Lagani
#'
#' @references
#' Pesarin, Fortunato, and Luigi Salmaso. Permutation tests for complex data: theory, applications and software.
#' John Wiley & Sons, 2010.
#' Nestoras Karathanasis, Ioannis Tsamardinos and Vincenzo Lagani. omicsNPC: applying the Non-Parametric Combination 
#' methodology to the integrative analysis of heterogeneous omics data. PlosONE 11(11): e0165545. doi:10.1371/journal.pone.0165545
#'
#' @examples
#' # Load the data
#' data("TCGA_BRCA_Batch_93")
#' # Setting dataTypes, the first two ExpressionSets include RNAseq data,
#' # the third ExpressionSet includes Microarray data.
#' dataTypes <- c("count", "count", "continuous")
#' # Setting methods to combine pvalues
#' combMethods = c("Fisher", "Liptak", "Tippett")
#' # Setting number of permutations
#' numPerms = 1000
#' # Setting number of cores
#' numCores = 1
#' # Setting omicsNPC to print out the steps that it performs.
#' verbose = TRUE
#' # Run omicsNPC analysis.
#' # The output contains a data.frame of p-values, where each row corresponds to a gene, 
#' # and each column corresponds to a method used in the analysis.
#' 
#' \dontrun{out <- omicsNPC(dataInput = TCGA_BRCA_Data, dataTypes = dataTypes,
#'                             combMethods = combMethods, numPerms = numPerms,
#'                             numCores = numCores, verbose = verbose)}

setGeneric(name="omicsNPC",
           def=function(dataInput,
                        dataMapping,
                        dataTypes = rep('continuous', length(dataInput)),
                        combMethods = c("Fisher", "Liptak", "Tippett"),
                        numPerms = 1000,
                        numCores = 1,
                        verbose = FALSE,
                        functionGeneratingIndex = NULL,
                        outcomeName = NULL,
                        allCombinations = FALSE,
                        dataWeights = rep(1, length(dataInput))/length(dataInput),
                        returnPermPvalues = FALSE,
                        ...)
               {standardGeneric("omicsNPC")}
)

setMethod(
    f="omicsNPC",
    signature=signature(dataInput = "list", dataMapping = 'data.frame'),
    definition=function(dataInput, dataMapping, dataTypes, combMethods, numPerms, 
                        numCores, verbose, functionGeneratingIndex, outcomeName, 
                        allCombinations, dataWeights, returnPermPvalues, ...){
        
      # Input check 

      # checking the data in input
      classNames <- sapply(dataInput, class)
      testClass <- sapply(X = classNames, FUN = match.arg, choices = c('ExpressionSet'))
      
      #if no names for the datasets, create fake ones
      if(is.null(names(dataInput))){
        names(dataInput) <- paste('dataset', 1:length(dataInput), sep = "_")
      }

      # data types
      if(class(dataTypes) == 'character'){
        testDataTypes <- sapply(X = dataTypes, FUN = match.arg, choices = c('count', 'continuous')) 
      }
      
      #weigths
      if(sum(dataWeights) != 1){
        stop('dataWeights must be a numeric vector such that sum(dataWeights) == 1')
      }
      
      # Create input for internal omicsNPC 
      
      # extract matrices from ExpressionSet
      dataMatrices <- lapply(dataInput, FUN = exprs)
      if(is.null(names(dataInput))){
        names(dataMatrices) <- paste('Dataset', 1:length(dataInput), sep = '_');
      }else{
        names(dataMatrices) <- names(dataInput);
      }

      # retrieving the whole design for each dataset
      designs <- lapply(X = dataInput, FUN = pData)
      if(is.null(names(dataInput))){
        names(designs) <- paste('Dataset', 1:length(dataInput), sep = '_');
      }else{
        names(designs) <- names(dataInput);
      }

      # functions to analyze data
      if(class(dataTypes) == 'character'){
        functionsAnalyzingData <- vector('list', length(dataTypes));
        for(i in 1:length(dataTypes)){
          if(dataTypes[i] == 'continuous'){
            functionsAnalyzingData[[i]] <- computeAssocContinuousData;
          }else{
            functionsAnalyzingData[[i]] <- computeAssocCountData;
          }
        }
      }else{
        functionsAnalyzingData <- dataTypes;
      }
      
      # function to permute data
      if(is.null(functionGeneratingIndex)){
        functionGeneratingIndex <- generate_iid_data_index 
      }
      
      # run omicsNPC 
      output <- omicsNPC_internal(dataMatrices = dataMatrices, 
                                  designs = designs,
                                  dataMapping = dataMapping,
                                  combMethods = combMethods,
                                  functionsAnalyzingData = functionsAnalyzingData,
                                  functionGeneratingIndex = functionGeneratingIndex,
                                  outcomeName = outcomeName,
                                  numPerms = numPerms,
                                  numCores = numCores,
                                  verbose = verbose,
                                  allCombinations = allCombinations,
                                  dataWeights = dataWeights,
                                  returnPermPvalues = returnPermPvalues,
                                  ...)

      #returning the results
      return(output)

    }
)

setMethod(
  f="omicsNPC",
  signature=signature(dataInput = "list", dataMapping = 'missing'),
  definition=function(dataInput, dataMapping, dataTypes, combMethods, numPerms, numCores, verbose, functionGeneratingIndex, outcomeName, returnPermPvalues, ...){
    
    # Input check 
    
    # checking the data in input
    if(class(dataInput) == 'list'){
      classNames <- sapply(dataInput, class)
      testClass <- sapply(X = classNames, FUN = match.arg, choices = c('ExpressionSet'))
    }
    
    #if no names for the datasets, create fake ones
    if(is.null(names(dataInput))){
      names(dataInput) <- paste('dataset', 1:length(dataInput), sep = "_")
    }
    
    # data types
    if(class(dataTypes) == 'character'){
      testDataTypes <- sapply(X = dataTypes, FUN = match.arg, choices = c('count', 'continuous')) 
    }
    
    #weigths
    if(sum(dataWeights) != 1){
      stop('dataWeights must be a numeric vector such that sum(dataWeights) == 1')
    }
    
    #creating the data mapping based on the rownames 
    #(we assume the all datasets have the same rownames encoding, e.g., probeset ids of the same platform)
    mappings <- vector('list', length(dataInput));
    for(i in 1:length(dataInput)){
      mappings[[i]] <- data.frame(id = rownames(dataInput[[i]]),
                                  reference = rownames(dataInput[[i]]));
    }
    dataMapping <- combiningMappings(mappings = mappings, reference = 'reference', retainAll = TRUE);
    
    # Create input for internal omicsNPC 
    
    # extract matrices from ExpressionSet
    dataMatrices <- lapply(dataInput, FUN = exprs)
    if(is.null(names(dataInput))){
      names(dataMatrices) <- paste('Dataset', 1:length(dataInput), sep = '_');
    }else{
      names(dataMatrices) <- names(dataInput);
    }
    
    # retrieving the whole design for each dataset
    designs <- lapply(X = dataInput, FUN = pData)
    if(is.null(names(dataInput))){
      names(designs) <- paste('Dataset', 1:length(dataInput), sep = '_');
    }else{
      names(designs) <- names(dataInput);
    }
    
    # functions to analyze data
    if(class(dataTypes) == 'character'){
      functionsAnalyzingData <- vector('list', length(dataTypes));
      for(i in 1:length(dataTypes)){
        if(dataTypes[i] == 'continuous'){
          functionsAnalyzingData[[i]] <- computeAssocContinuousData;
        }else{
          functionsAnalyzingData[[i]] <- computeAssocCountData;
        }
      }
    }else{
      functionsAnalyzingData <- dataTypes;
    }
    
    # function to permute data
    if(is.null(functionGeneratingIndex)){
      functionGeneratingIndex <- generate_iid_data_index 
    }
    
    # run omicsNPC 
    output <- omicsNPC_internal(dataMatrices = dataMatrices, 
                                designs = designs,
                                dataMapping = dataMapping,
                                combMethods = combMethods,
                                functionsAnalyzingData = functionsAnalyzingData,
                                functionGeneratingIndex = functionGeneratingIndex,
                                outcomeName = outcomeName,
                                numPerms = numPerms,
                                numCores = numCores,
                                verbose = verbose,
                                allCombinations = allCombinations,
                                dataWeights = dataWeights,
                                returnPermPvalues = returnPermPvalues,
                                ...)
    
    #returning the results
    return(output)
    
  }
)