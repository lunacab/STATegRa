#' @export
#' @import affy
#' @import edgeR
#' @import limma
#' @import foreach
#' @title omicsPC, applying the Parametric Combination (PC) on omics datasets
#' @aliases omicsPC,list,data.frame-method omicsPC,list,missing-method
#' @description
#' This function applies parametric combination methods on the integrative analysis 
#' of different omics data modalities.
#' It retrieves genes associated to a given outcome, taking into account all omics data.
#' First, each datatype is analyzed independently using the appropriate method.
#' omicsNPC analyses continuous data (microarray) using limma, while count data (e.g., RNAseq) 
#' are first preprocessed with using the "voom" function. The user can also specify their own 
#' function for computing deregulation / association 
#' The p-values from the single dataset analysis are then combined employing several combining functions, 
#' see references for further information.
#' @usage omicsPC(dataInput, dataMapping, dataTypes = rep('continuous', length(dataInput)), 
#'                combMethods = c("Fisher", "Liptak", "Tippett", "Benjamini", "Simes", "Sidak"), 
#'                verbose = FALSE, outcomeName = NULL, allCombinations = FALSE, ...)
#'
#' @param dataInput List of ExpressionSet objects, one for each data modality. 
#' @param dataMapping A data frame describing how to map measurements across datasets. See details for more information.
#' @param dataTypes Character vector with possible values: 'continuous', 'count'. Alternatively, a list of functions for assessing deregulation / association with an outcome
#' @param combMethods Character vector with possible values: "Fisher", "Liptak", "Tippett", "Benjamini", "Simes", "Sidak". If more than one is specified, all will be used.
#' @param verbose Logical, if set to TRUE omicsNPC prints out the step that it performs
#' @param outcomeName Name of the outcome of interest / experimental factor, as reported in the design matrices. If NULL, the last column of the design matrices is assumed to be the outcome of interest.
#' @param allCombinations Logical, if TRUE all combinations of omics datasets are considered
#' @param ... Additional arguments to be passed to the user-defined functions
#'
#' @return A list containing:
#' stats0 Partial deregulation / association statistics 
#' pvalues0 The partial p-values computed on each dataset
#' pvaluesPC The p-values computed through PC.
#'
#' @author Vincenzo Lagani
#'
#' @references
#' Peter H. Westfall. Combining P values. Encyclopedia of Biostatistics (2005).
#' Nestoras Karathanasis, Ioannis Tsamardinos and Vincenzo Lagani. omicsNPC: applying the Non-Parametric Combination 
#' methodology to the integrative analysis of heterogeneous omics data. PlosONE 11(11): e0165545. doi:10.1371/journal.pone.0165545
#'
#' @examples
#' # Load the data
#' data("TCGA_BRCA_Batch_93")
#' # Setting dataTypes, the first two ExpressionSets include RNAseq data,
#' # the third ExpressionSet includes Microarray data.
#' dataTypes <- c("count", "count", "continuous")
#' # Setting methods to combine pvalues. Notice that Benjamini and Simes are actually the same method
#' combMethods = c("Fisher", "Liptak", "Tippett", "Benjamini", "Simes", "Sidak")
#' # Setting omicsNPC to print out the steps that it performs.
#' verbose = TRUE
#' # Run omicsNPC analysis.
#' # The output list contains a data.frame of p-values, namely 'pvaluesPC', where each row 
#' # corresponds to a gene, and each column corresponds to a method used in the analysis.
#' 
#' \dontrun{out <- omicsPC(dataInput = TCGA_BRCA_Data, dataTypes = dataTypes,
#'                             combMethods = combMethods, verbose = verbose)}

setGeneric(name="omicsPC",
           def=function(dataInput,
                        dataMapping,
                        dataTypes = rep('continuous', length(dataInput)),
                        combMethods = c("Fisher", "Liptak", "Tippett", "Benjamini", "Simes", "Sidak"),
                        verbose = FALSE,
                        outcomeName = NULL,
                        allCombinations = FALSE,
                        ...)
               {standardGeneric("omicsPC")}
)

setMethod(
    f="omicsPC",
    signature=signature(dataInput = "list", dataMapping = 'data.frame'),
    definition=function(dataInput, dataMapping, dataTypes, combMethods, 
                        verbose, outcomeName, allCombinations, ...){
        
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
      }else{
        if(class(dataTypes) != 'list'){
          stop('dataTypes must be either a character vector or a list')
        }else{
          if(length(dataTypes) == 0 || class(dataTypes[[1]]) != 'function'){
            stop('if dataTypes is a list it must contain functions')
          }
        }
      }
      
      # Create input for internal omicsPC 
      
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
      
      # run omicsPC 
      output <- omicsPC_internal(dataMatrices = dataMatrices, 
                                  designs = designs,
                                  dataMapping = dataMapping,
                                  combMethods = combMethods,
                                  functionsAnalyzingData = functionsAnalyzingData,
                                  outcomeName = outcomeName,
                                  verbose = verbose,
                                  allCombinations = allCombinations,
                                  ...)

      #returning the results
      return(output)

    }
)

setMethod(
  f="omicsPC",
  signature=signature(dataInput = "list", dataMapping = 'missing'),
  definition=function(dataInput, dataMapping, dataTypes, combMethods, verbose, outcomeName, ...){
    
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
    }else{
      if(class(dataTypes) != 'list'){
        stop('dataTypes must be either a character vector or a list')
      }else{
        if(length(dataTypes) == 0 || class(dataTypes[[1]]) != 'function'){
          stop('if dataTypes is a list it must contain functions')
        }
      }
    }
    
    #creating the data mapping based on the rownames 
    #(we assume the all datasets have the same rownames encoding, e.g., probeset ids of the same platform)
    mappings <- vector('list', length(dataInput));
    for(i in 1:length(dataInput)){
      mappings[[i]] <- data.frame(id = rownames(dataInput[[i]]),
                                  reference = rownames(dataInput[[i]]));
    }
    dataMapping <- combiningMappings(mappings = mappings, reference = 'reference', retainAll = TRUE);
    
    # Create input for internal omicsPC 
    
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
    
    # run omicsPC 
    output <- omicsPC_internal(dataMatrices = dataMatrices, 
                                designs = designs,
                                dataMapping = dataMapping,
                                combMethods = combMethods,
                                functionsAnalyzingData = functionsAnalyzingData,
                                outcomeName = outcomeName,
                                verbose = verbose,
                                allCombinations = allCombinations,
                                ...)
    
    #returning the results
    return(output)
    
  }
)