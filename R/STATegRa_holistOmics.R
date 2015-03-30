# holistOmics -------------------------------------------------------
## INPUT:
##    ######...: Set of expression sets one for each omic dataset
##    dataInput: (list) List of expression sets one for each omic dataset
##    dataTypes: (character) string vector with possible values: 'RNA-seq', 'microarray'
##    design: (numeric) the outcome in numeric form
##    comb.method = c("Fisher", "Liptak", "Tippett"): (character) Method(s) for combining p-values
##    numPerm = 1000: (numeric) number of permutations to perform
##    numCores = 1: (numeric) number of CPU cores to use
## OUTPUT:
##    An object of class 'holistOmics'

#' @export
#' @import affy
#' @import edgeR
#' @import limma
#' @import foreach
#' @title HolistOmics an application of NPC on omics datasets
#' @aliases holistOmics,list,character-method
#' @description
#' This function applies the NonParametric Combination methodology on the integrative analysis of
#' different omics data modalities.
#' It retrieves diffential expressed genes between control and treatment, taking into account all omics data.
#' First, each datatype is analyzed independently using the appropriate method.
#' HolistOmics analyses static (one time point) RNAseq, using Voom+limma method, and static (one time point)
#' Microarray, using limma.
#' Then, the resulting p-values are combined employing Fisher, Liptak and Tippett combining functions.
#' Tippett function returns findings which are supported by at least one omics modality.
#' Liptak function returns findings which are supportd by most modalities.
#' Fisher function has an intermediate behavior between those of Tippett and Liptak.
#'
#' @usage holistOmics(dataInput, dataTypes, comb.method = c("Fisher", "Liptak", "Tippett"),
#'                numPerm = 1000, numCores = 1, verbose = FALSE)
#'
#' @param dataInput List of ExpressionSet objects, one for each data modality.
#' @param dataTypes Character vector with possible values: 'RNA-seq', 'microarray'
#' @param comb.method Character vector with possible values: 'Fisher', 'Liptak',
#' 'Tippett', if more than one is specified, all will be used.
#' @param numPerm Number of permutations
#' @param numCores Number of CPU cores to use
#' @param verbose Logical, if set to TRUE holistOmics prints out the step that it performs
#'
#' @return A data.frame
#'
#' @author Nestoras Karathanasis
#'
#' @references
#' Pesarin, Fortunato, and Luigi Salmaso. Permutation tests for complex data: theory, applications and software.
#' John Wiley & Sons, 2010.
#'
#' @examples
#' # Load the data
#' data("TCGA_BRCA_Batch_93")
#' # Setting dataTypes, the first two ExpressionSets include RNAseq data,
#' # the third ExpressionSet includes Microarray data.
#' dataTypes <- c("RNAseq", "RNAseq", "Microarray")
#' # Setting methods to combine pvalues
#' comb.method = c("Fisher", "Liptak", "Tippett")
#' # Setting number of permutations
#' numPerm = 1000
#' # Setting number of cores
#' numCores = 1
#' # Setting holistOmics to print out the steps that it performs.
#' verbose = TRUE
#' # Run holistOmics analysis.
#' # The output is a data.frame of p-values.
#' # Each row corresponds to a gene name. Each column corresponds to a method
#' # used in the analysis.
#' \dontrun{out <- holistOmics(dataInput = TCGA_BRCA_Data, dataTypes = dataTypes,
#'                             comb.method = comb.method, numPerm = numPerm,
#'                             numCores = numCores, verbose = verbose)}

setGeneric(name="holistOmics",
           def=function(dataInput, dataTypes,
                        comb.method = c("Fisher", "Liptak", "Tippett"),
                        numPerm = 1000, numCores = 1, verbose = FALSE)
               {standardGeneric("holistOmics")}
           )

setMethod(
    f="holistOmics",
    signature=signature(
        dataInput = "list",
        dataTypes = "character"),
    definition=function(dataInput, dataTypes, comb.method, numPerm, numCores, verbose){
        # Input check # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # class names
        classNames <- sapply(dataInput, class)
        testClass <- sapply(X = classNames, FUN = match.arg, choices = c('ExpressionSet'))

        # data types
        testDataTypes <- sapply(X = dataTypes, FUN = match.arg, choices = c('RNAseq', 'Microarray'))

        # Create input for internal holistOmics # # # # # # # # # # # # # # # # # # # # # # # #
        # extract matrices from ExpressionSet
        dataMatrices <- lapply(dataInput, FUN = exprs)
        # replicate designs
        designs <- lapply(X = dataInput, FUN = function(df){out <- pData(df)[[1]]})

        #function creating input
        createListsOfArguments <- function(datatype, argForRNA, argForMicr){
            if(datatype == "RNAseq"){
                output <- argForRNA
            }else{
                output <- argForMicr
            }
            return(output)
        }

        # functions to analyze data
        funAnalysis <- lapply(X = dataTypes, FUN = createListsOfArguments,
                              calcDifferenceRNAseqVoomLimma, calcDifferenceMicroarray)
        # args of respective functions
        argsAnalysis <- lapply(X = dataTypes, FUN = createListsOfArguments,
                               list("verbose" = FALSE), list("verbose" = FALSE))

        # functions to permute data
        funGenInd <- generate_static_data_index
        funPerm <- lapply(X = dataTypes, FUN = createListsOfArguments,
                          permute_static_data, permute_static_data)

        # run holistOmics # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        output <- holistOmics_internal(Data = dataMatrices, designs = designs,
                                       comb.method = comb.method,
                                       types.of.data = dataTypes,
                                       functions.Analyze.Data = funAnalysis,
                                       arguments.Analyze.Data = argsAnalysis,
                                       functions.generate.index = funGenInd,
                                       functions.permute.Data = funPerm,
                                       numPerm = numPerm,
                                       numCores = numCores,
                                       verbose = verbose)

        return(output)

    }
)
