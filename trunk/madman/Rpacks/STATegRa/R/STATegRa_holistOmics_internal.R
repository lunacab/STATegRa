# # holistOmics_internal
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

# Documentation to add in next release
# @title Central function of holistOmics
# @description
# This function applies the NonParametric Combination methodology on the integrative analysis
# of different omics data modalities. It retrieves diffential expressed genes between
# control and treatment, taking into account all omics data. First, each datatype
# is analyzed independently using the appropriate method.HolistOmics analyses static (one time point)
# RNAseq, using Voom+limma method, and static (one time point) Microarray, using limma.
# Then, the resulting p-values are combined employing Fisher, Liptak and Tippett combining functions.
# Tippett function returns findings which are supported by at least one omics modality.
# Liptak function returns findings which are supportd by most modalities.
# Fisher function has an intermediate behavior between those of Tippett and Liptak.
#
# @usage holistOmics_internal(Data, designs, comb.method, types.of.data, functions.Analyze.Data,
#                          arguments.Analyze.Data, functions.generate.index, functions.permute.Data,
#                          numPerm = 1000, numCores = 1, verbose = FALSE,outputfolder = NULL,
#                          save_intermediate = FALSE, extention = NULL)
#
# @param Data List of matrixes, one for each data modality.
# @param designs List of matrices. The first matrix, is used as the design of the first data matrix,
# the second matrix is used as the design of the second data matrix, etc.
# Depending on the data, designs could differ.
# @param comb.method Character vector with possible values: 'Fisher', 'Liptak', 'Tippett',
# if more than one is specified, all will be used.
# @param types.of.data Character vector labelling types of data. They do not influence analysis process.
# @param functions.Analyze.Data List containing functions to analyse data. First function
# analyse first element in Data list etc
# @param arguments.Analyze.Data List, possible arguments of functions. Each element in
# arguments is also a list
# @param functions.generate.index Function to generate index during permutation. The SAME index
# will be used for all data types
# @param functions.permute.Data List containing functions to permute data. First function
# permutes first element in Data list etc
# @param numPerm Number of permutations
# @param numCores Number of CPU cores to use
# @param verbose Logical, if set to TRUE holistOmics prints out the step that it performs
# @param outputfolder Character vector with folder name to save results, the folder to be created by the user.
# @param save_intermediate Logical indicating if the intermediate output will be saved.
# If TRUE a cache folder is created. If output is NULL cache is created in working directory.
# If output is not NULL cashe is created in output folder
# @param extention Characted vector. The extension will be added in all
# created files during a run
#
# @return A data.frame with pvalues
#
# @author Nestoras Karathanasis
#
# @references
# Pesarin, Fortunato, and Luigi Salmaso. Permutation tests for complex data: theory, applications and software. John Wiley & Sons, 2010.
#
# @examples
# # load data
# data(TCGA_BRCA_Batch_93)
# designs
# designs <- lapply(X = 1:3, FUN = function(i){out <- c(rep(0, 16), rep(1, 16))})
# Setting methods to combine individual pvalues
# comb.method <- c("Fisher", "Liptak", "Tippett")
# Setting functions to analyse data. Two functions are available.
# One which analyze RNAseq data, calcDifferenceRNAseqVoomLimma
# One which analyze Microarray data, calcDifferenceMicroarray.
# functions.Analyze.Data <- list(calcDifferenceRNAseqVoomLimma,
#                               calcDifferenceRNAseqVoomLimma,
#                               calcDifferenceMicroarray)
# Setting arguments of the above functions. Only argument verbose is allowed.
# Argument should be provided as lists.
# arguments.Analyze.Data <- list(list("verbose" = T), list("verbose" = T), list("verbose" = T))
# Setting function to generate index during permutation
# functions.generate.index <- generate_static_data_index
# Setting functions to permute data. One function is available.
# functions.permute.Data <- list(permute_static_data, permute_static_data, permute_static_data)
# Setting number of permutations to performed
# numPerm = 1000
# Setting number of cores to be used
# numCores = 1
# Setting verbose, if TRUE, the algorithm prints out the steps it performs.
# verbose = TRUE
# # The rest of the arguments are mainly for development purposes
# Setting folder to save results.
# outputfolder = NULL
# Setting save_intermediate, if TRUE the algorithm saves the intermediate results to the folder cache.
# If outputfolder is NULL cache is created to the working directory.
# If outputfolder is not NULL cache is created to the outputfolder directory.
# save_intermediate = FALSE
# Setting extension. Extention is a character vector which will be added in the end of all filenames
# created during a run.
# extention = NULL
# Run holistOmics analysis using holistOmics internal
# The output is a data.frame of p-values.
# Each row corresponds to a gene name. Each column corresponds to a method
# used in the analysis.
# out <- holistOmics_internal(Data, designs, comb.method, types.of.data,
#                                     functions.Analyze.Data, arguments.Analyze.Data,
#                                     functions.generate.index, functions.permute.Data,
#                                     numPerm = 1000, numCores = 1,
#                                     verbose = FALSE,
#                                     outputfolder = NULL,
#                                     save_intermediate = FALSE,
#                                     extention = NULL)}
#
holistOmics_internal <- function(Data, designs, comb.method, types.of.data,
                                 functions.Analyze.Data, arguments.Analyze.Data,
                                 functions.generate.index, functions.permute.Data,
                                 numPerm = 1000, numCores = 1,
                                 verbose = FALSE,
                                 outputfolder = NULL,
                                 save_intermediate = FALSE,
                                 extention = NULL){

    # folder to save intermediate results
    if(save_intermediate){
        if(!is.null(outputfolder)){
            intermediate_results_folder <- paste(outputfolder, "cache", sep = "/")
        }else{
            intermediate_results_folder <- "cache"
        }
        dir.create(intermediate_results_folder)
        message(paste("Intermediate resutls will be saved in ...", intermediate_results_folder))
    }

    # FUNCTION HANDLERS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # calculate DE
    calculateDE <- function(Data, designs, functions, arguments){
        # function
        executeFor <- function(cDE, Datas, fun, d, args){
            out <- fun[[cDE]](Datas[[cDE]], d[[cDE]],args[[cDE]])
            return(out)
        }
        results <- sapply(X = 1:length(Data), FUN = executeFor, Data, functions, designs, arguments)
        return(results)
    }

    # Permute Data
    permuteDATA <- function(designs, f_generate_index, f_perm){
        # function
        executeFor <- function(pD, edesigns, fun, args){
            out <- fun[[pD]](edesigns[[pD]], args[[pD]])
            return(out)
        }
        # use the first design table to produce the indexing
        index <- f_generate_index(classIn = designs[[1]])
        indexes <- rep(list(index), length(designs))
        results <- lapply(X = 1:length(designs), FUN = executeFor, designs, f_perm, indexes)
        return(results)
    }

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # ............................. STEP.1 Calculate initial pvalues .............................................
    if(verbose){
        message("Calculate initial pvalues on data")
    }
    pvalues0 <- calculateDE(Data, designs = designs, functions = functions.Analyze.Data, arguments =  arguments.Analyze.Data)
    if(save_intermediate){
        filename <- sprintf("%s/pvalues0%s.RData", intermediate_results_folder, extention)
        save(pvalues0, file = filename)
    }

    #..................... STEP.2 Building NULL distributions by permuting data ..........................................
    if(verbose){
        message("Building NULL distributions by permuting data")
    }

    # initialize lists of variables tables
    pvaluesPERM_program <- vector("list",numPerm+1)

    # fill pvalues table
    pvaluesPERM_program[[1]] <- pvalues0

    # set function for permutation
    executePermutation <- function(i){
        # permute designs
        permDesign <- permuteDATA(designs = designs, f_generate_index = functions.generate.index, f_perm = functions.permute.Data)
        # recalculate pvalues
        DEresultsi <- calculateDE(Data = Data, designs = permDesign, functions = functions.Analyze.Data, arguments = arguments.Analyze.Data)
        return(DEresultsi)
    }

    # run in parallel #
    i <- NULL;
    if (numCores != 1){
        # check if required packages are available
        if(requireNamespace("doSNOW", quietly = TRUE)){
            #number of cores to use
            nCores <- numCores

            #creating the cluster
            cl <- parallel::makeCluster(nCores)

            #registering the cluster
            doSNOW::registerDoSNOW(cl)

            # calculate DE expression
            pvalues_PERMi <- foreach(i = 1:numPerm) %dopar% executePermutation(i)

            #stopping the cluster
            parallel::stopCluster(cl)
        }else{
            message("Install doSNOW to run holistOmics in parallel...")
            message("holistOmics will run using ONE core")
            pvalues_PERMi <- foreach(i = 1:numPerm) %do% executePermutation(i)
        }
    }else {
        pvalues_PERMi <- foreach(i = 1:numPerm) %do% executePermutation(i)

    }
    pvaluesPERM_program[2:(numPerm+1)] <- pvalues_PERMi

    if(save_intermediate){
        filename <- sprintf("%s/pvaluesPERM_program%s.RData", intermediate_results_folder, extention)
        save(pvaluesPERM_program, file = filename)
    }

    #...................... STEP.3 Calculate p-values based on NULL distributions ............................
    if(verbose){
        message("Calculate pseudo p-values based on NULL distributions...")
    }
    # pvalues for each gene together
    pvaluesPerGeneProgram <- builtpvalueCalculatorInput(sample_list = pvaluesPERM_program)
    # calculate stats... 1 - pvalues
    stats <- lapply(X = pvaluesPerGeneProgram, FUN = p2t)
    # recalculate pvalues based on stats
    pvaluesPerGenePerm <- lapply(X = stats, FUN = t2p)
    pvaluesPERM_PERM <- builtpvalueCalculatorInput(sample_list = pvaluesPerGenePerm)

    #................... STEP .4 Calculate p-values for each gene .................................................
    if(verbose){
        message("p-values calculation...")
    }
    # NPC using Pesarin
    executeFor_NPC <- function(m){
        pvalues.NPC.method <- sapply(X = pvaluesPerGenePerm, FUN = calc.NPC.pvalues.Pesarin, comb.method[m])
        return(pvalues.NPC.method)
    }
    pvalues.NPC <- sapply(X = 1:length(comb.method), FUN = executeFor_NPC)
    colnames(pvalues.NPC) <- comb.method

    # separate pvalues for each datatype
    pvalues.reg <- pvaluesPERM_PERM[[1]]
    colnames(pvalues.reg) <- paste(types.of.data, "Perm", sep = "")
    # other programs pvalues
    colnames(pvalues0) <- paste(types.of.data, "Programs", sep = "")

    # all pvalues
    pvaluesALL <- data.frame(pvalues0, pvalues.reg, pvalues.NPC)
    # gene names as row names
    row.names(pvaluesALL) <- rownames(Data[[1]])

    # save pvalues
    if(!is.null(outputfolder)){
        filename <- sprintf("%s/pvaluesALL%s.RData", outputfolder, extention)
        save(pvaluesALL, file = filename)}

    return(pvaluesALL)

}



# functions to run holistOmics on STATIC data --------------------------------

# CALCULATE DIFFERENTIAL EXPRESSION ###################################
calcDifferenceMicroarray <- function(exprData, pheno, args){
    # extract args
    verbose <- args$verbose
    if(verbose){
        print("Calculate differential expression in microarray using LIMMA....")
    }
    #load the necessary libraries
    #library(limma)
    #matrix with two columns, 1st(Intercept) 1 + 2nd "phenoType" vals
    design <- model.matrix( ~ 1 + pheno)
    colnames(design) <- c("Affy.Intercept", "Affy.T")
    fit <- lmFit(exprData, design=design)
    fit2 <- eBayes(fit)
    dec.a <- decideTests(fit2, method="separate",adjust.method="fdr",p.value=0.05, lfc=0)
    results <- topTable(fit2, 2, number=Inf, sort.by='none')
    # add the regulation of each gene in the end of results
    results <- data.frame(results, regulation = dec.a[,2])
    pvalues <- results$P.Value
    return(pvalues)
}

calcDifferenceRNAseqVoomLimma <- function(dataRNAseq, design, args){
    # extract args
    verbose <- args$verbose
    if(verbose){
        print("Calculate differential expression of RNAseq using Voom + limma....")
    }
    nf = calcNormFactors(dataRNAseq, method = "TMM")
    voom.data = voom(dataRNAseq, design = model.matrix(~factor(design)),
                     lib.size = colSums(dataRNAseq) * nf)
    voom.data$genes = rownames(dataRNAseq)
    voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(design)))
    voom.fitbayes = eBayes(voom.fitlimma)
    pvalues = voom.fitbayes$p.value[, 2]
    stats = voom.fitbayes$t[,2]
    voom.adjpvalues = p.adjust(pvalues, method = "BH")
    output <- pvalues
    return(output)
}

# PERMUTE DATA #######################################################
permute_static_data <- function(classIn, index){
    conditionOut <- classIn[index]
    return(conditionOut)
}

generate_static_data_index <- function(classIn){
    index <- sample(x = 1:length(classIn), size = length(classIn),replace = FALSE)
    return(index)
}

# FUNCTIONS TO RUN HOLISTOMICS CANNOT BE CHANGED #########################
builtpvalueCalculatorInput <- function(sample_list){
    bind.ith.rows <- function(i) do.call(rbind, lapply(sample_list, "[", i, TRUE))
    nr <- nrow(sample_list[[1]])
    out <- lapply(1:nr, bind.ith.rows)

    return(out)
}

# CALCULATE P-VALUES ####################################################
# calculate pvalues using NPC
calc.NPC.pvalues.Pesarin <- function(table_Of_pvalues, method){
    if(method == "Fisher"){
        stats <- apply(table_Of_pvalues,1,function(x){-2*log(prod(x))})
        pvalue_temp <- t2p(stats)
        pvalue <- pvalue_temp[1]
    }else if(method == "Liptak"){
        stats <- apply(table_Of_pvalues,1,function(x){sum(qnorm(1-x))})
        pvalue_temp <- t2p(stats)
        pvalue <- pvalue_temp[1]
    }else if(method == "Tippett"){
        stats <- apply(table_Of_pvalues,1,min)
        pvalue_temp <- t2p(1/stats)
        pvalue <- pvalue_temp[1]
    }
    return(pvalue)
}

# statistics to pvalues
t2p <- function(table){
    # Pesarin
    if(is.null(dim(table))){table<-array(table,dim=c(length(table),1))}
    oth<-seq(1:length(dim(table)))[-1]
    B<-dim(table)[1]-1
    p<-dim(table)[2]
    if(length(dim(table))==3){C<-dim(table)[3]}
    rango<-function(x){
        r=1-rank(x[-1],ties.method="min")/B+1/B
        return(c(mean(x[-1]>=x[1]),r))
    }
    P=apply(table,oth,rango)
    return(P)
}

# pvalues to statistics
p2t <- function(pvalues){
    t <- 1 - pvalues
}
