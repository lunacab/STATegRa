####### MAIN FUNCTION #######

# caClass -----------------------------------------------------------------
## SLOTS:
##   InitialData: (list) List of ExpressionSets, one for each set of omics data
##   preprocessing: (vector) Vector indicating which prepropressing (scale/center/weight) is applied to initial data
##   preproData: (list) List of matrices, one for each set of preprocessed omics data
##   caMethod: (character) Method applied for Components Analysis
##   commonComps: (numeric) Number of common components
##   distComps: (vector) Number of distinctive components for each omics data set
##   scores: (list) List of matrices of common and distinctive scores
##   loadings: (list) List of matrices of common and distintive loadings
##   VAF: (list) list of matrices indicating VAF of each component in each block of data and in common
##   others: (list) List of other information that the selected Component Analysis Method should give.

#' @title caClass
#' @aliases caClass-class
#' @description Stores the results of any of the omicsPCA analyses.
#' @author Patricia Sebastian Leon
#' @slot InitialData List of ExpressionSets, one for each set of omics data
#' @slot Names Character vector giving names for the input data
#' @slot preprocessing Character vector describing the preprocessing applied to the data
#' @slot preproData List of matrices containing data after preprocessing
#' @slot caMethod Character giving the component analysis method name
#' @slot commonComps Numeric giving the number of common components
#' @slot distComps Numeric vector giving the number of distinctive components for each block
#' @slot scores List of matrices of common and distinctive scores
#' @slot loadings List of matrices of common and distinctive loadings
#' @slot VAF List of matrices indicating VAF (Variability Explained For) for each component in each block of data
#' @slot others List containing other miscellaneous information specific to different SCA methods
setClass(
    Class="caClass",
    slots=list(
        InitialData="list",
        Names="character",
        preprocessing="character",
        preproData="list",
        caMethod="character",
        commonComps="numeric",
        distComps="vector",
        scores="list",
        loadings="list",
        VAF="list",
        others="list"),
    prototype=list(
        InitialData=NULL,
        Names=NULL,
        preprocessing=c("none"),
        preproData=NULL,
        caMethod="none",
        commonComps=0,
        distComps=c(0,0),
        scores=NULL,
        loadings=NULL,
        VAF=NULL,
        others=NULL)
)

####### caClass MEMBER ACCESS METHODS #######

#' @export
#' @title Retrieve initial data from caClass objects
#' @aliases getInitialData,caClass-method
#' @description 
#' Generic function to retrieve the initial data used for by \code{\link{omicsCompAnalysis}} from a \code{\link{caClass-class}} object
#' @usage getInitialData(x, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param block Character indicating the block of data to be returned. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return The requested data block or blocks
#' 
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getInitialData(res)
#' getInitialData(res, block="expr")
setGeneric(
    name="getInitialData",
    def=function(x, block=NULL) { standardGeneric("getInitialData") }
)

setMethod(
    f="getInitialData",
    signature=signature(x="caClass"),
    definition=function(x, block=NULL){
        if (class(x) != "caClass") {
            stop("x is not a caclass object")
        } else {
            if (is.null(block)) {
                return(x@InitialData)
            } else {
                if (block == "1" | block == x@Names[1]) {
                    return(x@InitialData[[1]])
                } else if (block == "2" | block == x@Names[2]) {
                    return(x@InitialData[[2]])
                } else {
                    warning(paste(block, "is not an allowed option"))
                    return(x@InitialData)
                }
            }
        }
    }
)




#' @export
#' @title Retrieve information about preprocessing
#' @aliases getPreprocessing,caClass-method
#' @description 
#' Generic function to retrieve information about the preprocessing done by \code{\link{omicsCompAnalysis}} on a \code{\link{caClass-class}} object.
#' @usage getPreprocessing(x, process=FALSE, preproData=FALSE, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param process Logical indicating whether to return information about the processing done.
#' @param preproData Logical indicating whether to return the pre-processed data matrices.
#' @param block Character indicating the block of data to be returned. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return If both \code{process} and \code{preproData} are specified, a list containing (otherwise the individual item):
#' \describe{
#'      \item{process}{Character indicating the processing done}
#'      \item{preproData}{Matrix (or list of matrices, depending on \code{block}) containing pre-processed data}
#'  }
#'  
#' 
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getPreprocessing(res, process=TRUE)
#' getPreprocessing(res, preproData=TRUE, block="1")
setGeneric(
    name="getPreprocessing",
    def=function(x, process=FALSE, preproData=FALSE, block=NULL) {
        standardGeneric("getPreprocessing")
    }
)

setMethod(
    f="getPreprocessing",
    signature=signature(x="caClass"),
    definition=function(x, process, preproData, block) {
        if (class(x) != "caClass"){
            stop("x is not a caclass object")
        } else {
            if (process & preproData) {
                if (is.null(block)) {
                    res <- list(process=x@preprocessing,
                                preproData=x@preproData)
                } else {
                    if (block == "1" | block == x@Names[1]) {
                        res <- list(process=x@preprocessing,
                                    preproData=x@preproData[[1]])
                    } else if (block == "2" | block == x@Names[2]){
                        res <- list(process=x@preprocessing,
                                    preproData=x@preproData[[2]])
                    } else {
                        warning(paste(block, "is not an allowed option"))
                        res <- list(process=x@preprocessing,
                                    preproData=x@preproData)
                    }
                }
            } else if (process) {
                res <- x@preprocessing
            } else if (preproData) {
                if (is.null(block)) {
                    res <- x@preproData
                } else {
                    if (block == "1" | block == x@Names[1]) {
                        res <- x@preproData[[1]]
                    } else if (block == "2" | block == x@Names[2]) {
                        res <- x@preproData[[2]]
                    } else {
                        warning(paste(block, "is not an allowed option"))
                        res <- x@preproData
                    }
                }
            } else {
                res <- list(process=x@preprocessing, preproData=x@preproData)
            }
        }
        return(res)
    }
)



#' @export
#' @title Retrieve information about component analysis method
#' @aliases getMethodInfo,caClass-method
#' @description 
#' Generic function to retrieve information about the method used by \code{\link{omicsCompAnalysis}} on a \code{\link{caClass-class}} object.
#' @usage getMethodInfo(x, method=FALSE, comps=NULL, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param method Logical indicating whether to return the method name.
#' @param comps Character indicating which component number to return ("common", "distinctive" or "all")
#' @param block Character indicating the block of data for which the component count will be given. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return A list containing the requested information.
#'  
#' 
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getMethodInfo(res)
#' getMethodInfo(res, method=TRUE)
#' getMethodInfo(res, comps="all", block="expr")
setGeneric(
    name="getMethodInfo",
    def=function(x, method=FALSE, comps=NULL, block=NULL) {
        standardGeneric("getMethodInfo")
    }
)

setMethod(
    f="getMethodInfo",
    signature=signature(x="caClass"),
    definition=function(x, method, comps, block) {
        if (class(x) != "caClass"){
            stop("x is not a caclass object")
        } else {
            if (method & !is.null(comps)) {
                comps <- match.arg(comps, choices=c("common",
                                                    "distinctive",
                                                    "all"))
                if (comps == "common") {
                    res <- list(method=x@caMethod, comps=x@commonComps)
                } else if (comps == "distinctive") {
                    if (is.null(block)) {
                        res <- list(method=x@caMethod, comps=x@distComps)
                    } else {
                        if (block == "1" | block == x@Names[1]) {
                            res <- list(method=x@caMethod, comps=x@distComps[1])
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- list(method=x@caMethod, comps=x@distComps[2])
                        } else {
                            warning(paste(block, "is not an allowed option"))
                            res <- list(method=x@caMethod, comps=x@distComps)
                        }
                    }
                } else if (comps == "all"){
                    if (is.null(block)) {
                        res <- list(method=x@caMethod,
                                    comps=list(common=x@commonComps,
                                               distinctive=x@distComps))
                    } else {
                        if (block == "1" | block == x@Names[1]) {
                            res <- list(method=x@caMethod,
                                        comps=list(common=x@commonComps,
                                                   distinctive=x@distComps[1]))
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- list(method=x@caMethod,
                                        comps=list(common=x@commonComps,
                                                   distinctive=x@distComps[2]))
                        } else{
                            warning(paste(block, "is not an allowed option"))
                            res <- list(method=x@caMethod,
                                        comps=list(common=x@commonComps,
                                                   distinctive=x@distComps))
                        }
                    }
                }
            } else if (method) {
                res <- x@caMethod
            } else if (!is.null(comps)) {
                comps <- match.arg(comps, choices=c("common",
                                                    "distinctive",
                                                    "all"))
                if (comps == "common") {
                    res <- x@commonComps
                } else if (comps == "distinctive") {
                    if (is.null(block)) {
                        res <- x@distComps
                    } else {
                        if (block == "1" | block == x@Names[1]) {
                            res <- x@distComps[1]
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@distComps[2]
                        } else {
                            warning(paste(block, "is not an allowed option. The output is the list of both blocks"))
                            res <- x@distComps
                        }
                    }
                } else if (comps == "all") {
                    if (is.null(block)) {
                        res <- list(common=x@commonComps,
                                    distinctive=x@distComps)
                    } else {
                        if (block == "1" | block == x@Names[1]) {
                            res <- list(common=x@commonComps,
                                        distinctive=x@distComps[1])
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- list(common=x@commonComps,
                                        distinctive=x@distComps[2])
                        } else {
                            warning(paste(block, "is not an allowed option. The output is the list of both blocks"))
                            res <- list(common=x@commonComps,
                                        distinctive=x@distComps)
                        }
                    }
                }
            } else {
                res <- list(method=x@caMethod, 
                            comps=list(common=x@commonComps,
                                       distinctive=x@distComps))
            }
        }
        return(res)
    }
)


#' @export
#' @title Retrieve component analysis scores
#' @aliases getScores,caClass-method
#' @description 
#' Generic function to retrieve scores (common and distinctive) found by \code{\link{omicsCompAnalysis}} on a \code{\link{caClass-class}} object.
#' @usage getScores(x, part=NULL, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param part Character indicating whether "common" or "distinctive" scores should be displayed
#' @param block Character indicating the block of data for which the scores will be given. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return A list containing the requested information.
#'  
#' 
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getScores(res)
#' getScores(res, part="common")
#' getScores(res, part="distinctive", block="expr")
setGeneric(
    name="getScores",
    def=function(x, part=NULL, block=NULL) { standardGeneric("getScores") }
)

setMethod(
    f="getScores",
    signature=signature(x="caClass"),
    definition=function(x, part, block) {
        if (class(x) != "caClass") {
            stop("x is not a caclass object")
        } else {
            if (is.null(part)) {
                res <- x@scores
            } else {
                if (x@caMethod == "DISCO-SCA" | x@caMethod == "JIVE") {
                    part <- match.arg(part, choices=c("common", "distinctive"))
                    if (part == "common") {
                        res <- x@scores$common
                        if (!is.null(block)) {
                            warning(paste("Using", x@caMethod,
                                          "scores for common part are referring to both blocks together"))
                        }
                    } else if (part == "distinctive") {
                        if (is.null(block)) {
                            res <- x@scores$dist
                        } else if (block == "1" | block == x@Names[1]) {
                            res <- x@scores$dist[[1]]
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@scores[[2]]
                        } else {
                            warning(paste(block, "is not an allowed option."))
                            res <- x@scores
                        }
                    }
                } else if (x@caMethod == "O2PLS") {
                    if (part == "common") {
                        if (is.null(block)) {
                            res <- x@scores$common
                        } else if (block == "1" | block == x@Names[1]) {
                            res <- x@scores$common[[1]]
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@scores$common[[2]]
                        } else {
                            warning(paste(block, "is not an allowed option."))
                            res <- x@scores$common
                        }
                    } else if (part=="distinctive") {
                        if (is.null(block)) {
                            res <- x@scores$dist
                        } else if (block == "1" | block == x@Names[1]) {
                            res <- x@scores$dist[[1]]
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@scores$dist[[2]]
                        } else {
                            warning(paste(block, "is not an allowed option."))
                            res <- x@scores$dist
                        }
                    }
                }
            }
        }
        return(res)
    }
)


#' @export
#' @title Retrieve component analysis loadings
#' @aliases getLoadings,caClass-method
#' @description 
#' Generic function to retrieve loadings (common and distinctive) found by \code{\link{omicsCompAnalysis}} on a \code{\link{caClass-class}} object.
#' @usage getLoadings(x, part=NULL, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param part Character indicating whether "common" or "distinctive" loadings should be displayed
#' @param block Character indicating the block of data for which the loadings will be given. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return A list containing the requested information.
#'  
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getLoadings(res)
#' getLoadings(res, part="common", block="expr")
#' getLoadings(res, part="distinctive", block="expr")
setGeneric(
    name="getLoadings",
    def=function(x, part=NULL, block=NULL) { standardGeneric("getLoadings") }
)

setMethod(
    f="getLoadings",
    signature=signature(x="caClass"),
    definition=function(x, part, block){
        if (class(x) != "caClass") {
            stop("x is not a caclass object")
        } else {
            if (is.null(part)) {
                
            } else {
                part <- match.arg(part, choices=c("common", "distinctive"))
                if (part == "common"){
                    if (is.null(block)) {
                        res <- x@loadings$common
                    } else if (block == "1" | block == x@Names[1]) {
                        res <- x@loadings$common[[1]]
                    } else if (block == "2" | block == x@Names[2]) {
                        res <- x@loadings$common[[2]]
                    } else{
                        warning(paste(block, "is not an allowed option."))
                        res <- x@loadings$common
                    }
                } else if (part == "distinctive") {
                    if (is.null(block)) {
                        res <- x@loadings$dist
                    } else if (block == "1" | block == x@Names[1]) {
                        res <- x@loadings$dist[[1]]
                    } else if (block == "2" | block == x@Names[2]) {
                        res <- x@loadings$dist[[2]]
                    } else {
                        warning(paste(block, "is not an allowed option."))
                        res <- x@loadings$dist
                    }
                }
            }
        }
        return(res)
    }
)


#' @export
#' @title Retrieve information abotut VAF
#' @aliases getVAF,caClass-method
#' @description 
#' Generic function to retrieve VAF found by \code{\link{omicsCompAnalysis}} on a \code{\link{caClass-class}} object.
#' @usage getVAF(x, part=NULL, block=NULL)
#' @param x \code{\link{caClass-class}} object.
#' @param part Character indicating whether "common" or "distinctive" VAF should be displayed
#' @param block Character indicating the block of data for which the VAF will be given. It can be specified by the position of the block ("1" or "2") or the name assigned in the \code{\link{caClass-class}} object. If it is NULL both blocks are displayed.
#' 
#' @return A list containing the requested information.
#'  
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{omicsCompAnalysis}}, \code{\link{caClass-class}}
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA, pDataDescr=c("classname"))
#' # Omics components analysis
#' res <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "mirna"),
#'                          method="DISCOSCA", Rcommon=2, Rspecific=c(2, 2),
#'                          center=TRUE, scale=TRUE, weight=TRUE)
#' getVAF(res)
#' getVAF(res, part="common")
#' getVAF(res, part="distinctive", block="expr")
setGeneric(
    name="getVAF",
    def=function(x, part=NULL, block=NULL) { standardGeneric("getVAF") }
)

setMethod(
    f="getVAF",
    signature=signature(x="caClass"),
    definition=function(x, part, block){
        if (class(x) != "caClass") {
            stop("x is not a caclass object")
        } else {
            if (x@caMethod == "DISCO-SCA") {
                if (is.null(part)) {
                    res <- x@VAF
                } else {
                    part <- match.arg(part, choices=c("common", "distinctive"))
                    if (part == "common") {
                        res <- x@VAF$common
                    } else if(part=="distinctive") {
                        if (is.null(block)) {
                            res <- x@VAF$dist
                        } else if (block == "1" | block == x@Names[1]) {
                            res <- x@VAF$dist$Block1
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@VAF$dist$Block2
                        } else if (block=="cross") {
                            res <- x@VAF$dist$cross
                        } else {
                            warning(paste(block, "is not an allowed option."))
                            res <- x@VAF$common
                        }
                    }
                }
            } else if (x@caMethod == "JIVE") {
                if (is.null(part)) {
                    res <- x@VAF
                } else {
                    part <- match.arg(part, choices=c("common", "distinctive"))
                    if (part == "common") {
                        res <- x@VAF$common
                    } else if (part == "distinctive") {
                        if (is.null(block)) {
                            res <- x@VAF$dist
                        } else if (block == "1" | block == x@Names[1]) {
                            res <- x@VAF$dist$Block1
                        } else if (block == "2" | block == x@Names[2]) {
                            res <- x@VAF$dist$Block2
                        } else {
                            warning(paste(block, "is not an allowed option."))
                            res <- x@VAF$dist
                        }
                    }
                }
            } else if (x@caMethod == "O2PLS") {
                warning("VAF per component cannot be calculated with O2PLS approach because components are not orthogonal")
                res <- NULL
            }
        }
        return(res)
    }
)