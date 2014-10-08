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


# InitialData ------------------------------------------------------------

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

# preprocessing -----------------------------------------------------------

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

# methodInfo --------------------------------------------------------------

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

# scores ------------------------------------------------------------------

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

# loadings ----------------------------------------------------------------

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

# VAF ---------------------------------------------------------------------

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