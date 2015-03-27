#' @include STATegRa_omicsCLUST_bioMap.R

##################################
###### PART2: CLASS FOR THE DISTANCES
##################################


setClass("bioDistclass",
         slots=c(name = "character",
                 distance = "matrix",
                 map.name= "character",
                 map.metadata = "list",
                 params= "list") #
)

####### INITIALIZE

setMethod(
    f="initialize",
    signature="bioDistclass",
    definition=function(.Object, name = "",
                        distance = matrix,
                        map.name = "",
                        map.metadata = list(),
                        params = list())
    {
        .Object@name <- name
        .Object@distance <-  distance
        .Object@map.name <- map.name
        .Object@map.metadata <- map.metadata
        .Object@params <- params
        return(.Object)
    }
)

###### ACCESSING THEM

setMethod(
    f= "getName",
    signature= "bioDistclass",
    definition=function(object){
        return(object@name)
    }
)

setGeneric(
    name= "getDist",
    def=function(object){standardGeneric("getDist")}
)

setMethod(
    f= "getDist",
    signature= "bioDistclass",
    definition=function(object){
        return(object@distance)
    }
)

setGeneric(
    name= "getMapName",
    def=function(object){standardGeneric("getMapName")}
)

setMethod(
    f= "getMapName",
    signature= "bioDistclass",
    definition=function(object){
        return(object@map.name)
    }
)


setGeneric(
    name= "getMapInfo",
    def=function(object){standardGeneric("getMapInfo")}
)

setMethod(
    f= "getMapInfo",
    signature= "bioDistclass",
    definition=function(object){
        return(object@map.metadata)
    }
)

setGeneric(
    name= "getDistParam",
    def=function(object){standardGeneric("getDistParam")}
)

setMethod(
    f= "getDistParam",
    signature= "bioDistclass",
    definition=function(object){
        return(object@params)
    }
)

setGeneric(
    "setDist<-",function(object,value)
    {
        standardGeneric("setDist<-")
    }
)

setReplaceMethod(
    f="setDist",
    signature="bioDistclass",
    definition=function(object,value)
    {
        object@distance <-value
        return (object)
    }
)

###### CONSTRUCTOR

#' @export
#' @title bioDistclass
#' @description Class to manage mappings between genomic features.
#' @usage bioDistclass(name, distance, map.name, map.metadata, params)
#' @param name Name assigned to the object
#' @param distance Matrix giving the distance between features
#' @param map.name Charactering giving the name of the bioMap object used to compute the distance
#' @param map.metadata List of parameters used to generate the mapping
#' @param params List of parameters used to generate the distance
bioDistclass <- function(name, distance, map.name, map.metadata, params)
{
    new("bioDistclass", name=name, distance=distance, map.name=map.name, map.metadata=map.metadata, params=params)
}

######  bioAggregate function

setGeneric(
    name= "bioAggregate",
    def=function(referenceFeatures,
                 reference = "genes",
                 mapping,
                 referenceData = NULL,
                 surrogateData,
                 filtering = NULL,
                 aggregation = "sum",
                 maxitems = NULL,
                 selectionRule = "maxFC",
                 expfac = NULL,...){
        standardGeneric("bioAggregate")}
)

setMethod(
    f= "bioAggregate",
    signature(referenceFeatures = "character",
              mapping = "bioMap",
              surrogateData = "ExpressionSet",
              referenceData = "ExpressionSet"),
    definition=function(referenceFeatures,
                        reference = "genes",
                        mapping,
                        referenceData = NULL,
                        surrogateData,
                        filtering = NULL,   #c(a,b,c)
                        aggregation = "sum",
                        maxitems = NULL,
                        selectionRule = "sd",
                        expfac = NULL,...){

        #######################################
        ## STEP 0. Validity of input objects
        #######################################

        #######################################
        ## STEP 1. Filtering mapping
        #######################################
        #### By now we consider that it needs to provide
        ## a first element with the name of the column.
        ## a second element with the selected value

        filtmapping = getMap(mapping)
        # eliminate duplicates
        filtmapping.2<-unique(filtmapping)
        if(nrow(filtmapping.2)!=nrow(filtmapping))
        {
            warning("Duplicates found in mapping. Only 1 representation for each link is conserved")
        }
        filtmapping<-filtmapping.2
        rm(filtmapping.2)

        #######################################
        ## STEP 2. Selecting surrogate data
        #######################################

        # Set the column that acts as a surrogate
        surrogate = setdiff(colnames(filtmapping), reference)
        # Filter the mapping for those features included in the referenceFeatures
        filtmapping = filtmapping[as.character(filtmapping[,reference]) %in% referenceFeatures,]
        # Filter the mapping for those surrogate variables included in the surrogate profiling
        filtmapping = filtmapping[as.character(filtmapping[,surrogate]) %in% rownames(exprs(surrogateData)),]
        # Generate the surrogate profiling to be used.
        surData = exprs(surrogateData)[as.character(unique(filtmapping[,surrogate])),]  ### replaced assayData

        # Define the number of surrogate features associated to each main feature.
        numsurXref = table(filtmapping[,reference])
        # At least 1 surrogate.
        numsurXref = names(which(numsurXref > 0)) # reference features with >maxitems surrogates
        # Feature eliminated: forcing at least a number equal to the max number of items to be used.

        if (length(numsurXref) > 0 & !is.null(maxitems)) {
            ####### The first approach is to map by using data.
            ###############################

            if (selectionRule == "maxcor") {
                ## referenceData must be provided
                if(is.null(referenceData))
                {
                    stop("The maxcor selection rule requires to provide referenceData")
                }
                nover<-colnames(surData) [colnames(surData) %in% colnames(exprs(referenceData))]
                if(length(nover)<3)
                {
                    stop("The maxcor selection rule requires at least 3 samples names shared by referenceData and surrogateData.")
                }else{
                    for(i in 1:length(numsurXref))
                    {
                        whichgo<-which(filtmapping[,reference]==numsurXref[i])
                        if(length(whichgo)>maxitems)
                        {
                            rm(surData.f)
                            names.f<- as.character(filtmapping[whichgo,surrogate])
                            surData.f<-surData[unique(filtmapping[whichgo,surrogate]),nover]
                            refData.f<-as.numeric(exprs(referenceData)[numsurXref[i],nover])
                            corgo<-apply(surData.f,1,function(x){
                                cor.test(as.numeric(x),
                                         refData.f, method="pearson")$estimate})
                            corgo.order<-corgo[order(abs(corgo),decreasing=TRUE)]
                            filtmapping = filtmapping[-which((filtmapping[,reference]==numsurXref[i]) &
                                                                 (as.character(filtmapping[,surrogate]) %in%
                                                                      names(corgo.order)[(maxitems+1):(length(corgo.order))])),]
                        }

                    }

                }
            }

            ######### The second idea is to map by using the standard deviation criteria
            ###############################

            if(selectionRule != "maxcor"){
                for (ff in numsurXref) {
                    # how many elements
                    whichgo<-which(as.character(filtmapping[,reference]) == ff)
                    if(length(whichgo)>maxitems)
                    {
                        # data for >maxitems surrogates features associated to reference feature ff:
                        subconj = surData[as.character(filtmapping[whichgo,surrogate]),]
                        # to select surrogates to reomove for this ff reference:
                        subconj = sur2remove(subconj, maxitems = maxitems,
                                             selectionRule = selectionRule, expfac = expfac)
                        filtmapping = filtmapping[-which((filtmapping[,reference]==ff) &
                                                             (as.character(filtmapping[,surrogate]) %in% subconj)),]
                    }
                }
            }

            ## surrogate Data
            surData = exprs(surrogateData)[as.character(filtmapping[,surrogate]),] #assayData
        }

        #######################################
        ## STEP 3. Aggregating surrogate features
        #######################################
        if (aggregation != "no") {
            surData = aggregate(surData,
                                by = list("reference" = filtmapping[,reference]),
                                as.name(aggregation))
            rownames(surData) = surData[,"reference"]
            surData = surData[,-1]
        }

        #######################################
        ## STEP 4. Return surrogate data
        #######################################

        return(list(data = surData, mapping = filtmapping))
    }
)


#' @export
#' @title bioDist
#' @aliases bioDist,character,character,bioMap,ExpressionSet,ExpressionSet-method
#' @description 
#' Function to compute a bioDistclass object from profile data and a mapping. For details of the process see the user's guide, but briefly the process involves using the mapping to identify reference features appropriate to each surrogate feature (if any), aggregating the surrogate data into pseudo-data for each reference feature, and then calculating the correlation distance between the reference features according to the surrogate data.
#' @usage bioDist(referenceFeatures=NULL, reference=NULL, mapping=NULL, 
#'                referenceData=NULL, surrogateData=NULL, filtering=NULL, 
#'                noMappingDist=NA, distance="spearman", aggregation="sum", 
#'                maxitems=NULL, selectionRule="maxFC", expfac=NULL, 
#'                name=NULL, ...)
#' @param referenceFeatures subset of features to be considered for the computation of the distances. If NULL then the features are first gathered from the features in referenceData. If referenceData is not provided then the list of features are gathered from mapping (bioMap class) and using the reference.
#' @param reference A character indicating the variable that is being used as features to compute distance between
#' @param mapping The mapping between feature types
#' @param referenceData ExpressionSet object with the data from the reference features.
#' @param surrogateData ExpressionSet object with the data from the surrogate features.
#' @param filtering A filtering for the bioMap class. To be implemented.
#' @param noMappingDist Distance value to be used when a reference feature do not map to any surrogate feature.  If "max", maximum indirect distance among the rest of reference features is taken. If NA, distance weights are re-scaled so this surrogate association is not considered. If a number then the missing values are replaces with that value.
#' @param distance Distance between features to be computed. Possible values are "pearson", "kendall", "spearman", "euclidean", "maximum", "manhattan", "canberra", "binary" and "minkowski". Default is "spearman".
#' @param aggregation Action to perform when a reference feature maps to more than one surrogate feature. Options are "max", "sum", "mean" or "median"  and the the values are aggregated according to the chosen statistic.
#' @param maxitems The maximum number of surrogate features per reference feature to be used, selected according to "selectionRule" parameter. Default is 2.
#' @param selectionRule Rule to select the surrogate features to be used (the number is determined by "maxitems"). It can be one of the following: (1)  "maxcor" those presenting maximum correlation with corresponding main feature; in this case "referenceData" must be provided and the columns must overlap in at least 3 samples; (2) "maxmean": average across samples is computed and those features with higher mean are selected;  case (3) is simmilar to (2) but considering other statistics: "maxmedian", "maxdiff", "maxFC", "sd" , "ee".
#' @param expfac Not in use yet.
#' @param name Character that describes the nature of the bioDist class computed
#' @param ... extra arguments passed to \code{\link{dist}}, eg "p=value" for the power used if calculating minkowski distance
#' @return An object of class \code{bioDistclass} containing distances between the features in \code{surrogateData}.
#' @author David Gomez-Cabrero
#' @template omicsCLUST_examples_common
setGeneric(
    name= "bioDist",
    def=function(referenceFeatures = NULL,
                 reference = NULL,
                 mapping = NULL,
                 referenceData = NULL,
                 surrogateData = NULL,
                 filtering = NULL,
                 noMappingDist = NA,
                 distance = "spearman",
                 aggregation = "sum",
                 maxitems = NULL,
                 selectionRule = "maxFC",
                 expfac = NULL,
                 name=NULL,...){
        standardGeneric("bioDist")}
)


setMethod(
    f= "bioDist",
    signature(referenceFeatures = "character",
              reference = "character",
              mapping = "bioMap",
              surrogateData = "ExpressionSet",
              referenceData = "ExpressionSet"),
    definition=function(referenceFeatures=NULL,
                        reference=NULL,
                        mapping=NULL,
                        referenceData = NULL,
                        surrogateData = NULL,
                        filtering = NULL,
                        noMappingDist = NA,
                        distance = "spearman",
                        aggregation = "sum",
                        maxitems = NULL,
                        selectionRule = "maxFC",
                        expfac = NULL,
                        name = NULL,...){

        ## validity of input objects

        if(is.null(referenceFeatures))
        {
            warning("referenceFeatures set is not provided.")
            if(!is.null(referenceData))
            {
                referenceFeatures<-rownames(exprs(referenceData))
                warning("referenceFeatures is generated from referenceData features")
            }else{
                if(!is.null(mapping) & !is.null(reference))
                {
                    referenceFeatures<-unique(getMap(mapping)[,reference])
                    warning("referenceFeatures is generated from the bioMap object")
                }else{
                    stop("It is not possible to generate referenceFeatures.
               When referenceFeatures is not provided it is required either (1) referenceData or (2) a bioMap object (mapping) and the reference provided.")
                }
            }
        }


        ## reference = surrogate
        ### In this case the distance is compute throug the profiles of the main features.

        if (is.null(mapping)) {
            warning("When the bioMap object is not provided, the mapping makes use of the referenceData provided.")
            if(is.null(referenceData))
            {
                stop("There is no bioMap object (mapping) neither referenceData provided. \nIn this case it is not possible to compute the distance.")
            } else{
                surData = exprs(referenceData)[referenceFeatures,] #assayData
                if (!is.null(distance)) {
                    mydist = dist(t(surData), method = distance,...)
                }
            }
            ## reference != surrogate -> mapping has to be provided
        } else {

            myaggre = bioAggregate(referenceFeatures = referenceFeatures,
                                   reference = reference,
                                   mapping = mapping,
                                   referenceData = referenceData,
                                   surrogateData = surrogateData,
                                   filtering = filtering,
                                   aggregation = aggregation,
                                   maxitems = maxitems,
                                   selectionRule = selectionRule,
                                   expfac = expfac)
            surData = myaggre$data
            if (is.null(distance)) {  # no distance is computed (e.g. for k-means algorithm)
                mydist = NULL
            } else {  # distance is to be computed

                if (aggregation == "no") {
                    # all pairwise distances must be computed and aggregated afterwards:
                    # To be considered for development.
                    # By now we consider all the surrogates are computed first.

                } else {
                    # Compute distance
                    if(distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski"))
                    {
                        surDist = as.matrix(dist(surData, method = distance,...))
                    } else { #then it is correlation "pearson", "kendall", "spearman"
                        surDist = as.matrix(cor(t(surData), method = distance))
                    }

                    # Creating EMPTY distance matrix for reference features
                    mydist = matrix(noMappingDist,
                                    nrow = length(referenceFeatures),
                                    ncol = length(referenceFeatures))
                    rownames(mydist) = colnames(mydist) = referenceFeatures

                    # Complete empty distance matrix with references without surrogates
                    mydist[colnames(surDist),colnames(surDist)]<-surDist
                    if(noMappingDist=="max")
                    {
                        mydist[mydist=="max"]<-max(surDist)
                    }
                }
            }


        }

        return(new("bioDistclass",
                   name = name,
                   distance = mydist,
                   map.name = getName(mapping),
                   map.metadata = getMetadata(mapping),
                   params = list(filtering=filtering,distance=distance,aggregation=aggregation,
                                 maxitems=maxitems,selectionRule=selectionRule,expfac=expfac)))
    }
)


#############################################################################
# Auxiliary functions -----------------------------------------------------
############################################################################


## Function to select maxitems surrogates from a given reference feature
## according to selectionRule

sur2remove = function (data, maxitems = 3, selectionRule = "maxFC", expfac = NULL) {

    # Considering group information
    if (!is.null(expfac)) {
        # Averaging replicates within each experimental condition
        data = aggregate.con(data, groupby = expfac, aggregation = "rowMeans")
    }

    # Maximum fold-change
    if (selectionRule == "maxFC") {
        data = M.max.min(data)
        data = sort(data, decreasing = TRUE)  ## No specific actions when draws.
        data = names(data)[(maxitems+1):length(data)]
    }

    # Maximum difference
    if (selectionRule == "maxdiff") {
        data = rango(data)
        data = sort(data, decreasing = TRUE)  ## No specific actions when draws.
        data = names(data)[(maxitems+1):length(data)]
    }

    # Maximum standard deviation
    if (selectionRule == "sd") {
        data = apply(data,1,sd)
        data = sort(data, decreasing = TRUE)  ## No specific actions when draws.
        data = names(data)[(maxitems+1):length(data)]
    }

    # Maximum standard error
    if (selectionRule == "ee") {
        data1 = apply(data,1,sd)
        data2 = apply(data,1,mean)
        data = data1/data2
        data = sort(data, decreasing = TRUE)  ## No specific actions when draws.
        data = names(data)[(maxitems+1):length(data)]
    }

    return(data)  # names of selected surrogate features (maxitems)
}



## Function to aggregate replicates within a experimental conditions

# groupby is a vector containing the experimental group for each column of data
# aggregation: function to apply on a matrix and returning a vector

aggregate.con = function (data, groupby, aggregation = "rowMeans"){

    grupos = unique(groupby)

    agregados = lapply(grupos, function (gg) as.matrix(data[,grep(gg, groupby)]))

    agregados = sapply(agregados, as.name(aggregation))

    colnames(agregados) = grupos

    return (agregados)
}


## Function to compute maximum fold-change across experimental conditions
M.max.min = function (data){

    data = sinceros(data, k = NULL)

    M = log2(apply(as.matrix(data), 1, max, na.rm = TRUE) / apply(as.matrix(data), 1, min, na.rm = TRUE))
    names(M) = rownames(data)
    return(M)
}


## Functions to change 0s to a constant k
noceros = function (x, num = TRUE, k = 0) {

    nn = sum(x > k)

    if (num) {
        nn

    } else {
        if(nn > 0) { which(x > k) } else { NULL }
    }
}



sinceros = function (datos, k) {
    datos0 = datos

    if (is.null(k)) {

        mini0 = min(datos[noceros(datos, num = FALSE, k = 0)])

        kc = mini0/2

        datos0[datos0 == 0] = kc

    } else {

        datos0[datos0 == 0] = k

    }

    datos0
}


## Function to compute maximum difference across experimental conditions
rango = function (data){
    mirango = apply(as.matrix(data), 1, function (y) diff(range(y, na.rm = TRUE)))
    names(mirango) = rownames(data)
    return(mirango)
}



