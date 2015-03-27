
########  CLASS CONNECTING ELEMENTS
# we need `methods` for all the setClass/setMethod/setGeneric stuff, importing it arbitarily here
#' @import methods
setClass("bioMap",
         slots = c(name = "character",
                   metadata = "list",
                   map= "data.frame")
)
#### Vocabulary:
# name: is the name of the mapping, as a reference
# metadata: is a list of elements including references to the origin of the data
#            meta-data list must include the following variables:
#             type_v1: the type of variable 1
#             type_v2: the type of variable 2
#             source_database: the source data-base it was extracted from.
#             data_extraction: the date the mapping was gathered.
#
# map: is the actual mapping, data frame:
#           column 1: name "Var1" features of the variable 1
#           column 2: name "Var2" features of the variable 2
#           column 3 to x: the following columns are possible to exist, but not necessary:
#               "distance": distance between features
#               "location": exon, intron, etc...
#               "evidence": the origin of the evidence: "experimental", "computational",...
# Note that the vocabulary for data types is: transcript, gene, protein, metabolite, SNP , genomicRegion, etc.

####### INITIALIZE

setMethod(
    f="initialize",
    signature="bioMap",
    definition=function(.Object, name = "",
                        metadata = list(),
                        map = data.frame())
    {
        .Object@name <- name
        .Object@metadata <- metadata
        .Object@map <- map
        return(.Object)
    }
)

#######  ACCESSORS

setGeneric(
    name= "getName",
    def=function(object){standardGeneric("getName")}
)

setMethod(
    f= "getName",
    signature= "bioMap",
    definition=function(object){
        return(object@name)
    }
)

setGeneric(
    name= "getMetadata",
    def=function(object){standardGeneric("getMetadata")}
)

setMethod(
    f= "getMetadata",
    signature= "bioMap",
    definition=function(object){
        return(object@metadata)
    }
)

setGeneric(
    name= "getMap",
    def=function(object){standardGeneric("getMap")}
)

setMethod(
    f= "getMap",
    signature= "bioMap",
    definition=function(object){
        return(object@map)
    }
)

##########################################
##### SET
##########################################

setGeneric(
    "setName<-",function(object,value)
    {
        standardGeneric("setName<-")
    }
)

setReplaceMethod(
    f="setName",
    signature="bioMap",
    definition=function(object,value)
    {
        object@name <- value
        return (object)
    }
)

setGeneric(
    "setMetadata<-",function(object,value)
    {
        standardGeneric("setMetadata<-")
    }
)

setReplaceMethod(
    f="setMetadata",
    signature="bioMap",
    definition=function(object,value)
    {
        object@metadata <-value
        return (object)
    }
)

setGeneric(
    "setMap<-",function(object,value)
    {
        standardGeneric("setMap<-")
    }
)

setReplaceMethod(
    f="setMap",
    signature="bioMap",
    definition=function(object,value)
    {
        object@map <-value
        return (object)
    }
)

###### CONSTRUCTOR

#' @export
#' @title bioMap
#' @description Function to generate a bioMap object.
#' @usage bioMap(name, metadata, map)
#' @param name Name to assign the object
#' @param metadata A list with information of the mapping. Elements expected in the list are: (1) "type_v1" and "type_v2", refer to the nature of the features mapped; a vocabulary we recommend is "gene", "mRNA", "miRNA", "proteins", etc. (2) "source_database", provides information on the source of the mapping; from a specific data-base e.g. "targetscan.Hs.eg.db" to a genomic location mapping. (3) "data_extraction" stores information on the data the mapping was generated or downloaded.
#' @param map A data.frame object storing the mapping. The data.frame may inclue an unlimited number of columns, however the first column must be named "Var1" and refer to the elements of "type_v1" and simmilarly for the second column ("Var2", "type_v2").
#' @return An object of class \code{bioMap}
#' @author David Gomez-Cabrero
#' @examples 
#' data(STATegRa_S2)
#' map.gene.miRNA<-bioMap(name = "Symbol-miRNA",
#'                        metadata = list(type_v1="Gene",type_v2="miRNA",
#'                                        source_database="targetscan.Hs.eg.db",
#'                                        data_extraction="July2014"),
#'                        map=mapdata)
bioMap <- function(name, metadata, map)
{
    new("bioMap", name=name, metadata=metadata, map=map)
}