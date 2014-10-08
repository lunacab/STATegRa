
########  CLASS CONNECTING ELEMENTS

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

bioMap <- function(name, metadata, map)
{
    new("bioMap", name=name, metadata=metadata, map=map)
}