####### DATA MANAGEMENT #######

# CreateOmicsExpressionSet ------------------------------------------------
## Function to create a omics Expression datasets given the information separately in three matrices.
## Matrices about phenotype and features are optional
## INPUT:
##    Data: (matrix) Omics data
##    pData=NULL: (matrix) Data associated to the samples/phenotype (columns). rownames(pData)==colnames(exprs)
##    pDataDescr=NULL:  (matrix) Description of the phenotypic data
##    feaData=NULL: (matrix) Data associated to the variables/features annotation (rows). rownames(fData)==rownames(exprs)
##    feaDataDescr=NULL: (matrix) Description of the feature annotation
## OUTPUT:
##    ExpressionSet with the data provided
setGeneric(
    name="createOmicsExpressionSet",
    def=function(Data,pData=NULL,pDataDescr=NULL,feaData=NULL,feaDataDescr=NULL){
        standardGeneric("createOmicsExpressionSet")}
)

setMethod(
    f="createOmicsExpressionSet",
    signature=signature(
        Data="matrix"),
    definition=function(Data,pData,pDataDescr,feaData,feaDataDescr){
        ## Data
        exprs <- as.matrix(Data)
        ## Phenotypic data
        
        if(is.null(pData)){
            if (is.null(colnames(exprs))) {
                stop("If the data has no colnames, pData with names must be provided")
            }
            pData <- data.frame(Names=colnames(exprs))
            rownames(pData) <- colnames(exprs)
            pDataDescr <-c("Sample name")
        } else if (is.null(pDataDescr)) {
            pDataDescr <- colnames(pData)
        }
        pData <- as.data.frame(pData)
        if (is.null(rownames(pData)) | !all(rownames(pData)==colnames(exprs))){
            stop("pData rownames and exprs colnames must be the same")
        }
        metadata <- data.frame(labelDescription=pDataDescr)
        rownames(metadata) <- colnames(pData)
        phenoData <- new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
        ## Feature data
        if(is.null(feaData)){
            if (is.null(rownames(exprs))) {
                stop("If the data has no rownames, feaData with names must be provided")
            }
            feaData <- data.frame(Name=rownames(exprs))
            rownames(feaData) <- rownames(exprs)
            feaDataDescr <- c("Feature name")
        } else if (is.null(feaDataDescr)) {
            feaDataDescr <- colnames(feaData)
        }
        fData <- as.data.frame(feaData)
        if (is.null(rownames(fData)) | !all(rownames(fData)==rownames(exprs))){
            stop("annot rownames and exprs rownames must be the same")
        }
        metadata2 <- data.frame(labelDescription=feaDataDescr)
        rownames(metadata2) <- colnames(fData)
        featureData <- new("AnnotatedDataFrame",data=fData,varMetadata=metadata2)
        ## ExpressionSet
        es <- new("ExpressionSet",exprs=exprs,phenoData=phenoData,featureData=featureData)
        return(es)
    }
)


STATegRaUsersGuide <- function(view=TRUE)
    #	Find and optionally view edgeR STATegRa Guide
    #	David Gomez-Cabrero, following Gordon Smyth example on edgeR
    #	23 May 2012.
{
    f <- system.file("doc", "STATegRa.html", package="STATegRa")
    if(view) {
        browseURL(f)
    }
    return(f)
}
