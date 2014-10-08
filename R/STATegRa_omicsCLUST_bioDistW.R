
################################################################################
##### bioDistw CLASS
######## this one stores the weighted distance through surrogates features.
################################################################################

setClass(
    "bioDistWclass",
    slots=c(name = "character",
            distance = "matrix",
            map.name= "character",
            map.metadata = "list",
            params= "list",
            weights = "vector"))

####### INITIALIZE

setMethod(
    f="initialize",
    signature="bioDistWclass",
    definition=function(.Object, name = "",
                        distance = matrix,
                        map.name = "",
                        map.metadata = list(),
                        params = list(),
                        weights = vector())
    {
        .Object@name <- name
        .Object@distance <-  distance
        .Object@map.name <- map.name
        .Object@map.metadata <- map.metadata
        .Object@params <- params
        .Object@weights <- weights
        return(.Object)
    }
)

####### INITIALIZE

setGeneric(
    name= "getWeight",
    def=function(object){standardGeneric("getWeight")}
)

setMethod(
    f= "getWeight",
    signature= "bioDistWclass",
    definition=function(object){
        return(object@weights)
    }
)

setMethod(
    f= "getName",
    signature= "bioDistWclass",
    definition=function(object){
        return(object@name)
    }
)

setMethod(
    f= "getDist",
    signature= "bioDistWclass",
    definition=function(object){
        return(object@distance)
    }
)

# bioDistW CLASS --------------------------------------------------------

setGeneric(
    name= "bioDistW",
    def=function(referenceFeatures,
                 bioDistList,
                 weights){
        standardGeneric("bioDistW")}  ### update when finished!!!!
)


setMethod(
    f= "bioDistW",
    signature(referenceFeatures = "character",
              bioDistList = "list",
              weights="matrix"),
    definition=function(referenceFeatures,
                        bioDistList,
                        weights=1.){

        #### STEP 1:
        ## Validity of input objects

        # Check weights
        if(is.null(dim(weights)))
        {
            if(weights==1)
            {
                weigths<-matrix(1/length(bioDistList),1,length(bioDistList))
            }else
            {
                stop("weight needs to be defined")
            }
        }else{
            if(ncol(weights)!=length(bioDistList))
            {
                stop("The number of columns different from the number of list")
            }
            sumr<-apply(weights,1,sum)
            if(sum(sumr!=1)>0)
            {
                warning("For some row the total sum of weights is different than 1")
            }
        }

        # Number of features
        if(length(referenceFeatures)!=unique(length(referenceFeatures)))
        {
            stop("There are features names repeated")
        }

        # Check all features in all distances
        for(i in 1:length(bioDistList))
        {
            if(length(referenceFeatures)!=sum(unique(rownames(getDist(bioDistList[[i]]))) %in% referenceFeatures))
            {
                stop("There are features not represented in all distance matrices")
            }
        }
        # Set colnames
        colnames(weights)<-paste("V",1:ncol(weights),sep="-")
        for(i in 1:ncol(weights))
        {
            if(!is.na(getName(bioDistList[[i]])))
            {
                colnames(weights)[i]<-getName(bioDistList[[i]])
            }
        }

        #### STEP 2: COMPUTE ALL WEIGTHED DISTANCES
        ##
        listDistW<-vector("list",nrow(weights))
        for(i in 1:length(listDistW))
        {
            m1<-getDist(bioDistList[[1]])[referenceFeatures,referenceFeatures]*weights[i,1]
            for(j in 2:length(bioDistList))
            {
                m1 <- m1+ getDist(bioDistList[[j]])[referenceFeatures,referenceFeatures]*weights[i,j]
            }
            listDistW[[i]]<-new("bioDistWclass",
                                name = paste(colnames(weights),collapse="-"),
                                distance = m1,
                                weights = weights[i,])
        }

        listDistW
    })

########## COMPUTE THE DISTANCES FOR A GIVEN RELEVANT GENE
# listDistW<-bioDistWList

setGeneric(
    name= "bioDistFeature",
    def=function(Feature,
                 listDistW,
                 threshold.cor){
        standardGeneric("bioDistFeature")}
)


setMethod(
    f= "bioDistFeature",
    signature(Feature = "character",
              listDistW = "list",
              threshold.cor="numeric"),
    definition=function(Feature = "character",
                        listDistW = "list",
                        threshold.cor="numeric"){

        ##### STEP 1: SELECT OTHER FEATURES WITH HIGHER CORRELATION
        list.fea<-c()
        for(i in 1:length(listDistW))
        {#i<-1
            list.fea<-c(list.fea,colnames(getDist(listDistW[[i]]))[abs(getDist(listDistW[[i]])[Feature,])>threshold.cor])
        }
        list.fea<-unique(list.fea)
        list.fea<-list.fea[list.fea!=Feature]

        #### STEP 2:
        if(length(list.fea)>0)
        {
            RESULT<-matrix(NA,length(list.fea),length(listDistW))
            rownames(RESULT)<-list.fea
            colnames(RESULT)<-paste("Names_",1:length(listDistW),sep="")
            for(i in 1:length(listDistW) )
            {
                colnames(RESULT)[i]<-paste(getWeight(listDistW[[i]]),collapse="-")
                RESULT[list.fea,i]<-getDist(listDistW[[i]])[Feature,list.fea]
            }
        }else{
            warning("There is no Feature associated to selected feature given the threshold selection for any weighted distance")
        }

        RESULT
    }
)


########## COMPUTE THE DISTANCES BETWEEN THEM AND PLOT THEM.
# listDistW<-bioDistWList

setGeneric(
    name= "bioDistWPlot",
    def=function(referenceFeatures,
                 listDistW,
                 method.cor){
        standardGeneric("bioDistWPlot")}  ### update when finished!!!!
)


setMethod(
    f= "bioDistWPlot",
    signature(referenceFeatures = "character",
              listDistW = "list",
              method.cor="character"),
    definition=function(referenceFeatures,
                        listDistW,
                        method.cor="spearman"){

        # Compute Distances.
        dist.M<-matrix(1,length(listDistW),length(listDistW))
        for(i in 2:(nrow(dist.M)-1))
        {
            for(j in (i+1):nrow(dist.M))
            {
                dist.M[i,j]<-cor.test(getDist(listDistW[[i]])[referenceFeatures,referenceFeatures],
                                      getDist(listDistW[[j]])[referenceFeatures,referenceFeatures],
                                      method=method.cor)$estimate
                dist.M[j,i]<-dist.M[i,j]
            }
        }

        # Compute weights
        WEIGHT.ALL<-matrix(NA,length(listDistW),length(getWeight(listDistW[[i]])))
        colnames(WEIGHT.ALL)<-names(getWeight(listDistW[[i]]))
        rownames(WEIGHT.ALL)<-paste("Name",1:length(listDistW),sep="")

        for(i in 1:length(listDistW))
        {
            WEIGHT.ALL[i,]<-as.numeric(getWeight(listDistW[[i]]))
            rownames(WEIGHT.ALL)[i]<-paste(as.numeric(getWeight(listDistW[[i]])),collapse="-")
        }

        # Plot Distances.
        fit <- cmdscale(dist.M,eig=TRUE, k=2) # k is the number of dim
        fit # view results
        x <- fit$points[,1]
        y <- fit$points[,2]

        oldmfrow <- par("mfrow")
        par(mfrow = c(2, 2*(round(ncol(WEIGHT.ALL))-1)))
        plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
             main="MDS", pch=19)
        for(i in 1:ncol(WEIGHT.ALL))
        {
            plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
                 main=colnames(WEIGHT.ALL)[i],
                 pch=19,col=rgb(0,0,0,WEIGHT.ALL[,i]))
        }
        plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
             main="MDS", pch=19)

        textxy(x, y, rownames(WEIGHT.ALL))

        par(mfrow=oldmfrow)
    }

)


bioDistFeaturePlot<-function(data)
{
    colors = c(seq(-1,0,length=100),seq(0,1,length=100))
    my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 199)
    #par(mar=c(12,2,2,12))#,mai=c(12,2,2,12))
    heatmap.2(t(data),
              Rowv = FALSE,
              Colv=TRUE,
              distfun = function(x){dist(x,method="manhattan")},
              hclustfun = hclust,
              dendrogram = c("col"),
              symm = FALSE,
              cexRow=1,cexCol=1,margins=c(12,8),
              breaks=colors,
              col=my_palette,
              scale = c("none"),
              na.rm=TRUE,
              trace="none")
}



