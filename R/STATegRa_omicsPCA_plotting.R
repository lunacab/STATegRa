####### PLOTS #######

# g_legend ----------------------------------------------------------------
## Extract the legend for a general legend with grid.arrange
## code from https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
## INPUT:
##    a.gplot: (ggplot) A plot
## OUTPUT:
##    Legend of the plot
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}



# plotVAF (Associated to caClass) ------------------------------------------------------------
## INPUT:
##    object: (scaClass) Object with results of Simultaneous Components Analysis (SCA)
##    mainTitle: (character) Title of the plot
## OUTPUT:
##    (plot): VAF plot, one for each block and one for the joined data
#' @export
#' @import ggplot2
#' @import grid
#' @import Biobase
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @title Plot VAF (Variance Explained For) from Component Analysis
#' @aliases plotVAF,caClass-method
#' @description 
#' This function visualises the VAF results from component analysis. The input is a \code{\link{caClass-class}} object from \code{\link{omicsCompAnalysis}}. VAF cannot be calculated if mode "O2PLS" was used. The plots for modes "DISCOSCA" and "JIVE" are different since DISCO-SCA distinctive components have some VAF in the other block. This VAF can be interpreted as an error in the rotation.
#' @usage plotVAF(object, mainTitle="")
#' @param object \code{\link{caClass-class}} object containing component analysis results
#' @param mainTitle Plot title
#' @return \code{ggplot} object
#' @author Patricia Sebastian-Leon
#' @examples 
#' data("STATegRa_S3")
#' require(ggplot2)
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA,pDataDescr=c("classname"))
#' # Omics components analysis
#' discoRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                               method="DISCOSCA",Rcommon=2,Rspecific=c(2,2),
#'                               center=TRUE,scale=TRUE,weight=TRUE)
#' jiveRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                              method="JIVE",Rcommon=2,Rspecific=c(2,2),
#'                              center=TRUE,scale=TRUE,weight=TRUE)
#' 
#' # DISCO-SCA plotVAF
#' plotVAF(discoRes)
#' 
#' # JIVE plotVAF
#' plotVAF(jiveRes)
setGeneric(
    name="plotVAF",
    def=function(object,mainTitle=""){standardGeneric("plotVAF")}
)

setMethod(
    f="plotVAF",
    signature=signature(
        object="caClass"),
    definition=function(object,mainTitle){
        if(object@caMethod=="DISCO-SCA"){
            df_common <- data.frame(VAF=c(object@VAF$common),
                                    block=rep(object@Names,object@commonComps),
                                    comp=rep(1:object@commonComps,each=2),
                                    group=rep("common"),2*object@commonComps)
            df_dist <- data.frame(VAF=c(object@VAF$dist$Block1[1:object@distComps[1]],object@VAF$dist$Block2[1:object@distComps[2]],
                                        object@VAF$dist$cross[1,1:object@distComps[2]],object@VAF$dist$cross[2,1:object@distComps[1]]),
                                  comp=c(1:object@distComps[1],1:object@distComps[2],1:object@distComps[2],1:object@distComps[1]),
                                  block=c(rep(object@Names[1],object@distComps[1]),rep(object@Names[2],object@distComps[2]),
                                          rep(object@Names[1],object@distComps[2]),rep(object@Names[2],object@distComps[1])),
                                  group=c(rep(object@Names[1],object@distComps[1]),rep(object@Names[2],object@distComps[2]),
                                          rep(object@Names[2],object@distComps[2]),rep(object@Names[1],object@distComps[1])))
            upperlim <- max(df_common$VAF,df_dist$VAF)
            p_common <- ggplot(df_common,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("Proportion of VAF")+
                ggtitle("Common components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=10),legend.position="bottom")
            p_dist <- ggplot(df_dist,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("")+
                ggtitle("Distinctive components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=10),legend.position="bottom")
            mylegend <- g_legend(p_common)
            grid.arrange(arrangeGrob(p_common+ theme(legend.position="none"),p_dist+ theme(legend.position="none"),ncol=2,widths=c(1/3,2/3)),
                         mylegend,nrow=2,heights=c(10, 1),main=textGrob(label=mainTitle,gp=gpar(fontsize=20,cex=1)))
        } else if (object@caMethod=="JIVE"){
            df_common <- data.frame(VAF=c(object@VAF$common),
                                    block=rep("common",object@commonComps),
                                    comp=1:object@commonComps,
                                    group=rep("common"),object@commonComps)
            df_dist <- data.frame(VAF=c(c(object@VAF$dist$Block1),c(object@VAF$dist$Block2)),
                                  block=c(rep(object@Names[1],object@distComps[1]),rep(object@Names[2],object@distComps[2])),
                                  comp=c(1:object@distComps[1],1:object@distComps[2]),
                                  group=c(rep(object@Names[1],object@distComps[1]),rep(object@Names[2],object@distComps[2])))
            upperlim <- max(df_common$VAF,df_dist$VAF)
            p_common <- ggplot(df_common,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge",fill="#666666")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("Proportion of VAF")+
                ggtitle("Common components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=10),legend.position="bottom")
            p_dist <- ggplot(df_dist,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("")+
                ggtitle("Distinctive components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=10),legend.position="bottom")
            grid.arrange(arrangeGrob(p_common+ theme(legend.position="none"),p_dist+ theme(legend.position="none"),ncol=2,widths=c(1/3,2/3)),
                         nrow=2,heights=c(10, 1),main=textGrob(label=mainTitle,gp=gpar(fontsize=20,cex=1)))
        } else if(object@caMethod=="O2PLS"){
            warning("VAF per component cannot be calculated with O2PLS approach because components are not orthogonal")
        }
    }
)



# plotting ----------------------------------------------------------------
setGeneric(
  name="plotting",
  def=function(df,comps,xname,yname,color,what,type,combined,block,shape=NULL,labels=NULL,title,xlabel=NULL,
               ylabel=NULL,background=FALSE,palette=NULL,pointSize=4,labelSize=12,axisSize=8,titleSize=20,sizeValues=c(2,4),shapeValues=c(17,0)){standardGeneric("plotting")}
)


setMethod(
  f="plotting",
  signature=signature(
    df="data.frame",
    comps="numeric",
    xname="character",
    yname="character",
    color="character"),
  definition=function(df,comps,xname,yname,color,what,type,combined,block,shape,labels,title,xlabel,ylabel,background,palette,pointSize,
                      labelSize,axisSize,titleSize,sizeValues,shapeValues){
    
    if (what=="loadings"){
      df <- na.omit(df[,c(xname,yname)])
    }
    
    p <- ggplot(df,aes_string(xname,yname))
    
    if (what!="both"){
      if (!is.null(color) & !is.null(shape)){
        p <- p+geom_point(aes(color=color,shape=shape),size=pointSize)
      }else if (!is.null(color)){
          p <- p+geom_point(aes(color=color),size=pointSize)
      }else if (!is.null(shape)){
        p <- p+geom_point(aes(shape=shape),size=pointSize)
      }else{
        p <- p+geom_point(size=pointSize)
      }
    }else{
      p <- p + geom_point(data=df[df$class=="loadings",], na.rm=TRUE, aes(shape=class, size=class))+
        geom_point(data=df[df$class=="scores",],aes(color=color, shape=class, size=class))+
        scale_size_manual(name="",values=sizeValues)+
        scale_shape_manual(name="",values=shapeValues) 
    }

    if(!is.null(labels)){
      p <- p+geom_text(alpha=.5,size=pointSize,aes_string(label=labels),vjust=-.7)
    }  
    xlim <- ylim <- c(max(df[,xname],df[,yname],na.rm=TRUE),min(df[,xname],df[,yname], na.rm=TRUE))
    p <- p+xlim(xlim)+ylim(ylim)+coord_fixed()
    if (!background){
      p <- p+theme_bw()
    }
    if (is.null(title)){
      if(combined & type=="individual"){
        title <- paste(what,"plot","- combined",type,"components", sep=" ")
      }else if(!combined & type=="individual"){
        title <- paste(what,"plot","-",type,"components \nblock ",block, sep=" ")
      }else{
        title <- paste(what,"plot","-",type,"components", sep=" ")
      }
    }
    if (is.null(xlabel)){
      if(type=="individual"){
        xlabel <- paste("Component",comps[1],"- dist \n",block[1], sep=" ") 
      }else{
        xlabel <-  paste("Component",comps[1],"- common\n",block[1], sep=" ") 
      }
    }
    if (is.null(ylabel)){
      if(type=="common"){
        ylabel <- paste("Component",comps[2],"- common\n",block[2], sep=" ")
      }else{
        ylabel <- paste("Component",comps[2],"- dist \n",block[2], sep=" ")
      }
    }
    if (!is.null(palette)){
      p <- p+scale_colour_manual(values=palette)
    }
    p <- p+ggtitle(title)+xlab(xlabel)+ylab(ylabel)
    p <- p+theme(legend.position="bottom",
                 axis.text = element_text(size = axisSize),
                 axis.title = element_text(size = labelSize),
                 plot.title = element_text(size=titleSize, face="bold"))
    p
  }
)



# plotRes (Associated to caClass) --------------------------------------------
## Parameters:
##   what=c("scores","loadings","both"): plot scores, loadings or both
##   type=c("common","individual","both"): plot common components, individual components or a combination (common and individual)
##   combined=c(T,F): Logical indicating whether to combine components from different blocks
##   comps=c(comp_x,comp_y): components to plot. Indicate x and y, if combined=T, x=component_block1 and y=component_block2
#' @export
#' @import ggplot2
#' @import MASS
#' @importFrom  gridExtra arrangeGrob grid.arrange
#' @import Biobase
#' @title Plot component analysis results
#' @aliases plotRes,caClass,numeric,character,character,logical-method
#' @description 
#' Plot scatterplots of scores or loadings, for common and distinctive parts as well as combined plots.
#' @usage plotRes(object, comps=c(1, 2), what, type, combined, block=NULL, 
#'                color=NULL, shape=NULL, labels=NULL, title=NULL, xlabel=NULL, ylabel=NULL, background=TRUE, 
#'                palette=NULL, pointSize=4, labelSize=NULL, 
#'                axisSize=NULL, titleSize=NULL, sizeValues = c(2,4), shapeValues = c(17, 0))
#' @param object \code{\link{caClass-class}} containing component analysis results
#' @param comps If combined=FALSE, it indicates the x and y components of the type and block chosen. If \code{combined=TRUE}, it indicates the component to plot for the first block of information and the component for the second block of information to plot together. By default the components are set to c(1,2) if \code{combined=FALSE} and to c(1,1) if \code{combined=TRUE}.
#' @param what Either "scores", "loadings" or "both"
#' @param type Either "common", "individual" or "both"
#' @param combined Logical indicating whether to make a simple plot of two components from one block, or components from different blocks
#' @param block Which block to plot, either "1" or "2" or the name of the block.
#' @param color Character specifying a pData column from the original data to use to color points
#' @param shape Character specifying a pData column to select point shape
#' @param labels Character specifying a pData column from which to take point labels
#' @param title Main title
#' @param xlabel x-axis name
#' @param ylabel y-axis name
#' @param background Logical specifying whether to make a grey background
#' @param palette Vector giving the color palette for the plot
#' @param pointSize Size of plot points
#' @param labelSize Size of point labels if not NULL
#' @param axisSize Size of axis text
#' @param titleSize Size of title text
#' @param sizeValues Vector containing sizes for scores and loadings
#' @param shapeValues Vector indicating the shapes for scores and loadings
#' @return \code{ggplot} object
#' @author Patricia Sebastian-Leon
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,
#'                                pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA,pDataDescr=c("classname"))
#' # Omics components analysis
#' discoRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                               method="DISCOSCA",Rcommon=2,Rspecific=c(2,2),
#'                               center=TRUE,scale=TRUE,weight=TRUE)
#' jiveRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                              method="JIVE",Rcommon=2,Rspecific=c(2,2),
#'                              center=TRUE,scale=TRUE,weight=TRUE)
#' 
#' o2plsRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                               method="O2PLS",Rcommon=2,Rspecific=c(2,2),
#'                               center=TRUE,scale=TRUE,weight=TRUE)
#' 
#' # Scatterplot of scores variables associated to common components
#' 
#' # DISCO-SCA
#' plotRes(object=discoRes,comps=c(1,2),what="scores",type="common",
#'         combined=FALSE,block=NULL,color="classname",shape=NULL,labels=NULL,
#'         background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
#'         axisSize=NULL,titleSize=NULL)
#' # JIVE
#' plotRes(object=jiveRes,comps=c(1,2),what="scores",type="common",
#'        combined=FALSE,block=NULL,color="classname",shape=NULL,labels=NULL,
#'         background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
#'         axisSize=NULL,titleSize=NULL)
#' 
#' # O2PLS
#' # Scatterplot of scores variables associated to common components
#' # Associated to first block
#' p1 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
#'               combined=FALSE,block="expr",color="classname",shape=NULL,
#'               labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'               labelSize=NULL,axisSize=NULL,titleSize=NULL)
#' # Associated to second block
#' p2 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
#'               combined=FALSE,block="mirna",color="classname",shape=NULL,
#'               labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'               labelSize=NULL,axisSize=NULL,titleSize=NULL)
#' 
#' # Combined plot of scores variables assocaited to common components
#' plotRes(object=o2plsRes,comps=c(1,1),what="scores",type="common",
#'         combined=TRUE,block=NULL,color="classname",shape=NULL,
#'         labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'         labelSize=NULL,axisSize=NULL,titleSize=NULL)
#' 
#' # Loadings plot for individual components
#' # Separately for each block
#' p1 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
#'               combined=FALSE,block="expr",color="classname",shape=NULL,
#'               labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'               labelSize=NULL,axisSize=NULL,titleSize=NULL)
#' p2 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
#'               combined=FALSE,block="mirna",color="classname",shape=NULL,
#'               labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'               labelSize=NULL,axisSize=NULL,titleSize=NULL)
#' 
#' # Biplot: scores + loadings
#' plotRes(object=discoRes,comps=c(1,2),what="both",type="common",
#'         combined=FALSE,block="expr",color="classname",shape=NULL,
#'         labels=NULL,background=TRUE,palette=NULL,pointSize=4,
#'         labelSize=NULL,axisSize=NULL,titleSize=NULL)


setGeneric(
  name="plotRes",
  def=function(object,comps=c(1,2),what,type,combined,block=NULL,color=NULL,shape=NULL,labels=NULL,title=NULL,xlabel=NULL,ylabel=NULL,background=TRUE,palette=NULL,pointSize=4,
               labelSize=NULL,axisSize=NULL,titleSize=NULL,sizeValues=c(2,4),shapeValues=c(17,0)){
    standardGeneric("plotRes")}
)


setMethod(
  f="plotRes",
  signature=signature(
    object="caClass",
    comps="numeric",
    what="character",
    type="character",
    combined="logical"),
  definition=function(object,comps,what,type,combined,block,color,shape,labels,title,xlabel,ylabel,background,palette,pointSize,
                      labelSize,axisSize,titleSize,sizeValues,shapeValues){
    what <- match.arg(what,choices=c("scores","loadings","both"))
    type <- match.arg(type,choices=c("common","individual","both"))
    
    set.seed(12464863)
	
	#Combination of blocks can be done exclusibely for scores
    if (what!="scores" & combined==TRUE){
      warning("loadings are not able to be combined. Set combined=FALSE")
      combined <- FALSE
    }
    
    #Plotting loadings or both (scores + loadings) block definition is needed
    if (what!="scores" & is.null(block)){
      stop("'block' argument is required for loadings and loadings + score plot")
    }
	
    if (combined & object@caMethod%in%c("DISCO-SCA","JIVE") & type=="common" & what=="scores"){
     stop(paste(object@caMethod)," does not distinguish between common scores for each block. 
                    It is not possible to plot a combined plot in this case")
    }
    
    if (combined & object@caMethod=="O2PLS" & type!="both" & !is.null(block)){
      warning("'block' argument depreciated, if combined is TRUE both blocks are going to be used")
      block <- NULL
    }
    
    
    if (combined & object@caMethod%in%c("DISCO-SCA","JIVE") & type=="both" & what=="scores"){
      warning(paste(object@caMethod)," does not distinguish between common scores for each block")
    }
    
    #define auxiliar variable for ylabel when just 1 components is found for one
    common_aux <- FALSE
    dist1_aux <- FALSE
    dist2_aux <- FALSE

	
    #check and define type
    if (type=="common") type2 <- c("common","common")
    if (type=="individual") type2 <- c("dist","dist")
    if (type=="both") type2 <- c("common","dist")
    
    
    #check and define block
    if (is.null(block)){
      if (combined){
        block <- c("1","2")
        block_name <- object@Names
      }else if(type=="common" & object@caMethod%in%c("DISCO-SCA","JIVE")){
        block <- c("1","1")
        block_name <- c("","")
      }else{
        stop(block, " block argument must be provided")
      }
    }else{
      if (combined){
        if (type=="individual"){
          warning("'block' argument depreciated, if combined is TRUE both blocks are going to be used")
          block <- c("1","2")
          block_name <- object@Names
        }
        if (type=="both"){
          if (!block[1]%in%object@Names  & !block[1]%in%c("1","2")){
            stop("block ", block," does not exist")
          }
          if (block[1]%in%c("1","2")){
            bb1 <- which(c("1","2")==block)
            bb2 <- which(c("1","2")!=block)
            block_name <- c(object@Names[bb1], object@Names[bb2])
            block <- c(as.character(bb1),as.character(bb2))
          }
          if (block[1]%in%object@Names){
            bb1 <- which(object@Names==block)
            bb2 <- which(object@Names!=block)
            block_name <- c(object@Names[bb1], object@Names[bb2])
            block <- c(as.character(bb1),as.character(bb2))
          }
        }
      }else if (!block[1]%in%object@Names & !block[1]%in%c("1","2")){
        stop("block ", block," does not exist")
      }else{
        if (block[1]%in%c("1","2")){
          block <- c(block,block)
          bb1 <- which(c("1","2")==block)
          block_name <- c(object@Names[bb1], object@Names[bb1])
        }
        if (block[1]%in%object@Names){
          bb <- which(object@Names==block)
          block_name <- c(block, block)
          block <- c(as.character(bb),as.character(bb))
        }
      }
    }

    #matrix for plot (scores or loadings)
    if (what=="scores"){
      #generate a data frame with all the component data to plot
      if (object@caMethod%in%c("DISCO-SCA","JIVE")){
        df <- as.data.frame(cbind(object@scores$common,object@scores$common,object@scores$dist$Block1,object@scores$dist$Block2))
      }else{
        df <- as.data.frame(cbind(object@scores$common$Block1,object@scores$common$Block2,object@scores$dist$Block1,object@scores$dist$Block2))
      }
    }else if (what=="loadings"){
      nc1 <- matrix(NA,nrow=nrow(object@loadings$common$Block1),ncol=ncol(object@loadings$common$Block2))
      nc2 <- matrix(NA,nrow=nrow(object@loadings$common$Block2),ncol=ncol(object@loadings$common$Block1))
      nd1 <- matrix(NA,nrow=nrow(object@loadings$dist$Block1),ncol=ncol(object@loadings$dist$Block2))
      nd2 <- matrix(NA,nrow=nrow(object@loadings$dist$Block2),ncol=ncol(object@loadings$dist$Block1))
      nc <- nrow(nc1)+nrow(nc2)
      nd <- nrow(nd1)+nrow(nd2)
      nn1 <- matrix(NA,nrow=(nc-nd),ncol=ncol(object@loadings$dist$Block2))
      nn2 <- matrix(NA,nrow=(nc-nd),ncol=ncol(object@loadings$dist$Block1))
      df <- as.data.frame(cbind(rbind(object@loadings$common$Block1,nc2), rbind(object@loadings$common$Block2,nc1), rbind(object@loadings$dist$Block1,nd2,nn2), rbind(object@loadings$dist$Block2,nd1,nn1)))
      rownames(df) <- c(1:nrow(df))
    }else{ #both
      #scores
      if (object@caMethod%in%c("DISCO-SCA","JIVE")){
        dfs <- as.data.frame(cbind(object@scores$common,object@scores$common,object@scores$dist$Block1,object@scores$dist$Block2))
      }else{
        dfs <- as.data.frame(cbind(object@scores$common$Block1,object@scores$common$Block2,object@scores$dist$Block1,object@scores$dist$Block2))
      }
      #loadings
      nc1 <- matrix(NA,nrow=nrow(object@loadings$common$Block1),ncol=ncol(object@loadings$common$Block2))
      nc2 <- matrix(NA,nrow=nrow(object@loadings$common$Block2),ncol=ncol(object@loadings$common$Block1))
      nd1 <- matrix(NA,nrow=nrow(object@loadings$dist$Block1),ncol=ncol(object@loadings$dist$Block2))
      nd2 <- matrix(NA,nrow=nrow(object@loadings$dist$Block2),ncol=ncol(object@loadings$dist$Block1))
      nc <- nrow(nc1)+nrow(nc2)
      nd <- nrow(nd1)+nrow(nd2)
      nn1 <- matrix(NA,nrow=(nc-nd),ncol=ncol(object@loadings$dist$Block2))
      nn2 <- matrix(NA,nrow=(nc-nd),ncol=ncol(object@loadings$dist$Block1))
      dfl <- as.data.frame(cbind(rbind(object@loadings$common$Block1,nc2), rbind(object@loadings$common$Block2,nc1), rbind(object@loadings$dist$Block1,nd2,nn2), rbind(object@loadings$dist$Block2,nd1,nn1)))
      rownames(dfl) <- c(1:nrow(dfl))
      #combined dataframe
      df <- rbind(dfs,dfl)
      df$class <- c(rep("scores",nrow(dfs)),rep("loadings",nrow(dfl)))
    }
    if (what!="both"){
      colnames(df) <- c(paste("common1",1:object@commonComps,sep="_"),
                      paste("common2",1:object@commonComps,sep="_"),
                      paste("dist1",1:object@distComps[1],sep="_"),
                      paste("dist2",1:object@distComps[2],sep="_"))
    }else{
      colnames(df) <- c(paste("common1",1:object@commonComps,sep="_"),
                        paste("common2",1:object@commonComps,sep="_"),
                        paste("dist1",1:object@distComps[1],sep="_"),
                        paste("dist2",1:object@distComps[2],sep="_"),"class")      
    }
      
    #create auxiliar variables in case of just one component in common or distinctive structures
    if (object@commonComps==1){
      df$common1_2 <- runif(nrow(df),min(df$common1_1, na.rm=TRUE),max(df$common1_1, na.rm=TRUE))
      df$common2_2 <- runif(nrow(df),min(df$common2_1, na.rm=TRUE),max(df$common2_1, na.rm=TRUE))
      common_aux <- TRUE
    }
    if (object@distComps[1]==1){  #In Block 1
      df$dist1_2 <- runif(nrow(df),min(df$dist1_1, na.rm=TRUE),max(df$dist1_1, na.rm=TRUE))
      dist1_aux <- TRUE
    }     
    if (object@distComps[2]==1){  #In Block 2
      df$dist2_2 <- runif(nrow(df),min(df$dist2_1, na.rm=TRUE),max(df$dist2_1, na.rm=TRUE))
      dist2_aux <- TRUE
    }     

    if(type=="common" & common_aux){
      if(object@commonComps+1 < max(comps)) stop("Number of compsonents higher than allowed")
    }
    if(type=="common" & !common_aux){
      if(object@commonComps < max(comps)) stop("Number of compsonents higher than allowed")
    }
    
    
    if(type=="individual" & dist1_aux & block[1]=="1"){
      if(object@distComps[1]+1 < max(comps)) stop("Number of compsonents higher than allowed")
    }
    if(type=="individual" & !dist1_aux & block[1]=="1"){
      if(object@distComps[1] < max(comps)) stop("Number of compsonents higher than allowed")
    }
    
    if(type=="individual" & dist2_aux & block[1]=="2"){
      if(object@distComps[2]+1 < max(comps)) stop("Number of compsonents higher than allowed")
    }
    if(type=="individual" & !dist2_aux & block[1]=="2"){
      if(object@distComps[2] < max(comps)) stop("Number of compsonents higher than allowed")
    }
    
    
    if(type=="both" & common_aux){
      if(object@commonComps+1 < comps[1]) stop("Number of compsonents higher than allowed")
    }
    if(type=="both" & !common_aux){
      if(object@commonComps < comps[1]) stop("Number of compsonents higher than allowed")
    }
    if(type=="both" & dist1_aux & block[1]=="1"){
      if(object@distComps[1]+1 < comps[2]) stop("Number of compsonents higher than allowed")
    }
    if(type=="both" & !dist1_aux & block[1]=="1"){
      if(object@distComps[1] < comps[2]) stop("Number of compsonents higher than allowed")
    }
    if(type=="both" & dist2_aux & block[1]=="2"){
      if(object@distComps[2]+1 < comps[2]) stop("Number of compsonents higher than allowed")
    }
    if(type=="both" & !dist2_aux & block[1]=="2"){
      if(object@distComps[2] < comps[2]) stop("Number of compsonents higher than allowed")
    }
	
    #graphical parameters
    if (what=="scores"){     
      if (!is.null(color)){
        color <- match.arg(color,choices=colnames(pData(object@InitialData[[1]])))
      }else{
        color <- names(pData(object@InitialData[[1]]))[1]       
      }
      df <- data.frame(df,color=pData(object@InitialData[[1]])[,color])
      if (!is.null(shape)){
        shape <- match.arg(shape,choices=colnames(pData(object@InitialData[[1]])))
        df <- data.frame(df,shape=pData(object@InitialData[[1]])[,shape])
      }
      if (!is.null(labels)){
        labels <- match.arg(labels,choices=colnames(pData(object@InitialData[[1]])))
        df <- data.frame(df,labels=pData(object@InitialData[[1]])[,labels])
      }
    }else if (what=="loadings"){
      df <- data.frame(df,color="color")
      color <- "color"
    }else{ #both
      if (!is.null(color)){
        color <- match.arg(color,choices=colnames(pData(object@InitialData[[1]])))
      }else{
        color <- names(pData(object@InitialData[[1]]))[1]       
      }
      df <- data.frame(df,color=as.factor(c(as.character(pData(object@InitialData[[1]])[,color]),rep("color",nrow(dfl)))))
    }

    if(type=="common" & common_aux & sum(comps)>2){
      ylabel <- "Auxiliar to plot"
      warning("Only one common component found, auxiliar variable is defined to help visualization")
    }
    if(type=="individual" & (dist1_aux | dist2_aux)){
      if(block[1]=="1" & object@distComps[1]==1 & !combined){
        if(comps[1]==1) ylabel <- "Auxiliar to plot"
        if(comps[2]==1) xlabel <- "Auxiliar to plot"
        warning("Only one individual component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      } 
      if(block[1]=="2" & object@distComps[2]==1 & !combined){
        if(comps[1]==1) ylabel <- "Auxiliar to plot"
        if(comps[2]==1) xlabel <- "Auxiliar to plot"
        warning("Only one individual component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      }
      if(block[1]=="1" & object@distComps[1]==1 & combined & comps[1]>=2){
        xlabel <- "Auxiliar to plot"
        warning("Are the user sure to plot component ",comps[1]," for Block",block[1],"?\n Only one common component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      } 
      if(block[1]=="2" & object@distComps[2]==1 & combined & comps[2]>=2){
        ylabel <- "Auxiliar to plot"
        warning("Are the user sure to plot component ",comps[2]," for Block",block[1],"?\n Only one common component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      } 
    }
    if(type=="both" & (common_aux | dist1_aux | dist2_aux)){
      if(comps[1]>1 & common_aux){
        xlabel <- "Auxiliar to plot"
        warning("Only one common component found, auxiliar variable is defined to help visualization")
      }
      if(block[1]=="1" & dist1_aux & comps[2]>1){
        ylabel <- "Auxiliar to plot"
        warning("Only one individual component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      }
      if(block[1]=="2" & dist2_aux & comps[2]>1){
        ylabel <- "Auxiliar to plot"
        warning("Only one individual component found for Block",block[1],". Auxiliar variable is defined to help visualization.")
      }
    }
	
      #plot
      plotting(df,comps,xname=paste(type2[1],block[1],"_",comps[1],sep=""),yname=paste(type2[2],block[2],"_",comps[2],sep=""),color,what,type,combined,block=block_name,shape,labels,title,xlabel,
               ylabel,background,palette,pointSize,labelSize,axisSize,titleSize,sizeValues,shapeValues)  
  }      
)      
