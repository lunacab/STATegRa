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
                theme(legend.title=element_blank(),plot.title=element_text(size=18),legend.position="bottom")
            p_dist <- ggplot(df_dist,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("")+
                ggtitle("Distinctive components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=18),legend.position="bottom")
            mylegend <- g_legend(p_common)
            grid.arrange(arrangeGrob(p_common+ theme(legend.position="none"),p_dist+ theme(legend.position="none"),ncol=2,widths=c(1/3,2/3)),
                         mylegend,nrow=2,heights=c(10, 1),main=textGrob(label=mainTitle,gp=gpar(fontsize=20,cex=2)))
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
                theme(legend.title=element_blank(),plot.title=element_text(size=18),legend.position="bottom")
            p_dist <- ggplot(df_dist,aes(x=factor(comp),y=VAF,fill=factor(block)))+
                geom_bar(stat="identity",position="dodge")+
                facet_wrap(~ group)+
                xlab("Components")+
                ylab("")+
                ggtitle("Distinctive components")+
                ylim(c(0,upperlim))+
                theme(legend.title=element_blank(),plot.title=element_text(size=18),legend.position="bottom")
            grid.arrange(arrangeGrob(p_common+ theme(legend.position="none"),p_dist+ theme(legend.position="none"),ncol=2,widths=c(1/3,2/3)),
                         nrow=2,heights=c(10, 1),main=textGrob(label=mainTitle,gp=gpar(fontsize=20,cex=2)))
        } else if(object@caMethod=="O2PLS"){
            warning("VAF per component cannot be calculated with O2PLS approach because components are not orthogonal")
        }
    }
)

# plotting ----------------------------------------------------------------
setGeneric(
    name="plotting",
    def=function(df,xname,yname,color,shape=NULL,labels=NULL,title="",pointsize=4,xlabel=NULL,ylabel=NULL){standardGeneric("plotting")}
)

setMethod(
    f="plotting",
    signature=signature(
        df="data.frame",
        xname="character",
        yname="character",
        color="character"),
    definition=function(df,xname,yname,color,shape,labels,title,pointsize,xlabel,ylabel){
        xlim <- ylim <- c(max(df[,xname],df[,yname]),min(df[,xname],df[,yname]))
        p <- ggplot(df,aes_string(x=xname,y=yname))
        if (!is.null(shape)){
            p <- p+geom_point(aes_string(color=color,shape=shape),size=pointsize)
        }else{
            p <- p+geom_point(aes_string(color=color),size=pointsize)
        }
        if (is.null(xlabel) & is.null(ylabel)){
            p <- p+xlab("First component")+ylab("Second component")
        } else{
            p <- p+xlab(xlabel)+ylab(ylabel)
        }
        p <- p+xlim(xlim)+
            ylim(ylim)+
            ggtitle(title)+
            coord_fixed()+
            theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="bottom")
        if(!is.null(labels)){
            p <- p+geom_text(alpha=.5,size=pointsize,aes_string(label=labels),vjust=-.7)
        }
        return(p)
    }
)

# biplotRes (Associated to caClass) ------------------------------------------------------------------
#   Parameters:
#   type=c("common","individual"): plot common or individual components
#   comps=c(comp_x,comp_y): components to plot. Indicate x and y, if combined=T, x=common_comp and y=dist_comp
#   block=c(1,2): Block to plot, if common=T ignored except for O2PLS
# Default values
# comps <- c(1,2)
# sizeValues <- c(2,4)
# colorCol <- "classname"
# shapeValues <- c(17,0)
setGeneric(
    name="biplotRes",
    def=function(object,type,comps,block,title=NULL,colorCol=NULL,sizeValues=c(2,4),shapeValues=c(17,0),background=TRUE,pointSize=4,
                 labelSize=NULL,axisSize=NULL,titleSize=NULL){standardGeneric("biplotRes")}
)

setMethod(
    f="biplotRes",
    signature=signature(
        object="caClass",
        type="character",
        comps="numeric",
        block="character"),
    definition=function(object,type,comps,block,title,colorCol,sizeValues,shapeValues,background,pointSize,labelSize,axisSize,
                        titleSize){
        if (is.null(colorCol)){
            colorCol <- colnames(pData(object@InitialData[[1]]))[1]
        }
        if (type=="common"){
            if (object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                sc.df <- data.frame(values=rbind(object@scores$common[,comps],object@scores$common[,comps]),
                                    class="scores",
                                    block=c(rep("B1",nrow(object@scores$common)),rep("B2",nrow(object@scores$common))),
                                    color=rep(unlist(pData(object@InitialData[[1]])[,colorCol]),2),row.names=NULL)
            }else{
                sc.df <- data.frame(values=rbind(object@scores$common[[1]][,comps],object@scores$common[[2]][,comps]),
                                    class="scores",
                                    block=c(rep("B1",nrow(object@scores$common[[1]])),rep("B2",nrow(object@scores$common[[2]]))),
                                    color=rep(unlist(pData(object@InitialData[[1]])[,colorCol]),2))
            }
            ld.df <- data.frame(values=rbind(object@loadings$common[[1]][,comps],object@loadings$common[[2]][,comps]),
                                class="loadings",
                                block=c(rep("B1",nrow(exprs(object@InitialData[[1]]))),rep("B2",nrow(exprs(object@InitialData[[2]])))),
                                color=c(rep("B1",nrow(exprs(object@InitialData[[1]]))),rep("B2",nrow(exprs(object@InitialData[[2]])))))
            df <- rbind(sc.df,ld.df)
            if(object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                ## In this case information in B1 is the same than in B2
                if (is.null(title)){
                    title <- paste("Biplot associated to common components of",object@Names[1])
                }
                p <- ggplot(data=df,aes(x=values.1,y=values.2,shape=class,size=class))+
                    geom_point(data=df[df$class=="scores",],aes(color=color))+
                    geom_point(data=df[df$class=="loadings",])+
                    scale_size_manual(name="",values=sizeValues)+
                    scale_shape_manual(name="",values=shapeValues)+
                    theme(legend.title=element_blank(),legend.position="bottom")+
                    xlab(paste("Component",comps[1]))+
                    ylab(paste("Component",comps[2]))+
                    ggtitle(title)
            }else{
                if (block==1 | block==object@Names[1]){
                    BB <- df[df$block=="B1",]
                    if (is.null(title)){
                        title <- paste("Biplot associated to common components of",object@Names[1])
                    }
                    p <- ggplot(data=BB,aes(x=values.1,y=values.2,shape=class,size=class))+
                        geom_point(data=BB[BB$class=="scores",],aes(color=color))+
                        geom_point(data=BB[BB$class=="loadings",])+
                        scale_size_manual(name="",values=sizeValues)+
                        scale_shape_manual(name="",values=shapeValues)+
                        theme(legend.title=element_blank())+
                        xlab(paste("Component",comps[1]))+
                        ylab(paste("Component",comps[2]))+
                        ggtitle(title)
                } else if (block==2 | block==object@Names[2]){
                    BB <- df[df$block=="B2",]
                    if (is.null(title)){
                        title <- paste("Biplot associated to common components of",object@Names[2])
                    }
                    p <- ggplot(data=BB,aes(x=values.1,y=values.2,shape=class,size=class))+
                        geom_point(data=BB[BB$class=="scores",],aes(color=color))+
                        geom_point(data=BB[BB$class=="loadings",])+
                        scale_size_manual(name="",values=sizeValues)+
                        scale_shape_manual(name="",values=shapeValues)+
                        theme(legend.title=element_blank())+
                        xlab(paste("Component",comps[1]))+
                        ylab(paste("Component",comps[2]))+
                        ggtitle(title)
                }else{
                    warning(paste("block=",block,"is not allowed."))
                }
            }
        } else{ ## type="individual"
            sc.df <- data.frame(values=rbind(object@scores$dist[[1]][,comps],object@scores$dist[[2]][,comps],deparse.level=0),
                                class="scores",
                                block=c(rep("B1",nrow(object@scores$common)),rep("B2",nrow(object@scores$common))),
                                color=rep(unlist(pData(object@InitialData[[1]])[,colorCol]),2))
            ld.df <- data.frame(values=rbind(object@loadings$dist[[1]][,comps],object@loadings$dist[[2]][,comps]),
                                class="loadings",
                                block=c(rep("B1",nrow(exprs(object@InitialData[[1]]))),rep("B2",nrow(exprs(object@InitialData[[2]])))),
                                color=c(rep("B1",nrow(exprs(object@InitialData[[1]]))),rep("B2",nrow(exprs(object@InitialData[[2]])))))
            df <- rbind(sc.df,ld.df)
            if (block==1 | block==object@Names[1]){
                BB <- df[df$block=="B1",]
                if (is.null(title)){
                    title <- paste("Biplot associated to distinctive components of",object@Names[1])
                }
                p <- ggplot(data=BB,aes(x=values.1,y=values.2,shape=class,size=class))+
                    geom_point(data=BB[BB$class=="scores",],aes(color=color))+
                    geom_point(data=BB[BB$class=="loadings",])+
                    scale_size_manual(name="",values=sizeValues)+
                    scale_shape_manual(name="",values=shapeValues)+
                    theme(legend.title=element_blank())+
                    xlab(paste("Component",comps[1]))+
                    ylab(paste("Component",comps[2]))+
                    ggtitle(title)
            } else if (block==2 | block==object@Names[2]){
                BB <- df[df$block=="B2",]
                if (is.null(title)){
                    title <- paste("Biplot associated to distinctive components of",object@Names[2])
                }
                p <- ggplot(data=BB,aes(x=values.1,y=values.2,shape=class,size=class))+
                    geom_point(data=BB[BB$class=="scores",],aes(color=color))+
                    geom_point(data=BB[BB$class=="loadings",])+
                    scale_size_manual(name="",values=sizeValues)+
                    scale_shape_manual(name="",values=shapeValues)+
                    theme(legend.title=element_blank())+
                    xlab(paste("Component",comps[1]))+
                    ylab(paste("Component",comps[2]))+
                    ggtitle(title)
            }else{
                warning(paste("block=",block,"is not allowed."))
            }
        }
        ## Adding the other graphical parameters
        if (!is.null(p)){
            if (!background){
                p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank())
            }
            if (!is.null(labelSize)){
                p <- p+theme(axis.title=element_text(size=labelSize))
            }
            if (!is.null(axisSize)){
                p <- p+theme(axis.text=element_text(size=axisSize))
            }
            if (!is.null(titleSize)){
                p <- p+theme(plot.title = element_text(size=titleSize))
            }
            return(p)
        }
    }
)

# plotRes (Associated to caClass) --------------------------------------------
# Parameters:
#   what=c("scores","loadings"): plot scores or loadings
#   type=c("common","individual","both"): plot common or individual components (both for combined plots)
#   combined=c(T,F): Combine common and distinctive components?
#   comps=c(comp_x,comp_y): components to plot. Indicate x and y, if combined=T, x=common_comp and y=dist_comp
#   block=c(1,2): Block to plot, if common=T ignored except for O2PLS
setGeneric(
    name="plotRes",
    def=function(object,comps=c(1,2),what,type,combined,block,color=NULL,shape=NULL,labels=NULL,background=TRUE,palette=NULL,pointSize=4,
                 labelSize=NULL,axisSize=NULL,titleSize=NULL){
        standardGeneric("plotRes")}
)
setMethod(
    f="plotRes",
    signature=signature(
        object="caClass",
        comps="numeric",
        what="character",
        type="character",
        combined="logical",
        block="character"),
    definition=function(object,comps,what,type,combined,block,color,shape,labels,background,palette,pointSize,
                        labelSize,axisSize,titleSize){
        p <- NULL
        what <- match.arg(what,choices=c("scores","loadings"))
        type <- match.arg(type,choices=c("common","individual","both"))
        if(what=="scores"){
            ## Creating data.frame
            if(object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                if (object@commonComps==1){
                    if (object@distComps[1]==1 & object@distComps[2]==2){
                        scores.df <- data.frame(object@scores$common[,comps[1]],rep(0,nrow(object@scores$common)),object@scores$dist[[1]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[1]])),object@scores$dist[[2]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else if (object@distComps[1]==1){
                        scores.df <- data.frame(object@scores$common[,comps[1]],rep(0,nrow(object@scores$common)),object@scores$dist[[1]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[1]])),object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }else if (object@distComps[2]==1){
                        scores.df <- data.frame(object@scores$common[,comps[1]],rep(0,nrow(object@scores$common)),object@scores$dist[[1]][,comps],
                                                object@scores$dist[[2]][,comps[1]],rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else{
                        scores.df <- data.frame(object@scores$common[,comps[1]],rep(0,nrow(object@scores$common)),object@scores$dist[[1]][,comps],
                                                object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }
                }else{
                    if (object@distComps[1]==1 & object@distComps[2]==2){
                        scores.df <- data.frame(object@scores$common[,comps],object@scores$dist[[1]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[1]])),object@scores$dist[[2]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else if (object@distComps[1]==1){
                        scores.df <- data.frame(object@scores$common[,comps],object@scores$dist[[1]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[1]])),object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }else if (object@distComps[2]==1){
                        scores.df <- data.frame(object@scores$common[,comps],object@scores$dist[[1]][,comps],
                                                object@scores$dist[[2]][,comps[1]],rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else{
                        scores.df <- data.frame(object@scores$common[,comps],object@scores$dist[[1]][,comps],object@scores$dist[[2]][,comps],
                                                pData(object@InitialData[[1]]))
                    }
                }
                names(scores.df) <- c(paste("comm",1:2,sep="_"),paste("distB1",1:2,sep="_"),paste("distB2",1:2,sep="_"),
                                      names(pData(object@InitialData[[1]])))
            } else{ ## object@caMethod=="O2PLS"
                if (object@commonComps==1){
                    if (object@distComps[1]==1 & object@distComps[2]==2){
                        scores.df <- data.frame(object@scores$common[[1]][,comps[1]],rep(0,nrow(object@scores$common[[1]])),
                                                object@scores$common[[2]][,comps[1]],rep(0,nrow(object@scores$common[[2]])),
                                                object@scores$dist[[1]][,comps[1]],rep(0,nrow(object@scores$dist[[1]])),
                                                object@scores$dist[[2]][,comps[1]],rep(0,nrow(object@scores$dist[[2]])),
                                                pData(object@InitialData[[1]]))
                    }else if (object@distComps[1]==1){
                        scores.df <- data.frame(object@scores$common[[1]][,comps[1]],rep(0,nrow(object@scores$common[[1]])),
                                                object@scores$common[[2]][,comps[1]],rep(0,nrow(object@scores$common[[2]])),
                                                object@scores$dist[[1]][,comps[1]],rep(0,nrow(object@scores$dist[[1]])),
                                                object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }else if (object@distComps[2]==1){
                        scores.df <- data.frame(object@scores$common[[1]][,comps[1]],rep(0,nrow(object@scores$common[[1]])),
                                                object@scores$common[[2]][,comps[1]],rep(0,nrow(object@scores$common[[2]])),
                                                object@scores$dist[[1]][,comps],object@scores$dist[[2]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else{
                        scores.df <- data.frame(object@scores$common[[1]][,comps[1]],rep(0,nrow(object@scores$common[[1]])),
                                                object@scores$common[[2]][,comps[1]],rep(0,nrow(object@scores$common[[2]])),
                                                object@scores$dist[[1]][,comps],object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }
                }else{
                    if (object@distComps[1]==1 & object@distComps[2]==2){
                        scores.df <- data.frame(object@scores$common[[1]][,comps],object@scores$common[[2]][,comps],
                                                object@scores$dist[[1]][,comps[1]],rep(0,nrow(object@scores$dist[[1]])),
                                                object@scores$dist[[2]][,comps[1]],rep(0,nrow(object@scores$dist[[2]])),
                                                pData(object@InitialData[[1]]))
                    }else if (object@distComps[1]==1){
                        scores.df <- data.frame(object@scores$common[[1]][,comps],object@scores$common[[2]][,comps],
                                                object@scores$dist[[1]][,comps[1]],rep(0,nrow(object@scores$dist[[1]])),
                                                object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }else if (object@distComps[2]==1){
                        scores.df <- data.frame(object@scores$common[[1]][,comps],object@scores$common[[2]][,comps],
                                                object@scores$dist[[1]][,comps],object@scores$dist[[2]][,comps[1]],
                                                rep(0,nrow(object@scores$dist[[2]])),pData(object@InitialData[[1]]))
                    }else{
                        scores.df <- data.frame(object@scores$common[[1]][,comps],object@scores$common[[2]][,comps],object@scores$dist[[1]][,comps],
                                                object@scores$dist[[2]][,comps],pData(object@InitialData[[1]]))
                    }
                }
                names(scores.df) <- c(paste("commB1",1:2,sep="_"),paste("commB2",1:2,sep="_"),paste("distB1",1:2,sep="_"),
                                      paste("distB2",1:2,sep="_"),names(pData(object@InitialData[[1]])))
            }
            if (is.null(color)){
                color <- names(pData(object@InitialData[[1]]))[1]
            }
            ## Evaluating parameters
            if (combined){
                if (type=="common"){
                    if (object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                        warning(paste("It is not possible to combine common components in",object@caMethod,"approach."))
                    } else{ ## object@caMethod=="O2PLS"
                        p <- plotting(df=scores.df,xname="commB1_1",yname="commB2_2",color=color,
                                      title=paste(object@caMethod,"common scores:",object@Names[1],"vs",object@Names[2]),
                                      xlabel=paste("Component",comps[1],"of",object@Names[1]),
                                      ylabel=paste("Component",comps[2],"of",object@Names[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)
                    }
                } else if (type=="individual"){
                    p <- plotting(df=scores.df,xname="distB1_1",yname="distB2_2",color=color,
                                  title=paste(object@caMethod,"distinctive scores:",object@Names[1],"vs",object@Names[2]),
                                  xlabel=paste("Component",comps[1],"of",object@Names[1]),
                                  ylabel=paste("Component",comps[2],"of",object@Names[2]),
                                  pointsize=pointSize,shape=shape,labels=labels)
                } else{
                    if(block==1 | block==object@Names[1]){
                        if (object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                            p <- plotting(df=scores.df,xname="comm_1",yname="distB1_2",color=color,
                                          title=paste(object@caMethod,"common and distinctive scores:",object@Names[1]),
                                          xlabel=paste("Common component",comps[1]),
                                          ylabel=paste("Distinctive component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }else{ ## object@caMethod=="O2PLS"
                            p <- plotting(df=scores.df,xname="commB1_1",yname="distB1_2",color=color,
                                          title=paste(object@caMethod,"common and distinctive scores:",object@Names[1]),
                                          xlabel=paste("Common component",comps[1]),
                                          ylabel=paste("Distinctive component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }
                    }else if (block==2 | block==object@Names[2]){
                        if (object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                            p <- plotting(df=scores.df,xname="comm_2",yname="distB2_2",color=color,
                                          title=paste(object@caMethod,"common and distinctive scores:",object@Names[2]),
                                          xlabel=paste("Common component",comps[1]),
                                          ylabel=paste("Distinctive component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }else{ ## object@caMethod=="O2PLS"
                            p <- plotting(df=scores.df,xname="commB2_1",yname="distB2_2",color=color,
                                          title=paste(object@caMethod,"common and distinctive scores:",object@Names[2]),
                                          xlabel=paste("Common component",comps[1]),
                                          ylabel=paste("Distinctive component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }
                    }else{
                        warning(paste("block=",block,"is not allowed."))
                    }
                }
            } else{ ## combined=FALSE
                if (type=="common"){
                    if(object@caMethod=="DISCO-SCA" | object@caMethod=="JIVE"){
                        p <- plotting(df=scores.df,xname="comm_1",yname="comm_2",color=color,
                                      title=paste(object@caMethod,"common Scores"),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)
                    }else{ ## object@caMethod=="O2PLS"
                        if(block==1 | block==object@Names[1]){
                            p <- plotting(df=scores.df,xname="commB1_1",yname="commB1_2",color=color,
                                          title=paste(object@caMethod,"common scores:",object@Names[1]),
                                          xlabel=paste("Component",comps[1]),
                                          ylabel=paste("Component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }else if (block==2 | block==object@Names[2]){
                            p <- plotting(df=scores.df,xname="commB2_1",yname="commB2_2",color=color,
                                          title=paste(object@caMethod,"common scores:",object@Names[2]),
                                          xlabel=paste("Component",comps[1]),
                                          ylabel=paste("Component",comps[2]),
                                          pointsize=pointSize,shape=shape,labels=labels)
                        }else{
                            warning(paste("block=",block,"is not allowed."))
                        }
                    }
                } else if (type=="individual"){
                    if(block==1 | block==object@Names[1]){
                        p <- plotting(df=scores.df,xname="distB1_1",yname="distB1_2",color=color,
                                      title=paste(object@caMethod,"distinctive scores:",object@Names[1]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)
                    }else if (block==2 | block==object@Names[2]){
                        p <- plotting(df=scores.df,xname="distB2_1",yname="distB2_2",color=color,
                                      title=paste(object@caMethod,"distinctive scores:",object@Names[2]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)
                    }else{
                        warning(paste("block=",block,"is not allowed. "))
                    }
                }else{ ## type=="both"
                    warning ("type=\"both\" is only allowed in combined plots (combined=TRUE)")
                }
            }
        } else{ ## loadings
            loadings.df <- data.frame(rbind(object@loadings$common[[1]][,comps],object@loadings$common[[2]][,comps]),
                                      rbind(object@loadings$dist[[1]][,comps],object@loadings$dist[[2]][,comps]),
                                      color=c(rep("B1",nrow(exprs(object@InitialData[[1]]))),rep("B2",nrow(exprs(object@InitialData[[2]])))))
            names(loadings.df) <- c(paste("comm",1:2,sep="_"),paste("dist",1:2,sep="_"),"color")
            ## Evaluating parameters
            if (combined){
                if (type=="common"){
                    ## Scaling
                    n1 <- nrow(exprs(object@InitialData[[1]]))
                    n2 <- nrow(exprs(object@InitialData[[2]]))
                    if (n1>n2){
                        loadings.df[loadings.df$color=="B2",1:4] <- loadings.df[loadings.df$color=="B2",1:4]*sqrt(n1/n2)
                    }else if (n1<n2){
                        loadings.df[loadings.df$color=="B1",1:4] <- loadings.df[loadings.df$color=="B1",1:4]*sqrt(n2/n1)
                    }
                    p <- plotting(df=loadings.df,xname="comm_1",yname="comm_2",color="color",
                                  title=paste(object@caMethod,"common loadings",object@Names[1],"vs",object@Names[2]),
                                  xlabel=paste("Component",comps[1]),
                                  ylabel=paste("Component",comps[2]),
                                  pointsize=pointSize,shape=shape,labels=labels)
                } else if (type=="individual"){
                    ## Scaling
                    n1 <- nrow(exprs(object@InitialData[[1]]))
                    n2 <- nrow(exprs(object@InitialData[[2]]))
                    if (n1>n2){
                        loadings.df[loadings.df$color=="B2",1:4] <- loadings.df[loadings.df$color=="B2",1:4]*sqrt(n1/n2)
                    }else if (n1<n2){
                        loadings.df[loadings.df$color=="B1",1:4] <- loadings.df[loadings.df$color=="B1",1:4]*sqrt(n2/n1)
                    }
                    p <- plotting(df=loadings.df,xname="dist_1",yname="dist_2",color="color",
                                  title=paste(object@caMethod,"distinctive loadings",object@Names[1],"vs",object@Names[2]),
                                  xlabel=paste("Component",comps[1]),
                                  ylabel=paste("Component",comps[2]),
                                  pointsize=pointSize,shape=shape,labels=labels)
                } else{ ## type="both"
                    if (block==1 | block==object@Names[1]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B1",],xname="comm_1",yname="dist_2",color="color",
                                      title=paste(object@caMethod,"common and distinctive loadings:",object@Names[1]),
                                      xlabel=paste("common component",comps[1]),
                                      ylabel=paste("distinctive component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else if (block==2 | block==object@Names[2]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B2",],xname="comm_1",yname="dist_2",color="color",
                                      title=paste(object@caMethod,"common and distinctive loadings:",object@Names[2]),
                                      xlabel=paste("common component",comps[1]),
                                      ylabel=paste("distinctive component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else{
                        warning(paste("block=",block,"is not allowed."))
                    }
                }
            } else{ ## combined=FALSE
                if (type=="common"){
                    if (block==1 | block==object@Names[1]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B1",],xname="comm_1",yname="comm_2",color="color",
                                      title=paste(object@caMethod,"common loadings:",object@Names[1]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else if (block==2 | block==object@Names[2]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B2",],xname="comm_1",yname="comm_2",color="color",
                                      title=paste(object@caMethod,"common loadings:",object@Names[2]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else{
                        warning(paste("block=",block,"is not allowed."))
                    }
                } else if (type=="individual"){
                    if (block==1 | block==object@Names[1]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B1",],xname="dist_1",yname="dist_2",color="color",
                                      title=paste(object@caMethod,"distinctive loadings:",object@Names[1]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else if (block==2 | block==object@Names[2]){
                        p <- plotting(df=loadings.df[loadings.df$color=="B2",],xname="dist_1",yname="dist_2",color="color",
                                      title=paste(object@caMethod,"distinctive loadings:",object@Names[2]),
                                      xlabel=paste("Component",comps[1]),
                                      ylabel=paste("Component",comps[2]),
                                      pointsize=pointSize,shape=shape,labels=labels)+
                            theme(legend.position="none")
                    } else{
                        warning(paste("block=",block,"is not allowed."))
                    }
                } else{
                    warning ("type=\"both\" is only allowed in combined plots (combined=TRUE)")
                }
            }
        }
        ## Adding the other graphical parameters
        if (!is.null(p)){
            if (!background){
                p <- p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.background = element_blank())
            }
            if (!is.null(palette)){
                p <- p+scale_colour_manual(values=palette)
            }
            if (!is.null(labelSize)){
                p <- p+theme(axis.title=element_text(size=labelSize))
            }
            if (!is.null(axisSize)){
                p <- p+theme(axis.text=element_text(size=axisSize))
            }
            if (!is.null(titleSize)){
                p <- p+theme(plot.title = element_text(size=titleSize))
            }
            return(p)
        }
    }
)