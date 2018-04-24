###########################################
########### EXAMPLE OF THE OMICSPCA
###########################################
require(STATegRa)

# g_legend (not exported by STATegRa any more)
## code from https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

#########################
## PART 1. Load data

## Load data
data(STATegRa_S3)

ls()

## Create ExpressionSets
# Block1 - Expression data
B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,pDataDescr=c("classname"))
# Block2 - miRNA expression data
B2 <- createOmicsExpressionSet(Data=Block2.PCA,pData=ed.PCA,pDataDescr=c("classname"))

#########################
## PART 2.  Model Selection

require(grid)
require(gridExtra)
require(ggplot2)

## Select the optimal components
ms <- modelSelection(Input=list(B1,B2),Rmax=4,fac.sel="single%",varthreshold=0.03,center=TRUE,scale=TRUE,weight=TRUE)


#########################
## PART 3. Component Analysis

## 3.1 Component analysis of the three methods
discoRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),method="DISCOSCA",Rcommon=2,Rspecific=c(2,2),center=TRUE,
                              scale=TRUE,weight=TRUE)
jiveRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),method="JIVE",Rcommon=2,Rspecific=c(2,2),center=TRUE,
                             scale=TRUE,weight=TRUE)
o2plsRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),method="O2PLS",Rcommon=2,Rspecific=c(2,2),center=TRUE,
                              scale=TRUE,weight=TRUE)

## 3.2 Exploring scores structures

# Exploring DISCO-SCA scores structure
discoRes@scores$common ## Common scores
discoRes@scores$dist[[1]] ## Distinctive scores for Block 1
discoRes@scores$dist[[2]] ## Distinctive scores for Block 2
# Exploring O2PLS scores structure
o2plsRes@scores$common[[1]] ## Common scores for Block 1
o2plsRes@scores$common[[2]] ## Common scores for Block 2
o2plsRes@scores$dist[[1]] ## Distinctive scores for Block 1
o2plsRes@scores$dist[[2]] ## Distinctive scores for Block 2

## 3.3 Plotting VAF

# DISCO-SCA plotVAF
plotVAF(discoRes)

# JIVE plotVAF
plotVAF(jiveRes)


#########################
## PART 4. Plot Results

# Scores for common part. DISCO-SCA
plotRes(object=discoRes,comps=c(1,2),what="scores",type="common",
             combined=FALSE,block=NULL,color="classname",shape=NULL,labels=NULL,
             background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
             axisSize=NULL,titleSize=NULL)

# Scores for common part. JIVE
plotRes(object=jiveRes,comps=c(1,2),what="scores",type="common",
             combined=FALSE,block=NULL,color="classname",shape=NULL,labels=NULL,
             background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
             axisSize=NULL,titleSize=NULL)

# Scores for common part. O2PLS.
p1 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=o2plsRes,comps=c(1,2),what="scores",type="common",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
legend <- g_legend(p1)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             legend,heights=c(6/7,1/7))

# Combined plot of scores for common part. O2PLS.
plotRes(object=o2plsRes,comps=c(1,1),what="scores",type="common",
             combined=TRUE,block=NULL,color="classname",shape=NULL,
             labels=NULL,background=TRUE,palette=NULL,pointSize=4,
             labelSize=NULL,axisSize=NULL,titleSize=NULL)


# Scores for distinctive part. DISCO-SCA. (two plots one for each block)
p1 <- plotRes(object=discoRes,comps=c(1,2),what="scores",type="individual",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="scores",type="individual",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
legend <- g_legend(p1)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             legend,heights=c(6/7,1/7))

# Combined plot of scores for distinctive part. DISCO-SCA
plotRes(object=discoRes,comps=c(1,1),what="scores",type="individual",
             combined=TRUE,block=NULL,color="classname",shape=NULL,
             labels=NULL,background=TRUE,palette=NULL,pointSize=4,
             labelSize=NULL,axisSize=NULL,titleSize=NULL)

# Combined plot of scores for common and distinctive part. O2PLS (two plots one for each block)
p1 <- plotRes(object=o2plsRes,comps=c(1,1),what="scores",type="both",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=o2plsRes,comps=c(1,1),what="scores",type="both",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
legend <- g_legend(p1)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             legend,heights=c(6/7,1/7))

# Combined plot of scores for common and distinctive part. DISCO  (two plots one for each block)
p1 <- plotRes(object=discoRes,comps=c(1,1),what="scores",type="both",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,1),what="scores",type="both",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
legend <- g_legend(p1)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             legend,heights=c(6/7,1/7))

# Loadings for common part. DISCO-SCA. (two plots one for each block)
p1 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="common",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="common",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))


# Loadings for distinctive part. DISCO-SCA. (two plots one for each block)
p1 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="loadings",type="individual",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))


# Combined plot for loadings from common and distinctive part  (two plots one for each block)
p1 <- plotRes(object=discoRes,comps=c(1,1),what="loadings",type="both",
              combined=FALSE,block="expr",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,1),what="loadings",type="both",
              combined=FALSE,block="mirna",color="classname",shape=NULL,
              labels=NULL,background=TRUE,palette=NULL,pointSize=4,
              labelSize=NULL,axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))



## Plot scores and loadings togheter: Common components DISCO-SCA
p1 <- plotRes(object=discoRes,comps=c(1,2),what="both",type="common",
        combined=FALSE,block="expr",color="classname",shape=NULL,labels=NULL,
        background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
        axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="both",type="common",
              combined=FALSE,block="mirna",color="classname",shape=NULL,labels=NULL,
              background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
              axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))


## Plot scores and loadings togheter:  Common components O2PLS
p1 <- plotRes(object=o2plsRes,comps=c(1,2),what="both",type="common",
              combined=FALSE,block="expr",color="classname",shape=NULL,labels=NULL,
              background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
              axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=o2plsRes,comps=c(1,2),what="both",type="common",
              combined=FALSE,block="mirna",color="classname",shape=NULL,labels=NULL,
              background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
              axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))


## Plot scores and loadings togheter: Distintive components DISCO-SCA
p1 <- plotRes(object=discoRes,comps=c(1,2),what="both",type="individual",
              combined=FALSE,block="expr",color="classname",shape=NULL,labels=NULL,
              background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
              axisSize=NULL,titleSize=NULL)
p2 <- plotRes(object=discoRes,comps=c(1,2),what="both",type="individual",
              combined=FALSE,block="mirna",color="classname",shape=NULL,labels=NULL,
              background=TRUE,palette=NULL,pointSize=4,labelSize=NULL,
              axisSize=NULL,titleSize=NULL)
grid.arrange(arrangeGrob(p1+theme(legend.position="none"),
                         p2+theme(legend.position="none"),nrow=1),
             heights=c(6/7,1/7))



