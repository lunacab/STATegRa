###########################################
########### EXAMPLE OF THE OMICSCLUSTERING
###########################################
require(STATegRa)

#############################################
## PART 1: CREATING a bioMap CLASS
#############################################
####### This part creates or reads the map between features.
####### In the present example the map is downloaded from a resource.
#######   then the class is created.

#load("../data/STATegRa_S2.rda")
data(STATegRa_S2)

MAP.SYMBOL<-bioMap(name = "Symbol-miRNA",
                metadata =  list(type_v1="Gene",type_v2="miRNA",
                                 source_database="targetscan.Hs.eg.db",
                                 data_extraction="July2014"),
                map=mapdata)


#############################################
## PART 2: CREATING a bioDist CLASS
#############################################
##### In the second part given a set of main features and surrogate feautres,
#####    the profile of the main features is computed through the surrogate features.

# Load Data
data(STATegRa_S1)
  #load("../data/STATegRa.S1.Rdata")

## Create ExpressionSets
#  source("../R/STATegRa_omicsPCA_classes_and_methods.R")
# Block1 - Expression data
mRNA.ds <- createOmicsExpressionSet(Data=Block1,pData=ed,pDataDescr=c("classname"))
# Block2 - miRNA expression data
miRNA.ds <- createOmicsExpressionSet(Data=Block2,pData=ed,pDataDescr=c("classname"))

# Create Gene-gene distance computed through miRNA data
bioDistmiRNA<-bioDist(referenceFeatures = rownames(Block1),     
             reference = "Var1",
             mapping = MAP.SYMBOL,
             surrogateData = miRNA.ds,  ### miRNA data
             referenceData = mRNA.ds,  ### mRNA data
             maxitems=2,
             selectionRule="sd",
             expfac=NULL,
             aggregation = "sum",
             distance = "spearman",
             noMappingDist = 0,
             filtering = NULL,
             name = "mRNAbymiRNA")

require(Biobase)

# Create Gene-gene distance through mRNA data
bioDistmRNA<-bioDistclass(name = "mRNAbymRNA",
                 distance = cor(t(exprs(mRNA.ds)),method="spearman"),
                 map.name = "id",
                 map.metadata = list(),
                 params = list())

#############################################
## PART 3: CREATING a LISTOF WEIGTHED DISTANCES MATRICES: bioDistWList
#############################################

bioDistList<-list(bioDistmRNA,bioDistmiRNA)
weights<-matrix(0,4,2)
weights[,1]<-c(0,0.33,0.67,1)
weights[,2]<-c(1,0.67,0.33,0)#

bioDistWList<-bioDistW(referenceFeatures = rownames(Block1),
                       bioDistList = bioDistList,
                       weights=weights)
length(bioDistWList)

#############################################
## PART 4: DEFINING THE STRENGTH OF ASSOCIATIONS IN GENERAL
#############################################

bioDistWPlot(referenceFeatures = rownames(Block1) ,
               listDistW = bioDistWList,
               method.cor="spearman")

#############################################
## PART 5: DEFINING THE ASSOCIATIONS FOR A GIVEN GENE
#############################################

## IDH1

IDH1.F<-bioDistFeature(Feature = "IDH1" ,
                       listDistW = bioDistWList,
                       threshold.cor=0.7)
bioDistFeaturePlot(data=IDH1.F)

## PDGFRA

#PDGFRA.F<-bioDistFeature(Feature = "PDGFRA" ,
#                       listDistW = bioDistWList,
#                       threshold.cor=0.7)
#bioDistFeaturePlot(data=PDGFRA.F,name="../vignettes/PDGFRA.png")

## EGFR
#EGFR.F<-bioDistFeature(Feature = "EGFR" ,
#                         listDistW = bioDistWList,
#                         threshold.cor=0.7)
#bioDistFeaturePlot(data=EGFR.F,name="../vignettes/EGFR.png")

## MGMT
#MGMT.F<-bioDistFeature(Feature = "MGMT" ,
#                         listDistW = bioDistWList,
#                         threshold.cor=0.5)
#bioDistFeaturePlot(data=MGMT.F,name="../vignettes/MGMT.png")




