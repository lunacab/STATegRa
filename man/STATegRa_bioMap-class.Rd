\name{bioMap-class}
\docType{class}
\alias{class:bioMap}
\alias{bioMap-class}
\title{
bioMap class
}
\description{
Class to manage mappings between omic features or identities.
}
\section{Constructor}{
  bioMap
}
\author{
David Gomez-Cabrero
}
\seealso{
\code{\link{bioMap-constructor}}
}
\examples{
# Definition of a bioMap Object

data(STATegRa_S2)

map.gene.miRNA<-bioMap(name = "Symbol-miRNA",
                metadata =  list(type_v1="Gene",type_v2="miRNA",
                                 source_database="targetscan.Hs.eg.db",
                                 data_extraction="July2014"),
                map=mapdata)
}
