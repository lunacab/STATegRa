data(STATegRa_S3)

test_PCAselection_single <- function() {
    checkEquals(PCA.selection(Data=Block2.PCA, fac.sel="single%", varthreshold=0.03)$numComps, 4)
}

test_PCAselection_accum <- function() {
    checkEquals(PCA.selection(Data=Block2.PCA, fac.sel="%accum", varthreshold=0.85)$numComps, 2)
}

test_PCAselection_relabs <- function() {
    checkEquals(PCA.selection(Data=Block2.PCA, fac.sel="rel.abs", nvar=0.5)$numComps, 4)
}

test_PCAselection_fixednum <- function() {    
    checkEquals(PCA.selection(Data=Block2.PCA, fac.sel="fixed.num", PCnum=3)$numComps, 3)
}