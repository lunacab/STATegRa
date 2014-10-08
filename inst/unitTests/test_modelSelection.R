data(STATegRa_S3)

B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA, pDataDescr=c("classname"))
B2 <- createOmicsExpressionSet(Data=Block2.PCA, pData=ed.PCA, pDataDescr=c("classname"))

test_modelSelection_single <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="single%", varthreshold=0.03)
    checkEquals(res$common, 2)
    checkEquals(res$dist, c(2, 2))
}

test_modelSelection_accum <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="%accum", varthreshold=0.9)
    checkEquals(res$common, 2)
    checkEquals(res$dist, c(0, 0))
}

#test_modelSelection_relabs <- function() {
    # doesn't work
    #res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="rel.abs", nvar=0.5)
#}

test_modelSelection_fixednum <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="fixed.num", PCnum=3)
    checkEquals(res$common, 2)
    checkEquals(res$dist, c(1, 1))
}
