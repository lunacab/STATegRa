data(STATegRa_S3)

B1 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA, pDataDescr=c("classname"))
B2 <- createOmicsExpressionSet(Data=Block2.PCA, pData=ed.PCA, pDataDescr=c("classname"))

test_modelSelection_single <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="single%", varthreshold=0.03, center=TRUE, scale=FALSE, weight=TRUE, plot_common=FALSE, plot_dist=FALSE)
    checkEquals(res$common$commonComps, 2)
    checkEquals(res$dist[[1]]$numComps, 2)
    checkEquals(res$dist[[2]]$numComps, 2)
}


test_modelSelection_accum <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="%accum", varthreshold=0.9, center=TRUE, scale=FALSE, weight=TRUE, plot_common=FALSE, plot_dist=FALSE)
    checkEquals(res$common$commonComps, 2)
    checkEquals(res$dist[[1]]$numComps, 0)
    checkEquals(res$dist[[2]]$numComps, 0)
}

#test_modelSelection_relabs <- function() {
    # doesn't work
    #res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="rel.abs", nvar=0.5)
#}

test_modelSelection_fixednum <- function() {
    res <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="fixed.num", PCnum=3, center=TRUE, scale=FALSE, weight=TRUE, plot_common=FALSE, plot_dist=FALSE)
    checkEquals(res$common$commonComps, 2)
    checkEquals(res$dist[[1]]$numComps, 1)
    checkEquals(res$dist[[2]]$numComps, 1)
}
