data(STATegRa_S3)

test_selectCommonComps <- function() {
    checkEquals(selectCommonComps(Block1.PCA, Block2.PCA, Rmax=3)$common, 2)
}
