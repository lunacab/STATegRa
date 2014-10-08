data(STATegRa_S1)
data(STATegRa_S3)

test_createOmicsExpressionSet <- function() {
    
    
    createOmicsExpressionSet(Data = Block1)
    #.PCA has no rownames, so doesn't work without pData
    checkException(createOmicsExpressionSet(Data = Block1.PCA))
    
    # pDataDescr should be auto-guessed
    createOmicsExpressionSet(Data = Block1.PCA, pData = ed.PCA, pDataDescr = c("classname"))
    checkEquals(createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA),
                createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA, pDataDescr=c("classname")))
    
    # unmatching or incomplete pData to be rejected
    checkException(createOmicsExpressionSet(Data=Block1, pData=ed[-1,,drop=F]))
    
    # lacking test data for feaData
}
