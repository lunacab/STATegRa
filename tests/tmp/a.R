statistics <- matrix(statsNPC[, 1, ], nrow = dim(statsNPC)[1], ncol = dim(statsNPC)[3])
rownames(statistics) <- dimnames(statsNPC)[[1]]
colnames(statistics) <- dimnames(statsNPC)[[3]]
dim(statistics)

q <- quantile(perm_stats,0.9)
pi0 <- (2*sum(init_stats <= q)) / length(init_stats)

# dataInput = list(dataset1 = dataset1, dataset2 = dataset2)
# dataTypes = rep('continuous', length(dataInput))
# outcomeName = NULL
# combMethods = c("Fisher", "Liptak", "Tippett")
# allCombinations = FALSE
# verbose = FALSE
# functionGeneratingIndex = NULL
# dataWeights = rep(1, length(dataInput))/length(dataInput)
# returnPermPvalues = FALSE

# dataInput = list(dataset1 = dataset1, dataset2 = dataset2, dataset3 = dataset3)
# dataTypes = rep('continuous', length(dataInput))
# outcomeName = NULL
# combMethods = c("Fisher", "Liptak", "Tippett")
# allCombinations = TRUE
# verbose = FALSE
# functionGeneratingIndex = NULL
# dataWeights = rep(1, length(dataInput))/length(dataInput)
# returnPermPvalues = FALSE
