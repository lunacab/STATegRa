
rm(list = ls());
set.seed(12345)
source('../../R/STATegRa_combiningMappings.R')

mapping1 <- data.frame(gene = c('a', 'a', 'a', 'b', 'b', 'c'),
                       transcript = c('ta1', 'ta2', 'ta3', 'tb1', 'tb2', 'tc1'))
mapping2 <- data.frame(gene = c('a', 'a', 'b', 'd'),
                       methylation = c('ma1', 'ma2', 'mb1', 'md1'))
mapping3 <- data.frame(gene = c('a', 'c'),
                       SNP = c('sa1', 'sc1'))

a <- combiningMappings(mappings = list(expression = mapping1, cpg = mapping2))
b <- combiningMappings(mappings = list(expression = mapping1, cpg = mapping2), 
                       retainAll = TRUE)
print(a)
print(b)

a <- combiningMappings(mappings = list(expression = mapping1, 
                                       cpg = mapping2, 
                                       genome = mapping3))

b <- combiningMappings(mappings = list(expression = mapping1, 
                                       cpg = mapping2, 
                                       genome = mapping3),
                       retainAll = TRUE)
print(a)
print(b)
