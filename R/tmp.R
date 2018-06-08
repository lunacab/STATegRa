
library(limma)

dataMatrix <- matrix(rnorm(1000), 100, 10)
outcome <- c(rep(1,5), rep(0, 5))

design <- model.matrix(~outcome)

fit <- lmFit(object = dataMatrix, design = design)
fit <- eBayes(fit)
names(topTable(fit, coef = 2))

library(survival)

outcome <- Surv(time = runif(10), event = c(rep(1,5), rep(0, 5)))
survModel <- summary(coxph(outcome ~ ., data = dataMatrix))$coefficients
