
#function for generating random, iid permutation indices
generate_iid_data_index <- function(design){
  numSamples <- dim(design)[1];
  index <- sample(x = 1:numSamples, size = numSamples, replace = FALSE)
  return(index)
}

#this function computes association for continuous data
computeAssocContinuousData <- function(dataMatrix, design, outcomeName, useVoom = FALSE, returnPValues = FALSE, ...){
  
  #save(list = ls(), file= 'safety.RData')
  #dataMatrix <- dataMatrices[[1]]
  #design <- designs[[1]]
  
  #requirements
  require(limma);
  require(survival);
  require(edgeR)
  
  #ensuring the outcome is the last column of the design
  if(!is.null(outcomeName)){
    
    #check
    if(!any(colnames(design) == outcomeName)){
      stop(paste('Ensure', outcomeName, 'is present in all design matrices'));
    }
    
    #moving the outcome to last position
    outcomeTemp <- design[[outcomeName]];
    design[[outcomeName]] <- NULL;
    design <- cbind(design, outcomeTemp);
    names(design)[ncol(design)] <- outcomeName;
    
  }
  
  #we assume that the outcome is the last column
  outcomeId <- ncol(design);
  
  #what type of outcome?
  multiGroups <- FALSE;
  surv <- FALSE;
  if((is.factor(design[[outcomeId]]) & length(levels(design[[outcomeId]])) >= 3) ||
     (is.character(design[[outcomeId]]) & length(unique(design[[outcomeId]])) >= 3)){
    multiGroups <- TRUE;
  }
  if(is.Surv(design[[outcomeId]])){
    surv <- TRUE;
  }
  
  #linear models or survival models
  if(!surv){
    
    #creating the model matrix 
    numGroups <- length(unique(design[[outcomeId]]));
    modelFormula <- as.formula(paste('~', paste(colnames(design), collapse = '+')))
    #save(modelFormula, design, file = 'safety.Rdata')
    design <- model.matrix(modelFormula, design)
    numCols <- ncol(design)
    
    #if there are n groups, the last n-1 columns represent them
    numColGroups <- (ncol(design) - (numGroups - 1) + 1):ncol(design)
    
    #use voom?
    if(useVoom){
      nf = calcNormFactors(dataMatrix, method = "TMM")
      tmpGenes = rownames(dataMatrix)
      dataMatrix = voom(dataMatrix, design = design,
                        lib.size = colSums(dataMatrix) * nf)
      dataMatrix$genes = tmpGenes
      
    }
    
    #fitting the models
    fit <- lmFit(dataMatrix, design = design)
    fit2 <- eBayes(fit)
    
    #if not multiple groups
    if(!multiGroups){
      results <- topTable(fit2, numCols, number = Inf, sort.by = 'none')
      if(returnPValues){
        statistics <- results$P.Value
      }else{
        statistics <- results$t;
      }
    }else{
      results <- topTable(fit2, numColGroups, number = Inf, sort.by = 'none')
      if(returnPValues){
        statistics <- results$P.Value
      }else{
        statistics <- results$F;
      }
    }
    
    #returning statistics
    names(statistics) <- rownames(results);
    return(list(statistics = statistics, results = results));  
    
  }else{
    
    #use voom?
    if(useVoom){
      nf = calcNormFactors(dataMatrix, method = "TMM")
      tmpGenes = rownames(dataMatrix)
      if(outcomeId == 1){
        dataMatrix = voom(dataMatrix, lib.size = colSums(dataMatrix) * nf)
      }else{
        tmpDesign <- design[, -outcomeId];
        dataMatrix = voom(dataMatrix, design = tmpDesign,
                          lib.size = colSums(dataMatrix) * nf)
      }
      dataMatrix <- dataMatrix$E;
      rownames(dataMatrix) = tmpGenes;
    }
    
    #model function
    modelFunction <- function(x, design, outcomeId){
      
      #removing the survival outcome
      survOutcome <- design[[outcomeId]];
      design[[outcomeId]] <- NULL;
      
      #including x in the model
      if(outcomeId == 1){
        designX <- data.frame(x = x);; 
      }else{
        designX <- cbind(design, x); 
      }
      
      #return
      summary(coxph(survOutcome ~ ., data = designX))$coefficients['x', ]
      
    }
    
    #retrieving the statistics
    results <-t(apply(dataMatrix, 1, modelFunction, design, outcomeId));
    if(returnPValues){
      statistics <- results[, 'Pr(>|z|)']
    }else{
      statistics <- results[, 'z']
    }
    
    #returning statistics (log transformed p-values)
    names(statistics) <- rownames(dataMatrix);
    rownames(results) <- rownames(dataMatrix);
    return(list(statistics = statistics, results = results));  
    
  }
  
}

#this function computes association for count data
computeAssocCountData <- function(dataMatrix, design, outcomeName, ...){
  
  computeAssocContinuousData(dataMatrix, design, outcomeName, useVoom = TRUE, ...)
  
}
