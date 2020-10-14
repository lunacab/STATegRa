
#this function computes association for continuous data
ttCoxphContinuous <- function(dataMatrix, design, outcomeName, useVoom = FALSE, returnPValues = FALSE, ...){
  
  #requirements
  require(limma);
  require(survival);
  require(edgeR)
  
  # ensuring that there is a column named "age"
  if(!('age' %in% colnames(design))){
    stop('age is not among the covariates')
  }
  
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
  
  # ensuring we have surv outcome
  if(!is.Surv(design[[outcomeId]])){
    stop('The outcome is supposed to be of type "Surv" ')
  }
  
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
    summary(coxph(survOutcome ~ . + tt(age), data = designX, 
                  tt = function(y, t, ...){y * log(t + 20)}))$coefficients['x', ]
    
  }
  
  #retrieving the statistics
  results <- t(apply(dataMatrix, 1, modelFunction, design, outcomeId));
  if(returnPValues){
    statistics <- results[, 'Pr(>|z|)']
  }else{
    statistics <- results[, 'z']
  }
  
  # adding the adjusted p-values for consistency with limma
  results <- cbind(results, p.adjust(results[, 'Pr(>|z|)'], method = 'fdr'))
  colnames(results)[ncol(results)] <- 'adj.P.Val'
  #returning statistics (log transformed p-values)
  names(statistics) <- rownames(dataMatrix);
  rownames(results) <- rownames(dataMatrix);
  return(list(statistics = statistics, fullMatrix = results));  
  
}

#this function computes association for count data
ttCoxphCount <- function(dataMatrix, design, outcomeName, ...){
  
  ttCoxphContinuous(dataMatrix, design, outcomeName, useVoom = TRUE, ...)
  
}
