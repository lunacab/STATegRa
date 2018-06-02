#' @export
#' @title combiningMappings, combining several mappings for use in the omicsNPC function
#' @description
#' This function combines several annotation so that measurements across different datasets 
#' are mapped to the same reference elements (e.g., genes).
#' The annotations should all be either data frame / matrices, named vectors/lists, or bioMap objects.
#' See the examples for further details
#' @usage combiningMappings(mappings, reference = NULL, retainAll = FALSE)
#'
#' @param mappings List of annotations.
#' @param reference If the annotations are data frame, matrices or bioMap objects, 
#' the name of the column containing the reference elements
#' @param retainAll Logical, if set to TRUE combination not covering all available data modalities are retained.
#'
#' @return A data frame encoding the mapping across several dataset
#'
#' @author Vincenzo Lagani
#'
#' @references
#' Nestoras Karathanasis, Ioannis Tsamardinos and Vincenzo Lagani. omicsNPC: applying the Non-Parametric Combination 
#' methodology to the integrative analysis of heterogeneous omics data. Submitted to PlosONE. 
#'
#' @examples
#' #Example 1
#' #Mapping with data frames
#' mRNA <- data.frame(gene = rep(c('G1', 'G2', 'G3'), each = 2), probeset = paste('p', 1:6, sep = ''));
#' methylation <- data.frame(gene = c(rep('G1', 3), rep('G2', 4)),
#'                                  methy = paste('methy', 1:7, sep = ''));
#' miRNA <- data.frame(gene = c(rep('G1', 2), rep('G2', 1), rep('G3', 2)),
#'                                miR = c('miR1', 'miR2', 'miR1', 'miR1', 'miR2'));
#' mappings <- list(mRNA = mRNA, methylation = methylation, miRNA = miRNA);
#' combiningMappings(mappings = mappings, retainAll = TRUE)
#' 
#' #Example 2
#' #Mapping with character vectors
#' mRNA <- rep(c('G1', 'G2', 'G3'), each = 2);
#' names(mRNA) = paste('p', 1:6, sep = '');
#' methylation <- c(rep('G1', 3), rep('G2', 4));
#' names(methylation) = paste('methy', 1:7, sep = '');
#' miRNA <- c(rep('G1', 2), rep('G2', 1), rep('G3', 2));
#' names(miRNA) = c('miR1', 'miR2', 'miR1', 'miR1', 'miR2');
#' mappings <- list(mRNA = mRNA, methylation = methylation, miRNA = miRNA);
#' combiningMappings(mappings = mappings, retainAll = TRUE)
#' 

if(!require(data.table)){
  install.packages('data.table')
  library(data.table)
}

combiningMappings <- function(mappings, reference = NULL, retainAll = FALSE){
  
  #validity checking
  if(!(is.list(mappings) && is.logical(retainAll))){
    stop('Mappings must be a list, retainAll must be a logical')
  }
  
  #let's check that the list mappings contains either 
  # - data frame / matrices or 
  # - lists / character vectors or
  # - bioMaps objects
  classes <- sapply(mappings, class);
  if(!all(classes %in% c('data.frame', 'matrix'))){
    if(!all(classes %in% c('list', 'character'))){
      if(!all(classes == 'bioMap') )
      stop('Individual mappings must all be either (a) data frames / matrices or (b) named lists / character vectors or (c) bioMap objects') 
    }
  }
  
  #if no names for the mappings, create fake ones
  if(is.null(names(mappings))){
    names(mappings) <- paste('dataset', 1:length(mappings), sep = "_")
  }
  
  #if the mappings are lists, then let's have them as data frames
  if(all(classes %in% c('list', 'character'))){
    
    #converting each individual mapping
    for(i in 1:length(mappings)){
      mappings[[i]] <- data.frame(names(mappings[[i]]), as.character(mappings[[i]]));
      names(mappings[[i]]) <- c(paste('Measurements_dataset_', i, sep = ''), 'Reference_elements');
    }
    
  }
  
  #if the mappings are bioMap objects, then let's have them as data frames
  if(all(classes == 'bioMap')){
    
    #converting each individual mapping
    for(i in 1:length(mappings)){
      map <- getMap(mappings[[i]])[, c(1,2)];
      names(map)[1] <- getMetadata(mappings[[i]])$type_v1;
      names(map)[2] <- getMetadata(mappings[[i]])$type_v2;
      mappings[[i]] <- map;
    }
    
  }
  
  #let's check that each mapping has two columns
  numColumns <- sapply(mappings, function(x){dim(x)[2]});
  if(any(numColumns != 2)){
    stop('Individual mappings must have two columns')
  }
  
  #ensuring all columns are characters
  for(i in 1:length(mappings)){
    mappings[[i]][, 1] <- as.character(mappings[[i]][, 1])
    mappings[[i]][, 2] <- as.character(mappings[[i]][, 2])
  }
  #checking on reference
  if (is.null(reference)){
    
    #mappings must have exactly one column in common
    reference <- colnames(mappings[[1]]);
    for(i in 2:length(mappings)){
      reference <- intersect(reference, colnames(mappings[[i]]));
    }
    if(length(reference) != 1){
      stop('If no reference column is indicated, individual mappings must have their columns named and exactly one column must be in common across them')
    } 
    
  }else{
    
    #only one reference allowed
    if(length(reference) != 1){
      stop('Please specify only one reference column')
    } 
    
    #all mappings must have this column
    for(i in 1:length(mappings)){
      if(!any(colnames(mappings[[i]]) == reference)){
        stop(paste('All mappings must contain the reference column', reference));
      }
    }
  }
  
  #sanitazing the mapping names
  for(i in 1:length(mappings)){
    toChange <- which(colnames(mappings[[i]]) != reference);
    colnames(mappings[[i]])[toChange] <- paste(colnames(mappings[[i]])[toChange], 
                                               names(mappings)[i], sep = '_');
  }
  
  #eliminating rows that are all NA
  for(i in 1:length(mappings)){
    toKeep <- apply(mappings[[i]], 1, function(x){!all(is.na(x))})
    mappings[[i]] <- mappings[[i]][toKeep, ];
  }
  
  #creating the mapping!
  previousMapping <- as.data.frame(mappings[[1]], stringAsFactors = FALSE);
  for(i in 2:length(mappings)){
    
    #merging the two mappings
    currentMapping <- as.data.frame(mappings[[i]], stringAsFactors = FALSE);
    mapping <- merge(previousMapping, currentMapping, by = reference, all = retainAll);
    
    #if we shall retain all elements...
    if(retainAll){
    
      #elements of previousMapping that are not in mapping
      augmentedMapping <- previousMapping;
      augmentedMapping[ , setdiff(colnames(mapping), colnames(augmentedMapping))] <- NA;
      augmentedMapping <- augmentedMapping[ , colnames(mapping)]
      toAdd <- setdiffMatrix(augmentedMapping, mapping)
      
      #if needed, adding the rows toAdd
      if(length(toAdd) > 0){
        
        #rows to add
        rowsToAdd <- augmentedMapping[toAdd, ]
        
        #adding
        mapping <- rbind(mapping, rowsToAdd);
        
      }
      
      #elements of currentMapping that are not in mapping
      augmentedMapping <- currentMapping;
      augmentedMapping[ , setdiff(colnames(mapping), colnames(currentMapping))] <- NA;
      augmentedMapping <- augmentedMapping[ , colnames(mapping)]
      toAdd <- setdiffMatrix(augmentedMapping, mapping)
      
      #if needed, adding the rows toAdd
      if(length(toAdd) > 0){
        
        #rows to add
        rowsToAdd <- augmentedMapping[toAdd, ]
        
        #adding
        mapping <- rbind(mapping, rowsToAdd);
        
      }
        
    }
    
    #looping...
    previousMapping <- mapping;
    
  }
  
  #removing the reference from the mapping and 
  referenceElements <- mapping[[reference]];
  mapping[[reference]] <- NULL;
  
  # #eliminating duplicated
  # toKeep <- which(!duplicated(mapping));
  # mapping <- mapping[toKeep, ];
  # referenceElements <- referenceElements[toKeep];
  
  #making reference as the last column
  mapping[[reference]] <- referenceElements; #this will put the reference as last column
  
  #ensuring uniquiness
  newRownames <- apply(mapping, 1, function(x){paste(x, collapse = ':')});
  toKeep <- !duplicated(newRownames);
  mapping <- mapping[toKeep, ];
  rownames(mapping) <- NULL;
  
  #naming the columns of mapping
  colnames(mapping) <- c(names(mappings), reference);
  
  #let's have all columns of mappings as characters
  for(i in 1:(dim(mapping)[2])){
    mapping[[i]] <- as.character(mapping[[i]]);
  }
  
  #removing rows with only NA values
  isAllNAs <- which(apply(mapping, 1, function(x){all(is.na(x))}))
  if(length(isAllNAs) > 0){
    mapping <- mapping[-isAllNAs, ]
  }
  
  #removing rows with only one non-NA values
  numNAs <- apply(mapping, 1, function(x){sum(is.na(x))})
  toRemove <- which(numNAs == (dim(mapping)[2] - 2))
  if(length(toRemove) > 0){
    mapping <- mapping[-toRemove, ]
  }
  
  #ordering according to reference
  mapping <- mapping[order(mapping[[reference]]), ]
  
  #return
  return(mapping);
  
}

setdiffMatrix <- function(m1, m2){
  
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=colnames(m2))
  for(i in 1:dim(DT1)[2]){
    DT1[[i]] <- as.character(DT1[[i]])
  }
  setnames(DT2, c(names(DT2)[1:(dim(DT2)[2]-1)], "found"))
  intersection <- DT2[DT1, ];
  inM1notInM2 <- which(is.na(intersection$found))
  
}

