## Goal: Learn local Dynamic Bayesian Network (DBN)
##
LearnLocalDbn <- function(local.dbn.input.data, 
                               scoring.func) {
  
  ## Number of nodes in the local DBN
  num.nodes.local.dbn <- ncol(local.dbn.input.data)
  
  ## E.g., [1] "v1_t2" 
  tgt.node.name <- colnames(local.dbn.input.data)[ncol(local.dbn.input.data)]
  
  ## Number of shortlisted source nodes.
  ## The last col is for the target node.
  num.sl.src.nodes <- (ncol(local.dbn.input.data) - 1) 
  
  ## E.g., [1] "v1_t1" "v2_t1" "v3_t1"
  sl.src.node.names <- colnames(local.dbn.input.data)[1:num.sl.src.nodes]
  
  ########################################################################
  ## Begin: Find out the subset of source nodes with the best score
  ########################################################################
  
  ## Begin with the empty subset
  # curr.subset <- c()
  curr.subset.str <- as.list(rep(FALSE, num.sl.src.nodes))
  
    ########################################################################
    ## Begin: Compose model string for the source nodes
    ########################################################################
    src.model.str <- c()
    
    for (node.name in sl.src.node.names) {
      node.to.add <- paste('[', 
                           node.name, 
                           ']', 
                           sep = '')
      
      src.model.str <- paste(src.model.str, 
                              node.to.add, 
                              sep = '')
    }
    rm(node.name)
    ##> 'src.model.str' = '[v1_t1][v2_t1][v3_t1]'
    ########################################################################
    ## End: Compose model string for the source nodes
    ########################################################################
  
  ## Model string for the target node
  tgt.model.str <- paste('[', 
                         tgt.node.name, 
                         ']', 
                         sep = '')
  ##> 'tgt.model.str' = '[v1_t2]'
  
  ## Current model string
  curr.model.str <- paste(src.model.str, 
                          tgt.model.str, 
                          sep = '')
  ##> 'curr.model.str' = '[v1_t1][v2_t1][v3_t1][v1_t2]'
  
  ## Current network model
  curr.net <- bnlearn::model2network(curr.model.str)
  
  ## Initialize the current score
  curr.score <- NULL
  
    ########################################################################  
    ## Begin: Calc score of the empty subset
    ########################################################################
    if (scoring.func == 'BIC') {
      curr.score <- bnlearn::score(curr.net, 
                                   local.dbn.input.data, 
                                   type = 'bic')
    } else {
      stop('Only supports BIC scoring function till now')
    }
    ########################################################################
    ## End: Calc score of the empty subset
    ########################################################################
  
  best.score <- curr.score
  best.subset.str <- curr.subset.str
  
    ########################################################################  
    ## Begin: Calc scores of the non-empty subsets
    ########################################################################
  
    num.subsets <- (2^num.sl.src.nodes)
  
    for (subset.idx in 2:num.subsets) {
      
      ## 'find.next.subset.str()' is defined in this script
      curr.subset.str <- find.next.subset.str(curr.subset.str)
      ##> 'curr.subset.str' = list(TRUE, FALSE, TRUE)
      
      ## Model string for the target node
      tgt.model.str <- paste('[', 
                             tgt.node.name, 
                             '|', 
                             sep = '')
      ##> 'tgt.model.str' = '[v1_t2|'
      
      is.first <- TRUE
      
      for (node.idx in length(curr.subset.str)) {
        if (curr.subset.str[[node.idx]]) {
          
          if (is.first) {
            is.first <- FALSE
            
          } else {
            tgt.model.str <- paste(tgt.model.str, 
                                   ':', 
                                   sep = '')
          }
          
          tgt.model.str <- paste(tgt.model.str, 
                                 sl.src.node.names[node.idx], 
                                 sep = '')
        }
      }
      rm(node.idx)
      
      tgt.model.str <- paste(tgt.model.str, 
                             ']',
                             sep = '')
      ##> tgt.model.str = '[v1_t2|v1_t1:v1_t3]'
      
      ## Current model string
      curr.model.str <- paste(src.model.str, 
                              tgt.model.str, 
                              sep = '')
      ##> 'curr.model.str' = '[v1_t1][v2_t1][v3_t1][v1_t2|V1_t1:v1_t3]'
      
      ## Current network model
      curr.net <- bnlearn::model2network(curr.model.str)
      
      if (scoring.func == 'BIC') {
        curr.score <- bnlearn::score(curr.net, 
                                     local.dbn.input.data, 
                                     type = 'bic')
      } else {
        stop('Only supports BIC scoring function till now')
      }
      
      if (curr.score > best.score) {
        best.score <- curr.score
        best.subset.str <- curr.subset.str
      }
      
    }
    rm(subset.idx, num.subsets)
    ########################################################################
    ## End: Calc scores of the non-empty subsets
    ########################################################################    
  
  ########################################################################
  ## End: Find out the subset of source nodes with the best score
  ########################################################################
}

#######################################################################################
## Goal: Find next subset string, given the previous subset string.
## Example 1:
## If 'prev.subset.str' = list(FALSE, FALSE)
## then it returns 'list(FALSE, TRUE)'.
##
## Example 2:
## If 'prev.subset.str' = list(FALSE, TRUE)
## then it returns 'list(TRUE, FALSE)'.
##
## So this function basically adds a TRUE to the last element of the
## list 'prev.subset.str'.
#######################################################################################
find.next.subset.str <- function(prev.subset.str) {
  
  num.nodes.local.dbn <- length(prev.subset.str)
  
  ## Initialize the carry bit
  carry.bit <- FALSE
  
  for (node.idx in num.nodes.local.dbn:1) {
    
    ## Cuurent bit
    curr.bit <- prev.subset.str[[node.idx]]
    
    if (!curr.bit & !carry.bit) {
      prev.subset.str[[node.idx]] <- FALSE
      carry.bit <- FALSE
    } else if (curr.bit & carry.bit) {
      prev.subset.str[[node.idx]] <- FALSE
      carry.bit <- TRUE
    } else {
      prev.subset.str[[node.idx]] <- TRUE
      carry.bit <- FALSE
    }
  }
  rm(loop.idx)
  
  ## It is an in-place replacement i.e.
  ## the content of the previous subset
  ## string is modified to represent
  ## the next subset string
  return(prev.subset.str)
}