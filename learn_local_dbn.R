## Goal: Learn local Dynamic Bayesian Network (DBN)
##
LearnLocalDbn <- function(local.dbn.input.data, 
                               scoring.func) {
  
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
  curr.subset <- c()
  curr.subset.str <- as.list(FALSE, num.sl.src.nodes)
  
  ## Current model string
  curr.model.str <- paste('[', tgt.node.name, ']', sep = '')
  ##> [v1_t2] 
  
  ## Current network model
  curr.net <- bnlearn::model2network(curr.model.str)
  
  if (scoring.func == 'BIC') {
    curr.score <- bnlearn::score(curr.net, data, type = 'bic')
  }
  
  ########################################################################
  ## End: Find out the subset of source nodes with the best score
  ########################################################################
}

#######################################################################################
next.candidate.src.node.string <- calc.next.candidate.parent.string(curr.candidate.src.node.string, num.candidate.src.nodes)
next.candidate.src.node.set <- candidate.src.node.names[next.candidate.src.node.string]

calc.next.candidate.parent.string <- function(curr.candidate.src.node.string, num.candidate.src.nodes) {
  
  next.candidate.src.node.string <- vector(n)
  
  least.sig.bit <- curr.candidate.src.node.string[num.candidate.src.nodes]
  
  carry.bit <- FALSE
  
  if (least.sig.bit) {
    next.candidate.src.node.string[num.candidate.src.nodes] <- FALSE
    carry.bit <- TRUE
  } else {
    next.candidate.src.node.string[num.candidate.src.nodes] <- TRUE
  }
  
  for (loop.idx in (num.candidate.src.nodes - 1):1) {
    curr.bit <- curr.candidate.src.node.string[loop.idx]
    
    if (curr.bit & carry.bit) {
      next.candidate.src.node.string[loop.idx] <- FALSE
      carry.bit <- TRUE
    } else if (!curr.bit & !carry.bit) {
      next.candidate.src.node.string[loop.idx] <- FALSE
      carry.bit <- FALSE
    } else {
      next.candidate.src.node.string[loop.idx] <- TRUE
      carry.bit <- FALSE
    }
    
  }
  rm(loop.idx)
}