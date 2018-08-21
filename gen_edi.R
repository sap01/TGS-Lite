## Goal: Generate synthetic networks and synthetic datasets using 
## R package 'EDISON' version 1.1.1
## 
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params:
## See the input params of function 'generateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##
## Output params:
## See the output params of function 'generateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##################################################################################################
GenEdiNet <- function(lambda_2=0.45, q=10, min_phase_length=1, k_bar=5, l=10, 
                       lambda_3=2, spacing=1, gauss_weights=TRUE, same=FALSE, 
                       change_method='sequential',
                       fixed=FALSE, cps=NULL) {
  
  ##------------------------------------------------------------
  ## Begin: Load the required packages
  ##------------------------------------------------------------
  library(EDISON)
  ##------------------------------------------------------------
  ## End: Load the required packa                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ges
  ##------------------------------------------------------------
  
  ## Generate synthetic networks with the given parameters
  edi.net <- EDISON::generateNetwork(lambda_2 = lambda_2, 
                                     q = q, 
                                     min_phase_length = min_phase_length,
                                     k_bar = k_bar, 
                                     l = l, 
                                     lambda_3 = lambda_3, 
                                     spacing = spacing, 
                                     gauss_weights = gauss_weights,
                                     same = same, 
                                     change_method = change_method, 
                                     fixed = fixed, 
                                     cps = cps)
  
  return(edi.net)
}
##################################################################################################

##################################################################################################
## Goal: Generate synthetic time-series gene expression data time-varying gene regulatory networks
## using R package 'EDISON' version 1.1.1
##
## Input params: 
## 'edi.net': A set of time-varying gene regulatory networks.
## See input param 'net' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'noise': See input param 'noise' of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
## 'num.time.series': Number of time series to be generated.
##
## Output params:
## See the output params of function 'simulateNetwork()' 
## in R package 'EDISON' version 1.1.1.
##################################################################################################
GenEdiData <- function(edi.net, noise, num.time.series) {
  
  ##------------------------------------------------------------
  ## Begin: Load the required packages
  ##------------------------------------------------------------
  library(EDISON)
  ##------------------------------------------------------------
  ## End: Load the required packa                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ges
  ##------------------------------------------------------------
  
  ## TODO(sap): create a .tsv file comaptible with TGS.R
  ## S = 46
  
  ## Simulate data using given networks
  edi.data <- EDISON::simulateNetwork(noise = 0.05, net = edi.net)  
  
}
##################################################################################################
  
