## Aditi M. Bhangale
## Last updated: 9 October 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to calculate the standardised Frobenius norm (SRMR)

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

# SRMR formula
SRMR <- function(S, sigma.hat, M = NULL, mu.hat = NULL, include.mean = F) {
  if(is.matrix(S) && is.matrix(sigma.hat)) {
    if (!all(dim(S) == dim(sigma.hat)))  stop("`S` and `sigma.hat` must have the same dimensions")
    if(!(nrow(S) == ncol(S) | nrow(sigma.hat) == ncol(sigma.hat))) stop("please specify square matrices only")
    
    p <- nrow(S) # number of observed variables
    
    # squared fitted residuals
    res1 <- lapply(1:p, function(i) lapply(1:i, function(j) (S[i,j]-sigma.hat[i,j])^2/(S[i,i]*S[j,j])))
    num1 <- do.call("sum", do.call("c", res1)) # numerator 1 (covariance structure)
    # above is same as num1 <- mean(do.call("c", do.call("c", res1)))
    
    if (include.mean) {
      if (!(length(M) > 1L & length(mu.hat) > 1L)) stop("`M` and `mu.hat` must be vectors if `S` and `sigma.hat` are square matrices")
      if (!(length(M) == length(mu.hat)))  stop("`M` and `mu.hat` must have the same length")
      
      res2 <- lapply(1:p, function(i) (M[i]-mu.hat[i])^2/(S[i,i]*S[i,i]))
      num2 <- do.call("sum", res2) # numerator 2 (mean structure)
      
      num <- num1 + num2
    } else {
      num <- num1
    }
  } else {
    ## relevant for SRMR for yy part of the matrix, since if n_y = 1L, 
    ## S and M will simply contain the variance and mean values, respectively
    p <- length(S)
    
    num1 <- sqrt((S-sigma.hat)^2/(S*S)) 
    
    if (include.mean) {
      if (!(length(M) == 1L & length(mu.hat) == 1L)) stop("`M` and `mu.hat` must be a numeric value if `S` and `sigma.hat` are 1x1 matrices")
      num2 <- sqrt((M-mu.hat)^2/(S*S)) 
      
      num <- num1 + num2
    } else {
      num <- num1
    }
  }
  
  den <- ifelse(include.mean, p*(p+3)*0.5, p*(p+1)*0.5) # denominator
  
  SRMR <- sqrt(num/den)
  
  return(SRMR)
}

## output of this function for the full covariance matrix will match `lavResiduals(fit)$summary["srmr","cov"]`
## when `include.mean = F` and `lavResiduals(fit)$summary["srmr","total"]` when `include.mean = T`. 
## the function works for symmetric matrices
## however, it will not work for non-symmetric matrices (and hence, does not work for the yx component)

# all the elements in yx are unique, so i have to compute the Frobenius norm and standardise it
SRMR.yx <- function(S, Syx, sigmayx.hat, M = NULL, mu.hat = NULL, ynames, xnames, include.mean = F) { # TODO add include.mean argument
  if(is.matrix(Syx) & is.matrix(sigmayx.hat) & is.matrix(S)) {
    
    den <- nrow(Syx) * ncol(Syx) # number of unique elements
    
    # squared fitted residuals
    res1 <- lapply(ynames, function(y) 
      lapply(xnames, function(x) (Syx[y,x]-sigmayx.hat[y,x])^2/(S[y,y]*S[x,x])))
    num1 <- do.call("sum", do.call("c", res1)) # numerator
    
    if (include.mean) {
      #FIXME fix code below this
      res2 <- lapply(1:p, function(i) lapply(1:i, function(j) (M[i]-mu.hat[i])^2/(S[i,i]*S[i,i])))
      num2 <- do.call("sum", do.call("c", res2)) # numerator 2 (mean structure)
    } else {
     num <- num1 
    }
  } else {
    ## relevant when n_y = 1L
    
    den <- length(Syx)
    
    var.y <- S[ynames, ynames]
    
    res <- lapply(xnames, function(x) (Syx[x]-sigmayx.hat[x])^2/(var.y*S[x,x]))
    num <- do.call("sum", res) # numerator
  }
  
  SRMR <- sqrt(num/den) 
  
  return(SRMR)
}

