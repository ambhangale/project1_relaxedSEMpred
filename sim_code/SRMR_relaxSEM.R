## Aditi M. Bhangale
## Last updated: 31 March 2026

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
      
      res2 <- lapply(1:p, function(i) (M[i]-mu.hat[i])^2/S[i,i]^2)
      num2 <- do.call("sum", res2) # numerator 2 (mean structure)
      
      num <- num1 + num2
    } else {
      num <- num1
    }
  } else {
    ## relevant for SRMR for yy part of the matrix, since if n_y = 1L, 
    ## S and M will simply contain the variance and mean values, respectively
    p <- length(S)
    
    num1 <- (S-sigma.hat)^2/(S*S) 
    
    if (include.mean) {
      if (!(length(M) == 1L & length(mu.hat) == 1L)) stop("`M` and `mu.hat` must be a numeric value if `S` and `sigma.hat` are 1x1 matrices")
      num2 <- (M-mu.hat)^2/(S*S)
      
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
## however, it will not work for non-symmetric matrices (and hence, does not work for the xy component)

# all the elements in xy are unique, so i have to compute the Frobenius norm and standardise it
SRMR.xy <- function(S, Sxy, sigmaxy.hat, M = NULL, mu.hat = NULL, xnames, ynames, include.mean = F) {
  # if(is.matrix(Sxy) & is.matrix(sigmaxy.hat) & is.matrix(S)) {
  if(length(ynames) > 1L) {
    # squared fitted residuals
    res1 <- lapply(ynames, function(y) 
      lapply(xnames, function(x) (Sxy[x,y]-sigmaxy.hat[x,y])^2/(S[y,y]*S[x,x])))
    num1 <- do.call("sum", do.call("c", res1)) # numerator
    
    if (include.mean) {
      res2 <- lapply(1:length(M), function(i) (M[i]-mu.hat[i])^2/S[i,i]^2)
      num2 <- do.call("sum", res2) # numerator 2 (mean structure)
      
      num <- num1 + num2
      den <- nrow(Sxy)*ncol(Sxy) + length(M) # number of unique elements
    } else {
      num <- num1 
      den <- nrow(Sxy) * ncol(Sxy) # number of unique elements
    }
  } else {
    ## relevant when n_y = 1L and Sxy is a vector
    var.y <- S[ynames, ynames]
    
    res1 <- lapply(xnames, function(x) (Sxy[x]-sigmaxy.hat[x])^2/(var.y*S[x,x]))
    num1 <- do.call("sum", res1) # numerator
    
    if (include.mean) {
      res2 <- lapply(1:length(M), function(i) (M[i]-mu.hat[i])^2/S[i,i]^2)
      num2 <- do.call("sum", res2) # numerator 2 (mean structure)
      
      num <- num1 + num2
      den <- length(Sxy) + length(M)
    } else {
      num <- num1
      den <- length(Sxy)
    }
  }
  
  SRMR <- sqrt(num/den) 
  
  return(SRMR)
}

