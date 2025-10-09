## Aditi M. Bhangale
## Last updated: 8 October 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to calculate the standardised Frobenius norm (SRMR)

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

# SRMR formula
SRMR <- function(S, sigma.hat) {
  if(is.matrix(S) & is.matrix(sigma.hat)) {
    if (!all(dim(S) == dim(sigma.hat)))  stop("`S` and `sigma.hat` must have the same dimensions")
    if(!(nrow(S) == ncol(S) | nrow(sigma.hat) == ncol(sigma.hat))) stop("please specify square matrices only")
    
    p <- nrow(S) # when S is not a matrix
    den <- p*(p+1)*0.5 # denominator
    
    # squared fitted residuals
    res <- lapply(1:p, function(i) lapply(1:i, function(j) (S[i,j]-sigma.hat[i,j])^2/(S[i,i]*S[j,j])))
    num <- do.call("sum", do.call("c", res)) # numerator
    # above is same as num <- mean(do.call("c", do.call("c", res)))
    
    SRMR <- sqrt(num/den)
  } else {
    SRMR <- sqrt((S-sigma.hat)^2/(S*S))
  }
  return(SRMR)
}

## output of this function for the full covariance matrix will match `lavResiduals(fit)$summary["srmr",]`
## the function works for symmetric matrices
## however, it will not work for non-symmetric matrices (and hence, does not work for the yx component)

# all the elements in yx are unique, so i have to compute the Frobenius norm and standardise it
SRMR.yx <- function(S, Syx, sigmayx.hat, ynames, xnames) {
  if(is.matrix(Syx) & is.matrix(sigmayx.hat) & is.matrix(S)) {
    
    den <- nrow(Syx) * ncol(Syx) # number of unique elements
    
    # squared fitted residuals
    num <- 0
    for (y in ynames) {
      for (x in xnames) {
        num <- num + (Syx[y, x] - sigmayx.hat[y, x])^2 / (S[y, y] * S[x, x])
        print(num)
      }
    }
  } else {
    den <- length(Syx)
    
    var.y <- S[ynames, ynames]
    
    res <- lapply(xnames, function(x) (Syx[x]-sigmayx.hat[x])^2/(var.y*S[x,x]))
    num <- do.call("sum", res) # numerator
  }
  
  SRMR <- sqrt(num/den)
  
  return(SRMR)
}

