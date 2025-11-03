## Aditi M. Bhangale
## Last updated: 3 November 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### data generation file -- create data from scratch -- following Yves' approach

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

library(mvtnorm)
# library(lavaan)
library(portableParallelSeeds) # remotes::install_github("wjakethompson/portableParallelSeeds")

# create seed object for data generation and partition generation----
allSeeds <- seedCreator(nReps = 5e3, streamsPerRep = 2, seed = 10824)

# possible to run sampIDs 1:5000. adjust nReps if we want to go higher.
# 2 streams. stream 1 to be used for data generation and stream 2 to be used for partitioning
#----

# n_x = 24; n_eta_x = 3; n_y = 9; n_eta_y = 1
# beta = 0.3; lambda = 0.7; psi.cov = 0.2; r = 0.3
# miss.strength = "strong"

# misspecification types----
xxrescov <- function(n_x, n_eta_x, THETA, THETA.star, miss.strength) {
  
  nRC <- round(n_x*(n_x-1)*0.5/6) # number of residual covariances to introduce
  uniq.el <- combn(1:n_x, 2, simplify = F) # unique elements in THETA_xx
  
  # residual covariances to introduce
  RC <- sample(uniq.el, nRC)
  
  for (rc in 1:length(RC)) {
    THETA[paste0("x",RC[[rc]][1]), paste0("x",RC[[rc]][2])] <- 
      THETA[paste0("x",RC[[rc]][2]), paste0("x",RC[[rc]][1])] <- ifelse(miss.strength == "weak", 
                                              0.3*min(THETA.star), 
                                              0.6*min(THETA.star))
  }
  
  return(THETA)
}

xxcrossload <- function(n_x, n_eta_x, lambda, LAMBDA, miss.strength) {
  if (n_eta_x != 3L) stop("xx misspecification not supported for n_eta_x = 1L")
  if (!(n_x == 12L || n_x == 24L)) stop("xxcrossload misspecificaiton only supported for n_x == 12L | 24L for now")
  
  nCL <- ifelse(n_x == 12L, 1, 2) # number of cross-loadings to introduce per factor
  
  for (xx in 1:n_eta_x) {
    if(xx == 1L) {
      LAMBDA[paste0("x", sample((n_x/n_eta_x + 1):n_x, nCL)), 
             paste0("eta_x", xx)] <- ifelse(miss.strength == "weak", 0.5*lambda, 0.9*lambda)
    } else if (xx == 2L) {
      LAMBDA[paste0("x", sample(c(1:(n_x/n_eta_x),(n_x-n_x/n_eta_x+1):n_x), nCL)), 
             paste0("eta_x", xx)] <- ifelse(miss.strength == "weak", 0.5*lambda, 0.9*lambda)
    } else if (xx == 3L) {
      LAMBDA[paste0("x", sample(1:(n_x-n_x/n_eta_x), nCL)), 
             paste0("eta_x", xx)] <- ifelse(miss.strength == "weak", 0.5*lambda, 0.9*lambda)
    }
  }
  
  return(LAMBDA)
}

xydirect <- function(n_x, n_eta_x, n_y, n_eta_y, LAMBDA, B, PSI, THETA, 
                     beta, miss.strength, lvnames, obsnames) {
  
  if (n_eta_x == 1L) {
    # add new column to LAMBDA - create a new single-indicator factor
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x, 
                         dimnames = list(obsnames,"mis.eta_x1"))
    mis.LAMBDA["x1", "mis.eta_x1"] <- 1L
    mis.LAMBDA.val <- LAMBDA["x1", "eta_x1"] # this lambda value will now be in the B matrix
    LAMBDA["x1", "eta_x1"] <- 0L
    mis.LAMBDA <- cbind(LAMBDA, mis.LAMBDA)
    
    # add residual variance for the new single-indicator factor
    ## factor variance of the new single-indicator factor will be the residual indicator variance
    mis.PSI.val <- THETA["x1", "x1"]
    mis.PSI <- rbind(cbind(PSI, rep(0, length(lvnames))), 
                     c(rep(0, length(lvnames)), mis.PSI.val))
    dimnames(mis.PSI) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    # set residual variance to 0 (due to new single-indicator factor)
    THETA["x1", "x1"] <- 0L
    
    # add regression slope of outcome factor on new single-indicator factor
    ## rows are outcomes, columns are predictors
    mis.B.val <- ifelse(miss.strength == "weak", 
                           0.5*beta,
                           0.9*beta)
    mis.B <- rbind(cbind(B, c(0, mis.B.val)), 
                      c(mis.LAMBDA.val, rep(0, length(lvnames))))
    dimnames(mis.B) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    
  } else if (n_eta_x == 3L) {
    mis.lvnames <- c(lvnames, paste0("mis.eta_x", 
                                     c(1, 
                                       (n_x/n_eta_x + 1), 
                                       (n_x-n_x/n_eta_x+1))))
    
    # add new columns to LAMBDA -- create single-indicator factors
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x,
                         dimnames = list(obsnames, mis.lvnames[grep("mis.", mis.lvnames)]))
    mis.LAMBDA["x1", "mis.eta_x1"] <- 1L
    mis.LAMBDA[paste0("x", (n_x/n_eta_x + 1)), 
               paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 1L
    mis.LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1)), 
               paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 1L
    mis.LAMBDA.vals <- c(LAMBDA["x1", "eta_x1"],
                         LAMBDA[paste0("x", (n_x/n_eta_x + 1)), "eta_x2"],
                         LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1)), "eta_x3"])
    names(mis.LAMBDA.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    LAMBDA["x1", "eta_x1"] <- LAMBDA[paste0("x", (n_x/n_eta_x + 1)), "eta_x2"] <- 
      LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1)), "eta_x3"] <- 0L
    mis.LAMBDA <- cbind(LAMBDA, mis.LAMBDA)
    
    # add residual variance for new single-indicator factors
    ## factor variances of the new single-indicator factors will be the residual indicator variances
    mis.PSI.vals <- c(THETA["x1", "x1"],
                        THETA[paste0("x", (n_x/n_eta_x + 1)), paste0("x", (n_x/n_eta_x + 1))],
                        THETA[paste0("x", (n_x-n_x/n_eta_x+1)),  paste0("x", (n_x-n_x/n_eta_x+1))])
    names(mis.PSI.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    
    mis.PSI <- matrix(0, nrow = 2*n_eta_x + n_eta_y, ncol = 2*n_eta_x + n_eta_y,
                      dimnames = list(mis.lvnames, mis.lvnames))
    mis.PSI[lvnames, lvnames] <- PSI
    mis.PSI["mis.eta_x1", "mis.eta_x1"] <- mis.PSI.vals["mis.eta_x1"]
    mis.PSI[paste0("mis.eta_x", (n_x/n_eta_x + 1)), 
            paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 
      mis.PSI.vals[paste0("mis.eta_x", (n_x/n_eta_x + 1))]
    mis.PSI[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1)), 
            paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 
      mis.PSI.vals[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))]
    
    # set residual variances to 0
    THETA["x1", "x1"] <- 
      THETA[paste0("x", (n_x/n_eta_x + 1)), paste0("x", (n_x/n_eta_x + 1))] <- 
      THETA[paste0("x", (n_x-n_x/n_eta_x+1)),  paste0("x", (n_x-n_x/n_eta_x+1))] <- 0L
    
    # add regression slopes of outcome factor on new single-indicator factors
    ## rows are outcomes, columns are predictors
    mis.B.vals <- rep(ifelse(miss.strength == "weak", 0.5*beta, 0.9*beta), n_eta_x)
    names(mis.B.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    
    mis.B <- matrix(0, nrow = 2*n_eta_x + n_eta_y, ncol = 2*n_eta_x + n_eta_y,
                       dimnames = list(mis.lvnames, mis.lvnames))
    mis.B[lvnames, lvnames] <- B
    
    mis.B["mis.eta_x1", "eta_x1"] <- 
      mis.LAMBDA.vals[grep("\\_x1\\b", mis.lvnames[grep("mis.", mis.lvnames)])]
    mis.B[paste0("mis.eta_x", (n_x/n_eta_x + 1)), "eta_x2"] <- 
      mis.LAMBDA.vals[grep(paste0("x", (n_x/n_eta_x + 1)), mis.lvnames[grep("mis.", mis.lvnames)])]
    mis.B[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1)), "eta_x3"] <-
      mis.LAMBDA.vals[grep(paste0("x", (n_x-n_x/n_eta_x+1)), mis.lvnames[grep("mis.", mis.lvnames)])]
    
    ####
    mis.B["eta_y1", "mis.eta_x1"] <- mis.B.vals["mis.eta_x1"]
    mis.B["eta_y1", paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 
      mis.B.vals[paste0("mis.eta_x", (n_x/n_eta_x + 1))]
    mis.B["eta_y1", paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 
      mis.B.vals[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))]
  }
  
  # recompute phi
  Iden <- diag(1, nrow = nrow(mis.B))
  mis.PHI <- solve(Iden - mis.B) %*% mis.PSI %*% t(solve(Iden - mis.B))
  
  return(list(mis.LAMBDA = mis.LAMBDA, mis.THETA = THETA, 
              mis.PSI = mis.PSI, mis.B = mis.B, mis.PHI = mis.PHI))
}

#----

# generate random covariance matrices----
genCovmat <- function(n_x, n_eta_x, n_y,  n_eta_y = 1L, 
                      beta = 0.3, lambda = 0.7, psi.cov = 0.2, r = 0.3,
                      misspecify, miss.part = NULL, miss.strength = NULL) { 
  
  # check if n_x is divisible by n_eta_x, otherwise stop
  if (n_x %% n_eta_x != 0) stop("`n_eta_x` must be divisble by `n_x`")
  
  # check if n_eta_y = 1L, otherwise stop
  if (n_eta_y != 1L) stop("only `n_eta_y == 1L` supported for now")
  
  # check if n_eta_x = 1 or 3, otherwise stop
  if (!(n_eta_x == 1L || n_eta_x == 3L)) stop("only `n_eta_x == 1L or 3L` supported for now")
  
  # check if each factor has 4 or 8 indicators, otherwise warning
  if (!(n_x/n_eta_x == 4L || n_x/n_eta_x == 8L)) warning("all misspecification calculations only implemented for n_x/n_eta_x == 4L | 8L for now")
  
  obsnames <- c(paste0("x", 1:n_x), paste0("y", 1:n_y)) # indicator labels
  lvnames <- c(paste0("eta_x", 1:n_eta_x), paste0("eta_y", 1:n_eta_y)) # factor labels
  
  # B; matrix of regression coefficients (structural relations)
  B <- matrix(0, nrow = n_eta_x + n_eta_y, 
              ncol = n_eta_x + n_eta_y, 
              dimnames = list(lvnames, lvnames)) # rows are outcomes, columns are predictors
  B[lvnames[grep("_y", lvnames)], lvnames[grep("_x", lvnames)]] <- beta
  
  # PSI; factor covariance matrix, not accounting for structural relations
  PSI <- diag(1, nrow = n_eta_x + n_eta_y)
  dimnames(PSI) <- list(lvnames, lvnames)
  PSI[lvnames[grep("_y", lvnames)], lvnames[grep("_y", lvnames)]] <- 1-beta^2
  for (x in lvnames[grep("_x", lvnames)]) {
    for (xx in lvnames[grep("_x", lvnames)]) {
      if (x != xx) {
        PSI[x,xx] <- PSI[xx,x] <- psi.cov # fix factor covariances to some non-zero value
      }
    }
  }
  
  # LAMBDA; factor loading matrix (measurement part)
  LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x + n_eta_y, 
                   dimnames = list(obsnames, lvnames)) # rows are outcomes, columns are predictors
  
  ## for x factors
  if (n_eta_x == 1L) {
    LAMBDA[paste0("x", 1:n_x), "eta_x1"] <- lambda
  } else if (n_eta_x == 3L) {
    for (xx in 1:n_eta_x) { # only conditions wherein all factors have equal number of indicators supported for now
      if(xx == 1L) {
        LAMBDA[paste0("x", 1:(n_x/n_eta_x)), paste0("eta_x", xx)] <- lambda
      } else if (xx == 2L) {
        LAMBDA[paste0("x", (n_x/n_eta_x + 1):(n_x-n_x/n_eta_x)), paste0("eta_x", xx)] <- lambda
      } else if (xx == 3L) {
        LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1):n_x), paste0("eta_x", xx)] <- lambda
      }
    }
  }
  
  ## for y factor
    LAMBDA[paste0("y",1:n_y), "eta_y1"] <- lambda
  
  # PHI; calculate from LAMBDA and B, factor covariances incorporating structural part
  Iden <- diag(1, nrow = n_eta_x + n_eta_y)
  PHI <- solve(Iden - B) %*% PSI %*% t(solve(Iden - B))
  
  # THETA; residual variances of indicators
  THETA.dash <- diag(LAMBDA%*%PHI%*%t(LAMBDA))
  THETA.star <- THETA.dash*(1/r - 1)
  THETA <- diag(THETA.star, nrow = n_x + n_y)
  dimnames(THETA) <- list(obsnames, obsnames)
  
  # if the y-part of the model is a single-indicator factor model
  ## fix factor variance to residual indicator variance, then recompute PHI;
  ## fix residual indicator variance to 0L
  ## fix factor loading for single-indicator factor to 1L
  if (n_y == 1L) {
    PSI["eta_y1", "eta_y1"] <- THETA.star["y1"]
    PHI <- solve(Iden - B) %*% PSI %*% t(solve(Iden - B))
    THETA["y1", "y1"] <- 0L
    LAMBDA["y1", "eta_y1"] <- 1L
  }
  
  # introduce misspecification
  if(misspecify == T) {
    if (is.null(miss.part) | is.null(miss.strength)) {
      stop("specify which part of the model to misspecify and strength of misspecification")
    } else {
      if(!(miss.strength == "weak" || miss.strength == "strong")) stop("specify a valid misspecification strength in `miss.strength`") 
      if(!(miss.part %in% c("xx:rescov", "xx:crossload", "xy:direct", "both:cov", "both:load"))) stop("specify a valid misspecification type in `miss.part`")
      
      ## set seed for randomly selecting parts of model to misspecify
      set.seed(as.numeric(paste0(n_x,n_eta_x,n_y,n_eta_y))) 
      
      if (miss.part == "xx:rescov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, THETA = THETA, 
                          THETA.star = THETA.star, miss.strength = miss.strength)
      } else if (miss.part == "xx:crossload") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = lambda, 
                              LAMBDA = LAMBDA, miss.strength = miss.strength)
      } else if (miss.part == "xy:direct") {
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
      } else if (miss.part == "both:cov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, THETA = THETA, 
                          THETA.star = THETA.star, miss.strength = miss.strength)
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
        
      } else if (miss.part == "both:load") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = lambda, 
                              LAMBDA = LAMBDA, miss.strength = miss.strength)
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
        
      } else stop("specify valid `miss.part`")
    }
  }
  
  # population covariance matrix
  SIGMA.pop <- LAMBDA %*% PHI %*% t(LAMBDA) + THETA
  
  attr(SIGMA.pop, "n_x") <- n_x
  attr(SIGMA.pop, "n_eta_x") <- n_eta_x
  attr(SIGMA.pop, "n_y") <- n_y
  attr(SIGMA.pop, "n_eta_y") <- n_eta_y
  attr(SIGMA.pop, "misspecify") <- misspecify
  attr(SIGMA.pop, "miss.part") <- ifelse(!is.null(miss.part), miss.part, NA)
  attr(SIGMA.pop, "miss.strength") <- ifelse(!is.null(miss.strength), miss.strength, NA)
  
  return(SIGMA.pop)
}

# test function
## correctly specified
# genCovmat(n_x = 4, n_eta_x = 1, n_y = 1, misspecify = F)
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 1, misspecify = F)
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 4, misspecify = F)

## "xx:rescov"
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 1, misspecify = T, 
#           miss.part = "xx:rescov", miss.strength = "strong")
# genCovmat(n_x = 12, n_eta_x = 3, n_y = 4, misspecify = T, 
#           miss.part = "xx:rescov", miss.strength = "weak")

## "xx:crossload"
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 1, misspecify = T, 
#           miss.part = "xx:crossload", miss.strength = "strong")
# genCovmat(n_x = 12, n_eta_x = 3, n_y = 4, misspecify = T, 
#           miss.part = "xx:crossload", miss.strength = "weak")

## "xy:direct"
# genCovmat(n_x = 4, n_eta_x = 1, n_y = 1, misspecify = T, 
#           miss.part = "xy:direct", miss.strength = "strong")
# genCovmat(n_x = 12, n_eta_x = 3, n_y = 1, misspecify = T, 
#           miss.part = "xy:direct", miss.strength = "strong")
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 4, misspecify = T, 
#           miss.part = "xy:direct", miss.strength = "weak")

## "both:cov"
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 1, misspecify = T, 
#           miss.part = "both:cov", miss.strength = "strong")
# genCovmat(n_x = 12, n_eta_x = 3, n_y = 4, misspecify = T, 
#           miss.part = "both:cov", miss.strength = "weak")

# "both:load"
# genCovmat(n_x = 24, n_eta_x = 3, n_y = 1, misspecify = T, 
#           miss.part = "both:load", miss.strength = "strong")
# genCovmat(n_x = 12, n_eta_x = 3, n_y = 4, misspecify = T,
#           miss.part = "both:load", miss.strength = "weak")

#----

# generate data from covariance matrices----
gendat <- function(sampID = NULL, nCal, nPred, covmat, seed = NULL) {
  
  if (!is.null(sampID)) {
    # set seed based on `sampID`
    setSeeds(projSeeds = allSeeds, run = sampID)
    useStream(1) # use the first stream of the `sampID` rep for data generation
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    } else {
      stop("specify `seed` argument if `sampID` is left blank") 
    }
  }
  
  # generate mean-centered data with the covariance matrix generated from `genCovmat()`
  calibration <- as.matrix(rmvnorm(n = nCal, sigma = covmat,
                                   pre0.9_9994 = T)) # calibration set
  prediction  <- as.matrix(rmvnorm(n = nPred, sigma = covmat,
                                   pre0.9_9994 = T)) # prediction set
  colnames(calibration) <- colnames(prediction) <- colnames(covmat) # assign names to variables
  
  datlist <- list(calibration = calibration, prediction = prediction)
  
  ## add attributes in case you need them later
  attr(datlist, "sampID")     <- ifelse(!is.null(sampID), sampID, NA)
  attr(datlist, "nCal")       <- nCal
  attr(datlist, "nPred")      <- nPred
  attr(datlist, "seed")       <- ifelse(!is.null(seed), seed, NA)
  attr(datlist, "n_x") <- attr(covmat, "n_x")
  attr(datlist, "n_eta_x") <- attr(covmat, "n_eta_x")
  attr(datlist, "n_y") <- attr(covmat, "n_y")
  attr(datlist, "n_eta_y") <- attr(covmat, "n_eta_y")
  attr(datlist, "misspecify") <- attr(covmat, "misspecify")
  attr(datlist, "miss.part") <- attr(covmat, "miss.part")
  attr(datlist, "miss.strength") <- attr(covmat, "miss.strength")
  
  return(datlist)
  
}

# sig <- genCovmat(n_x = 12, n_eta_x = 3, n_y = 4, misspecify = T,
#                  miss.part = "both:load", miss.strength = "weak")
# dat <- gendat(sampID = 1, nCal = 250, covmat = sig)

#----

