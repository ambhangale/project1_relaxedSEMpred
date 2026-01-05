## Aditi M. Bhangale
## Last updated: 5 January 2026

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
  theta.val <- ifelse(miss.strength == "weak",
                      0.3*min(THETA.star),
                      0.6*min(THETA.star)) # value of the residual covariance
  
  if (n_eta_x == 1L) {
    RC <- if (n_x == 4L) list(c(2,3)) else list(c(2,3), c(4,5), c(6,7), c(2,7)) # residual covariances to introduce
    
    for (rc in 1:length(RC)) {
      THETA[paste0("x",RC[[rc]][1]), paste0("x",RC[[rc]][2])] <-
        THETA[paste0("x",RC[[rc]][2]), paste0("x",RC[[rc]][1])] <- theta.val
    }
  } else {
    # if n_eta_x = 3L,
    ## if n_x = 12L, rescov between 2nd and 3rd indicators of each factor
    ## if n_x = 24L, rescov between 2nd,3rd,7th,8th indicators of each factor
    
    RC <- if (n_x == 12L) c(2,3) else c(2,3,7,8)
    
    for(rc in RC) {
      THETA[paste0("x",rc),paste0("x",(n_x/n_eta_x+rc))] <- 
        THETA[paste0("x",(n_x/n_eta_x+rc)),paste0("x",rc)] <-
        THETA[paste0("x",rc),paste0("x",(n_x-n_x/n_eta_x+rc))] <- 
        THETA[paste0("x",(n_x-n_x/n_eta_x+rc)),paste0("x",rc)] <-
        THETA[paste0("x",(n_x/n_eta_x+rc)),paste0("x",(n_x-n_x/n_eta_x+rc))] <-
        THETA[paste0("x",(n_x-n_x/n_eta_x+rc)),paste0("x",(n_x/n_eta_x+rc))] <- theta.val
    }
    
  }
  return(THETA)
}

xxcrossload <- function(n_x, n_eta_x, lambda, LAMBDA, miss.strength) {
  if (n_eta_x != 3L) stop("xx misspecification not supported for n_eta_x = 1L")
  if (!(n_x == 12L || n_x == 24L)) stop("xxcrossload misspecificaiton only supported for n_x == 12L | 24L for now")
  
  lambda.val <- ifelse(miss.strength == "weak", 0.5*lambda, 0.9*lambda)
  
  CL <- c(2,7)
  
  for (xx in 1:n_eta_x) {
    if (xx == 1L) {
      LAMBDA[paste0("x", n_x/n_eta_x+CL[1]),  paste0("eta_x", xx)] <- lambda.val
      if (n_x == 24L) LAMBDA[paste0("x", n_x-n_x/n_eta_x+CL[2]), paste0("eta_x", xx)] <- lambda.val
    } else if (xx == 2L) {
      LAMBDA[paste0("x", n_x-n_x/n_eta_x+CL[1]), paste0("eta_x", xx)] <- lambda.val
      if (n_x == 24L) LAMBDA[paste0("x", CL[2]),  paste0("eta_x", xx)] <- lambda.val
    } else if (xx == 3L) {
      LAMBDA[paste0("x", CL[1]),  paste0("eta_x", xx)] <- lambda.val
      if (n_x == 24L) LAMBDA[paste0("x", n_x/n_eta_x+CL[2]),  paste0("eta_x", xx)] <- lambda.val
    }
  }
  
  return(LAMBDA)
}

xydirect <- function(n_x, n_eta_x, n_y, n_eta_y, LAMBDA, B, PSI, THETA, 
                     beta, r, obs.var, miss.strength, lvnames, obsnames) {

  mis.B.val <- ifelse(miss.strength == "weak", 0.5*beta, 0.9*beta) # structural coefficient for direct effects
  
  if (n_eta_x == 1L) {
    DE <- if (n_x == 4L) 1 else c(1,8)
    mis.lvnames <- c(lvnames,paste0("mis.eta_x", DE))
    mis.LAMBDA.val <- c() # vector of lambda values
    
    # create mis.LAMBDA matrix for new single indicator factor(s)
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, 
                         ncol = n_eta_x + n_eta_y + length(DE), 
                         dimnames = list(obsnames, mis.lvnames))
    mis.LAMBDA[rownames(LAMBDA), colnames(LAMBDA)] <- LAMBDA 
    
    # create mis.PSI and mis.B matrices
    mis.PSI <- mis.B <- matrix(0, nrow = n_eta_x + n_eta_y + length(DE), 
                               ncol = n_eta_x + n_eta_y + length(DE), 
                               dimnames = list(mis.lvnames, mis.lvnames))
    
    mis.PSI[rownames(PSI), colnames(PSI)] <- PSI 
    mis.B[rownames(B), colnames(B)] <- B 
    
    for (de in DE) {
      mis.LAMBDA[paste0("x", de), paste0("mis.eta_x", de)] <- 1L
      mis.LAMBDA.val[which(DE==de)] <- mis.LAMBDA[paste0("x", de), paste0("eta_x1")] # this lambda value will now be in the B matrix
      mis.LAMBDA[paste0("x", de), paste0("eta_x1")] <- 0L
      
      # set single-indicator factor variance to residual variance of the indicator
      # then set residual variance of indicator to zero
      mis.PSI[paste0("mis.eta_x", de), paste0("mis.eta_x", de)] <- 
        THETA[paste0("x", de), paste0("x", de)]
      THETA[paste0("x", de), paste0("x", de)] <- 0L
      
      # add regression slope of outcome factor on new single-indicator factor
      # and the factor loading on the single-indicator factors
      ## rows are outcomes, columns are predictors
      mis.B[paste0("mis.eta_x", de), "eta_x1"] <- mis.LAMBDA.val[which(DE==de)]
      mis.B["eta_y1", paste0("mis.eta_x", de)] <- mis.B.val 
    }
    
  } else if (n_eta_x == 3L) {
    DE <- if (n_x == 12L) c(1,4) else c(1,4,5)
    mis.lvnames <- c(lvnames, 
                     paste0("mis.eta_x", c(DE, (n_x/n_eta_x+DE), (n_x-n_x/n_eta_x+DE))))
    mis.LAMBDA.val <- vector(mode = "list", length = n_eta_x)
    
    # create mis.LAMBDA matrix for new single indicator factor(s)
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, 
                         ncol = n_eta_x + n_eta_y + n_eta_x*length(DE), 
                         dimnames = list(obsnames, mis.lvnames))
    mis.LAMBDA[rownames(LAMBDA), colnames(LAMBDA)] <- LAMBDA 
    
    # create mis.PSI and mis.B matrices
    mis.PSI <- mis.B <- matrix(0, nrow = n_eta_x + n_eta_y + n_eta_x*length(DE), 
                               ncol = n_eta_x + n_eta_y + n_eta_x*length(DE), 
                               dimnames = list(mis.lvnames, mis.lvnames))
    
    mis.PSI[rownames(PSI), colnames(PSI)] <- PSI 
    mis.B[rownames(B), colnames(B)] <- B 
    
    for (de in DE) {
      for (xx in 1:n_eta_x) {
        if (xx == 1L) {
          mis.LAMBDA[paste0("x", de), paste0("mis.eta_x", de)] <- 1L
          mis.LAMBDA.val[[xx]][which(DE==de)] <- mis.LAMBDA[paste0("x", de), paste0("eta_x",xx)] 
          mis.LAMBDA[paste0("x", de), paste0("eta_x", xx)] <- 0L
          
          mis.PSI[paste0("mis.eta_x", de), paste0("mis.eta_x", de)] <- 
            THETA[paste0("x", de), paste0("x", de)]
          THETA[paste0("x", de), paste0("x", de)] <- 0L
          
          mis.B[paste0("mis.eta_x", de), paste0("eta_x",xx)] <- mis.LAMBDA.val[[xx]][which(DE==de)]
          mis.B["eta_y1", paste0("mis.eta_x", de)] <- mis.B.val 
        } else if (xx == 2L) {
          mis.LAMBDA[paste0("x", (n_x/n_eta_x+de)), paste0("mis.eta_x", (n_x/n_eta_x+de))] <- 1L
          mis.LAMBDA.val[[xx]][which(DE==de)] <- mis.LAMBDA[paste0("x", (n_x/n_eta_x+de)), paste0("eta_x",xx)] 
          mis.LAMBDA[paste0("x", (n_x/n_eta_x+de)), paste0("eta_x", xx)] <- 0L
          
          mis.PSI[paste0("mis.eta_x", (n_x/n_eta_x+de)), paste0("mis.eta_x", (n_x/n_eta_x+de))] <- 
            THETA[paste0("x", (n_x/n_eta_x+de)), paste0("x", (n_x/n_eta_x+de))]
          THETA[paste0("x", (n_x/n_eta_x+de)), paste0("x", (n_x/n_eta_x+de))] <- 0L
          
          mis.B[paste0("mis.eta_x", (n_x/n_eta_x+de)), paste0("eta_x",xx)] <- mis.LAMBDA.val[[xx]][which(DE==de)]
          mis.B["eta_y1", paste0("mis.eta_x", (n_x/n_eta_x+de))] <- mis.B.val 
        } else if (xx == 3L) {
          mis.LAMBDA[paste0("x", (n_x-n_x/n_eta_x+de)), paste0("mis.eta_x", (n_x-n_x/n_eta_x+de))] <- 1L
          mis.LAMBDA.val[[xx]][which(DE==de)] <- mis.LAMBDA[paste0("x", (n_x-n_x/n_eta_x+de)), paste0("eta_x",xx)] 
          mis.LAMBDA[paste0("x", (n_x-n_x/n_eta_x+de)), paste0("eta_x", xx)] <- 0L
          
          mis.PSI[paste0("mis.eta_x", (n_x-n_x/n_eta_x+de)), paste0("mis.eta_x", (n_x-n_x/n_eta_x+de))] <- 
            THETA[paste0("x", (n_x-n_x/n_eta_x+de)), paste0("x", (n_x-n_x/n_eta_x+de))]
          THETA[paste0("x", (n_x-n_x/n_eta_x+de)), paste0("x", (n_x-n_x/n_eta_x+de))] <- 0L
          
          mis.B[paste0("mis.eta_x", (n_x-n_x/n_eta_x+de)), paste0("eta_x",xx)] <- mis.LAMBDA.val[[xx]][which(DE==de)]
          mis.B["eta_y1", paste0("mis.eta_x", (n_x-n_x/n_eta_x+de))] <- mis.B.val 
        }
      }
    }
  }
  
  # recompute PHI
  Iden <- diag(1, nrow = nrow(mis.B))
  mis.PHI <- solve(Iden - mis.B) %*% mis.PSI %*% t(solve(Iden - mis.B))
  
  return(list(mis.LAMBDA = mis.LAMBDA, mis.THETA = THETA, 
              mis.PSI = mis.PSI, mis.B = mis.B, mis.PHI = mis.PHI))
}

#----

# generate random covariance matrices----
genCovmat <- function(n_x, n_eta_x, n_y,  n_eta_y = 1L, 
                      beta = 0.3, lambda = 0.7, psi.cov = 0.2, r = 0.3, obs.var = 1,
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
  PSI[lvnames[grep("_y", lvnames)], lvnames[grep("_y", lvnames)]] <- ifelse(n_y == 1L, r*obs.var, 1-beta^2)
  ## fix factor variance to r*obs.var if y is a single-indicator factor 
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
    LAMBDA[paste0("y",1:n_y), "eta_y1"] <- ifelse(n_y == 1L, 1L, lambda)
    ## fix factor loading for single-indicator factor to 1L if y is a single-indicator factor 
  
  # PHI; calculate from LAMBDA and B, factor covariances incorporating structural part
  Iden <- diag(1, nrow = n_eta_x + n_eta_y)
  PHI <- solve(Iden - B) %*% PSI %*% t(solve(Iden - B))
  
  # THETA; residual variances of indicators
  THETA.dash <- diag(LAMBDA%*%PHI%*%t(LAMBDA))
  THETA.star <- THETA.dash*(1/r - 1)
  THETA <- diag(THETA.star, nrow = n_x + n_y)
  dimnames(THETA) <- list(obsnames, obsnames)
  if(n_y == 1L) THETA["y1", "y1"] <- (1-r)*obs.var
  ## fix residual indicator variance to (1-r)*obs.var if y is a single-indicator factor 
  
  # introduce misspecification
  if(misspecify == T) {
    if (is.null(miss.part) | is.null(miss.strength)) {
      stop("specify which part of the model to misspecify and strength of misspecification")
    } else {
      if(!(miss.strength == "weak" || miss.strength == "strong")) stop("specify a valid misspecification strength in `miss.strength`") 
      if(!(miss.part %in% c("xx:rescov", "xx:crossload", "xy:direct", "both:cov", "both:load"))) stop("specify a valid misspecification type in `miss.part`")
      
      if (miss.part == "xx:rescov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, THETA = THETA, 
                          THETA.star = THETA.star, miss.strength = miss.strength)
      } else if (miss.part == "xx:crossload") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = lambda, 
                              LAMBDA = LAMBDA, miss.strength = miss.strength)
      } else if (miss.part == "xy:direct") {
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta = beta, r = r, obs.var = obs.var,
                                miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames) 
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
      } else if (miss.part == "both:cov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, THETA = THETA, 
                          THETA.star = THETA.star, miss.strength = miss.strength)
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta = beta, r = r, obs.var = obs.var, 
                                miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
        
      } else if (miss.part == "both:load") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = lambda, 
                              LAMBDA = LAMBDA, miss.strength = miss.strength)
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                LAMBDA = LAMBDA, B = B, PSI = PSI, 
                                THETA = THETA, beta = beta, r = r, obs.var = obs.var, 
                                miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames) 
        LAMBDA <- direct.vals$mis.LAMBDA
        PHI <- direct.vals$mis.PHI
        THETA <- direct.vals$mis.THETA
        
      } else stop("specify valid `miss.part`")
    }
  }
  if (!(all(eigen(THETA)$values >= -3.14e-14))) warning ("THETA matrix may not be positive (semi-)definite")
  ## above, use an arbitrary small number because eigendecomposition is subject to errors on real-world computers.
  ## even though the true eigenvalue is 0, it may be computed as a negligible negative value due to these errors.
  if (!(all(eigen(PHI)$values >= 0))) warning ("PHI matrix may not be positive (semi-)definite")
  
  # population covariance matrix
  SIGMA.pop <- LAMBDA %*% PHI %*% t(LAMBDA) + THETA
  if (!(all(eigen(SIGMA.pop)$values >= 0))) stop ("SIGMA matrix is not positive (semi-)definite")
  
  attr(SIGMA.pop, "n_x") <- n_x
  attr(SIGMA.pop, "n_eta_x") <- n_eta_x
  attr(SIGMA.pop, "n_y") <- n_y
  attr(SIGMA.pop, "n_eta_y") <- n_eta_y
  attr(SIGMA.pop, "beta") <- beta
  attr(SIGMA.pop, "lambda") <- lambda
  attr(SIGMA.pop, "psi.cov") <- psi.cov
  attr(SIGMA.pop, "r") <- r
  attr(SIGMA.pop, "obs.var") <- obs.var
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

