## Aditi M. Bhangale
## Last updated: 20 August 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### data generation file -- create data from scratch

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

# n_x = 15; n_eta_x = 3; n_y = 9; n_eta_y = 1
# set.seed(as.numeric(paste0(n_eta_x, n_x, n_eta_y, n_y)))

# misspecification types----
xxrescov <- function(n_x, n_eta_x, theta, miss.strength) {
  if (miss.strength == "weak") {
  # rescov between the second indicator of the first factor 
  # and the second indicator of the second factor
  theta["x2", paste0("x", (n_x/n_eta_x + 2))] <- 
    theta[paste0("x", (n_x/n_eta_x + 2)), "x2"] <- 
    sign(runif(1, min = -1, max = 1)) * runif(1, min = 1, max = 4)
  
  # rescov between the third indicator of the second factor 
  # and the third indicator of the third factor
  theta[paste0("x", (n_x/n_eta_x + 3)), paste0("x", (n_x-n_x/n_eta_x+3))] <- 
    theta[paste0("x", (n_x-n_x/n_eta_x+3)), paste0("x", (n_x/n_eta_x + 3))] <- 
    sign(runif(1, min = -1, max = 1)) * runif(1, min = 1, max = 4)
  } else if (miss.strength == "strong") {
    # rescov between the second indicator of the first factor 
    # and the second indicator of the second factor
    theta["x2", paste0("x", (n_x/n_eta_x + 2))] <- 
      theta[paste0("x", (n_x/n_eta_x + 2)), "x2"] <- 
      sign(runif(1, min = -1, max = 1)) * runif(1, min = 8, max = 12)
    
    # rescov between the third indicator of the second factor 
    # and the third indicator of the third factor
    theta[paste0("x", (n_x/n_eta_x + 3)), paste0("x", (n_x-n_x/n_eta_x+3))] <- 
      theta[paste0("x", (n_x-n_x/n_eta_x+3)), paste0("x", (n_x/n_eta_x + 3))] <- 
      sign(runif(1, min = -1, max = 1)) * runif(1, min = 8, max = 12)
  } else stop("specify a valid misspecification strength in `miss.strength`")
  
  return(theta)
}

xxcrossload <- function(n_x, n_eta_x, lambda, miss.strength) {
  if(miss.strength == "weak") {
    # loading between first factor and second indicator of second factor
    lambda[paste0("x", (n_x/n_eta_x + 2)), "eta_x1"] <- runif(1, min = 0.2, max = 0.4)
    
    # loading between second factor and third indicator of third factor
    lambda[paste0("x", (n_x-n_x/n_eta_x+3)), "eta_x2"] <- runif(1, min = 0.2, max = 0.4)
      
  } else if (miss.strength == "strong") {
    # loading between first factor and second indicator of second factor
    lambda[paste0("x", (n_x/n_eta_x + 2)), "eta_x1"] <- runif(1, min = 0.8, max = 1)
    
    # loading between second factor and third indicator of third factor
    lambda[paste0("x", (n_x-n_x/n_eta_x+3)), "eta_x2"] <- runif(1, min = 0.8, max = 1)
  } else stop("specify a valid misspecification strength in `miss.strength`")

  return(lambda)
}

xydirect <- function(n_x, n_eta_x, n_y, n_eta_y, lambda, beta, psi, theta, 
                     miss.strength, lvnames, obsnames) {
  if (n_eta_x == 1) {
    # add new column to lambda - create a new single-indicator factor
    mis.lambda <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x, 
                         dimnames = list(obsnames,"mis.eta_x1"))
    mis.lambda["x1", "mis.eta_x1"] <- 1
    mis.lambda <- cbind(lambda, mis.lambda)
    
    # set residual variance to 0 (due to new single-indicator factor)
    theta["x1", "x1"] <- 0
    
    # add residual variance for the new single-indicator factor
    mis.psi.val <- ifelse(miss.strength == "weak", 
                          runif(1, min = 5, max = 20),
                          runif(1, min = 20, max = 50))
    mis.psi <- rbind(cbind(psi, rep(0, length(lvnames))), 
                     c(rep(0, length(lvnames)), mis.psi.val))
    dimnames(mis.psi) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    # add regression slope of outcome factor on new single-indicator factor
    mis.beta.val <- ifelse(miss.strength == "weak", 
                        sign(runif(1, min = -1, max = 1)) * runif(1, min = 0.15, max = 0.25),
                             sign(runif(1, min = -1, max = 1)) * runif(1, min = 0.25, max = 0.40))
    mis.beta <- rbind(cbind(beta, rep(0, length(lvnames))), 
                  c(rep(0, length(lvnames)), mis.beta.val))
    dimnames(mis.beta) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    
  } else if (n_eta_x == 3) {
    mis.lvnames <- c(lvnames, paste0("mis.eta_x", 
                                     c(1, 
                                       (n_x/n_eta_x + 1), 
                                       (n_x-n_x/n_eta_x+1))))
    
    # add new columns to lambda -- create single-indicator factors
    mis.lambda <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x,
                         dimnames = list(obsnames, mis.lvnames[grep("mis.", mis.lvnames)]))
    mis.lambda["x1", "mis.eta_x1"] <- 1
    mis.lambda[paste0("x", (n_x/n_eta_x + 1)), 
               paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 1
    mis.lambda[paste0("x", (n_x-n_x/n_eta_x+1)), 
               paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 1
    
    mis.lambda <- cbind(lambda, mis.lambda)
    
    # set residual variances to 0
    theta["x1", "x1"] <- 0
    theta[paste0("x", (n_x/n_eta_x + 1)), paste0("x", (n_x/n_eta_x + 1))] <- 0
    theta[paste0("x", (n_x-n_x/n_eta_x+1)),  paste0("x", (n_x-n_x/n_eta_x+1))] <- 0
    
    # add residual variance for new single-indicator factors
    mis.psi.vals <- if (miss.strength == "weak") {
      runif(n_eta_x, min = 5, max = 20)
    } else if (miss.strength == "strong") {
      runif(n_eta_x, min = 20, max = 50) 
    }
    names(mis.psi.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    
    mis.psi <- matrix(0, nrow = 2*n_eta_x + n_eta_y, ncol = 2*n_eta_x + n_eta_y,
                      dimnames = list(mis.lvnames, mis.lvnames))
    mis.psi[lvnames, lvnames] <- psi
    mis.psi["mis.eta_x1", "mis.eta_x1"] <- mis.psi.vals["mis.eta_x1"]
    mis.psi[paste0("mis.eta_x", (n_x/n_eta_x + 1)), 
            paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- mis.psi.vals[paste0("mis.eta_x", (n_x/n_eta_x + 1))]
    mis.psi[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1)), 
            paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- mis.psi.vals[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))]
    
    # add regression slopes of outcome factor on new single-indicator factors
    mis.beta.vals <- if (miss.strength == "weak") {
      sign(runif(n_eta_x, min = -1, max = 1)) * runif(n_eta_x, min = 0.15, max = 0.25)
    } else if (miss.strength == "strong") {
      sign(runif(n_eta_x, min = -1, max = 1)) * runif(n_eta_x, min = 0.25, max = 0.40)
    }
    names(mis.beta.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    
    mis.beta <- matrix(0, nrow = 2*n_eta_x + n_eta_y, ncol = 2*n_eta_x + n_eta_y,
                      dimnames = list(mis.lvnames, mis.lvnames))
    mis.beta[lvnames, lvnames] <- beta
    mis.beta["eta_y1", "mis.eta_x1"] <- mis.beta.vals["mis.eta_x1"]
    mis.beta["eta_y1", paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 
      mis.beta.vals[paste0("mis.eta_x", (n_x/n_eta_x + 1))]
    mis.beta["eta_y1", paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 
      mis.beta.vals[paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))]
    
  }
  
  # recompute phi
  Iden <- diag(1, nrow = nrow(mis.beta))
  mis.phi <- solve(Iden - mis.beta) %*% mis.psi %*% t(solve(Iden - mis.beta))
  
  return(list(mis.lambda = mis.lambda, mis.theta = theta, 
              mis.psi = mis.psi, mis.B = mis.beta, mis.phi = mis.phi))
}

#----

# generate random covariance matrices----
genCovmat <- function(n_x, n_eta_x, n_y,  n_eta_y, misspecify,
                      miss.part = NULL, 
                      miss.strength = NULL) {
  
  #TODO add check for if n_x is divisible by n_eta_x and n_y is divisible by n_eta_y
  #TODO add check for if n_eta_y is 1 -- only n_eta_y = 1 supported for now
  
  
  obsnames <- c(paste0("x", 1:n_x), paste0("y", 1:n_y)) # indicator labels
  lvnames <- c(paste0("eta_x", 1:n_eta_x), paste0("eta_y", 1:n_eta_y)) # factor labels
  
  set.seed(as.numeric(paste0(n_eta_x, n_x, n_eta_y, n_y)))
  
  # factor loadings
  LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x + n_eta_y, 
               dimnames = list(obsnames, lvnames))
  
  if (n_eta_x == 1L) {
    LAMBDA[1:n_x, paste0("eta_x1")] <- runif(n_x, min = 0.2, max = 1)
    } else if (n_eta_x == 3L) {
      for (xx in 1:n_eta_x) { # only conditions wherein all factors have equal number of indicators supported for now
        if(xx == 1L) {
          LAMBDA[paste0("x", 1:(n_x/n_eta_x)), paste0("eta_x", xx)] <- 
            runif(n_x/n_eta_x, min = 0.2, max = 1)
        } else if (xx == 2L) {
          LAMBDA[paste0("x", (n_x/n_eta_x + 1):(n_x-n_x/n_eta_x)), paste0("eta_x", xx)] <- 
            runif(n_x/n_eta_x, min = 0.2, max = 1)
        } else if (xx == 3L) {
          LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1):n_x), paste0("eta_x", xx)] <- 
            runif(n_x/n_eta_x, min = 0.2, max = 1)
        }
      }
    } else stop("only `n_eta_x = 1 or 3` supported for now")
  
  if (n_eta_y == 1L) {
    LAMBDA[paste0("y",1:n_y), "eta_y1"] <- runif(n_y, min = 0.2, max = 1)
    } else stop("only `n_eta_y = 1` supported for now")
  
  # factor covariances - ignoring structural part
  PSI <- crossprod(matrix(rnorm((n_eta_x + n_eta_y)^2, mean = 3, sd = 4),
                          nrow = n_eta_x + n_eta_y)) # create some arbitrary PSI matrix
  dimnames(PSI) <- list(lvnames, lvnames)
  # PSI; cov2cor(PSI)
  
  # regression coefficients - structural part
  B <- matrix(0, nrow = n_eta_x + n_eta_y, 
              ncol = n_eta_x + n_eta_y, 
              dimnames = list(lvnames, lvnames)) # rows are outcomes, columns are predictors
  B[paste0("eta_y", 1:n_eta_y), paste0("eta_x", 1:n_eta_x)] <- 
    sign(runif(n_eta_x, min = -1, max = 1)) * runif(n_eta_x, min = 0.15, max = 0.40)
  
  # factor covariances - incorporating information from structural part
  Iden <- diag(1, nrow = n_eta_x + n_eta_y)
  PHI <- solve(Iden - B) %*% PSI %*% t(solve(Iden - B))
  # PHI; cov2cor(PHI)
  
  # residual (co)variances of observed variables
  THETA <- diag(runif(n_x + n_y, min = 10, max = 50))
  dimnames(THETA) <- list(obsnames, obsnames)
  
  # introduce misspecification
  if(misspecify == T) {
    if (is.null(miss.part) | is.null(miss.strength)) {
      stop("specify which part of the model to misspecify and strength of misspecification")
    } else {
      # miss.part = c("xx:rescov", "xx:crossload", "xy:direct", 
      #               "both:cov", "both:load")
      if (miss.part == "xx:rescov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, theta = THETA, miss.strength = miss.strength)
      } else if (miss.part == "xx:crossload") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = LAMBDA, miss.strength = miss.strength)
      } else if (miss.part == "xy:direct") {
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                lambda = LAMBDA, beta = B, psi = PSI, 
                                theta = THETA, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.lambda
        PHI <- direct.vals$mis.phi
        THETA <- direct.vals$mis.theta
      } else if (miss.part == "both:cov") {
        THETA <- xxrescov(n_x = n_x, n_eta_x = n_eta_x, theta = THETA, miss.strength = miss.strength)
        
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                lambda = LAMBDA, beta = B, psi = PSI, 
                                theta = THETA, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.lambda
        PHI <- direct.vals$mis.phi
        THETA <- direct.vals$mis.theta
      } else if (miss.part == "both:load") {
        LAMBDA <- xxcrossload(n_x = n_x, n_eta_x = n_eta_x, lambda = LAMBDA, miss.strength = miss.strength)
        
        direct.vals <- xydirect(n_x = n_x, n_eta_x = n_eta_x, n_y = n_y, n_eta_y = n_eta_y, 
                                lambda = LAMBDA, beta = B, psi = PSI, 
                                theta = THETA, miss.strength = miss.strength,
                                lvnames = lvnames, obsnames = obsnames)
        LAMBDA <- direct.vals$mis.lambda
        PHI <- direct.vals$mis.phi
        THETA <- direct.vals$mis.theta
      } else stop("specify valid `miss.part`")
    }
  }
  
  # population covariance matrix
  SIGMA.pop <- LAMBDA %*% PHI %*% t(LAMBDA) + THETA 
  R.pop <- cov2cor(SIGMA.pop)
  
  attr(SIGMA.pop, "n_x") <- n_x
  attr(SIGMA.pop, "n_eta_x") <- n_eta_x
  attr(SIGMA.pop, "n_y") <- n_y
  attr(SIGMA.pop, "n_eta_y") <- n_eta_y
  attr(SIGMA.pop, "misspecify") <- misspecify
  attr(SIGMA.pop, "miss.part") <- ifelse(!is.null(miss.part), miss.part, NA)
  attr(SIGMA.pop, "miss.strength") <- ifelse(!is.null(miss.strength), miss.strength, NA)
  
  return(list(SIGMA.pop = SIGMA.pop, R.pop = R.pop))
}

# foo <- genCovmat(n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1, misspecify = T,
#                  miss.part = "both:load", miss.strength = "strong")$SIGMA.pop
# dim(foo)

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

# sigma <- genCovmat(n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1, misspecify = T, 
#                     miss.part = "both:cov", miss.strength = "strong")$SIGMA.pop
# 
# dat <- gendat(sampID = 1, nCal = 250, nPred = 1e5, covmat = sigma)

#----