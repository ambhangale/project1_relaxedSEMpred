## Aditi M. Bhangale
## Last updated: 2 October 2025

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

# n_x = 9; n_eta_x = 3; n_y = 9; n_eta_y = 1
# beta = 0.3; lambda = 0.7; psi.cov = 0.2; r = 0.3
# miss.strength = "strong"

# misspecification types----
xxrescov <- function(n_x, n_eta_x, THETA, THETA.star, miss.strength) {
  if (miss.strength == "weak") {
    # rescov between the second indicator of the first factor 
    # and the second indicator of the second factor
    # AND
    # rescov between the third indicator of the second factor 
    # and the third indicator of the third factor
    
    THETA["x2", paste0("x", (n_x/n_eta_x + 2))] <-  THETA[paste0("x", (n_x/n_eta_x + 2)), "x2"] <- 
      THETA[paste0("x", (n_x/n_eta_x + 3)), paste0("x", (n_x-n_x/n_eta_x+3))] <- 
      THETA[paste0("x", (n_x-n_x/n_eta_x+3)), paste0("x", (n_x/n_eta_x + 3))] <- 
      0.3*min(THETA.star)
    
  } else if (miss.strength == "strong") {
    # rescov between the second indicator of the first factor 
    # and the second indicator of the second factor
    # AND
    # rescov between the third indicator of the second factor 
    # and the third indicator of the third factor
    
    THETA["x2", paste0("x", (n_x/n_eta_x + 2))] <- THETA[paste0("x", (n_x/n_eta_x + 2)), "x2"] <- 
      THETA[paste0("x", (n_x/n_eta_x + 3)), paste0("x", (n_x-n_x/n_eta_x+3))] <- 
      THETA[paste0("x", (n_x-n_x/n_eta_x+3)), paste0("x", (n_x/n_eta_x + 3))] <- 
      0.6*min(THETA.star)
  } else stop("specify a valid misspecification strength in `miss.strength`")
  
  return(THETA)
}

xxcrossload <- function(n_x, n_eta_x, lambda, LAMBDA, miss.strength) {
  if(miss.strength == "weak") {
    # loading between first factor and second indicator of second factor
    # AND
    # loading between second factor and third indicator of third factor
    LAMBDA[paste0("x", (n_x/n_eta_x + 2)), "eta_x1"] <- 
      LAMBDA[paste0("x", (n_x-n_x/n_eta_x+3)), "eta_x2"] <- 0.5*lambda
    
  } else if (miss.strength == "strong") {
    # loading between first factor and second indicator of second factor
    # AND
    # loading between second factor and third indicator of third factor
    
    LAMBDA[paste0("x", (n_x/n_eta_x + 2)), "eta_x1"] <- 
      LAMBDA[paste0("x", (n_x-n_x/n_eta_x+3)), "eta_x2"] <- 0.9*lambda
    
  } else stop("specify a valid misspecification strength in `miss.strength`")
  
  return(LAMBDA)
}

xydirect <- function(n_x, n_eta_x, n_y, n_eta_y, LAMBDA, B, PSI, THETA, 
                     beta, miss.strength, lvnames, obsnames) {
  
  # FIXME --- can only add new factor for second indicator of each factor!!!!
  
  if (n_eta_x == 1) {
    # add new column to LAMBDA - create a new single-indicator factor
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x, 
                         dimnames = list(obsnames,"mis.eta_x1"))
    mis.LAMBDA["x1", "mis.eta_x1"] <- 1
    mis.LAMBDA <- cbind(LAMBDA, mis.LAMBDA)
    
    # set residual variance to 0 (due to new single-indicator factor)
    THETA["x1", "x1"] <- 0
    
    # add residual variance for the new single-indicator factor
    mis.PSI.val <- 1L # all exogenous factors have factor variance of 1
    #FIXME should these single indcator factor (above). have a smaller factor variance?
    mis.PSI <- rbind(cbind(PSI, rep(0, length(lvnames))), 
                     c(rep(0, length(lvnames)), mis.PSI.val))
    dimnames(mis.PSI) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    # add regression slope of outcome factor on new single-indicator factor
    mis.B.val <- ifelse(miss.strength == "weak", 
                           0.5*beta,
                           0.9*beta)
    mis.B <- rbind(cbind(B, rep(0, length(lvnames))), 
                      c(rep(0, length(lvnames)), mis.BETA.val))
    dimnames(mis.B) <- list(c(lvnames,"mis.eta_x1"), c(lvnames,"mis.eta_x1"))
    
    
  } else if (n_eta_x == 3) {
    mis.lvnames <- c(lvnames, paste0("mis.eta_x", 
                                     c(1, 
                                       (n_x/n_eta_x + 1), 
                                       (n_x-n_x/n_eta_x+1))))
    
    # add new columns to LAMBDA -- create single-indicator factors
    mis.LAMBDA <- matrix(0, nrow = n_x + n_y, ncol = n_eta_x,
                         dimnames = list(obsnames, mis.lvnames[grep("mis.", mis.lvnames)]))
    mis.LAMBDA["x1", "mis.eta_x1"] <- 1
    mis.LAMBDA[paste0("x", (n_x/n_eta_x + 1)), 
               paste0("mis.eta_x", (n_x/n_eta_x + 1))] <- 1
    mis.LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1)), 
               paste0("mis.eta_x", (n_x-n_x/n_eta_x+1))] <- 1
    
    mis.LAMBDA <- cbind(LAMBDA, mis.LAMBDA)
    
    # set residual variances to 0
    THETA["x1", "x1"] <- 0
    THETA[paste0("x", (n_x/n_eta_x + 1)), paste0("x", (n_x/n_eta_x + 1))] <- 0
    THETA[paste0("x", (n_x-n_x/n_eta_x+1)),  paste0("x", (n_x-n_x/n_eta_x+1))] <- 0
    
    # add residual variance for new single-indicator factors
    mis.PSI.vals <- rep(1L, n_eta_x) # all exogenous variables have factor variance of 1
    #FIXME should these single indcator factors (above). have a smaller factor variance?
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
    
    # add regression slopes of outcome factor on new single-indicator factors
    mis.B.vals <- rep(ifelse(miss.strength == "weak", 0.5*beta, 0.9*beta), n_eta_x)
    names(mis.B.vals) <- mis.lvnames[grep("mis.", mis.lvnames)]
    
    mis.B <- matrix(0, nrow = 2*n_eta_x + n_eta_y, ncol = 2*n_eta_x + n_eta_y,
                       dimnames = list(mis.lvnames, mis.lvnames))
    mis.B[lvnames, lvnames] <- B
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
                      misspecify, miss.part = NULL, miss.strength = NULL) { #FIXME add more arguments as necessary
  
  # check if n_x is divisible by n_eta_x, otherwise stop
  if (n_x %% n_eta_x != 0) stop("`n_eta_x` must be divisble by `n_x`")
  
  # check if n_eta_y = 1L, otherwise stop
  if (n_eta_y != 1L) stop("only `n_eta_y == 1L` supported for now")
  
  # check if n_eta_x = 1 or 3, otherwise stop
  if (!(n_eta_x == 1L || n_eta_x == 3L)) stop("only `n_eta_x == 1L or 3L` supported for now")
  
  
  obsnames <- c(paste0("x", 1:n_x), paste0("y", 1:n_y)) # indicator labels
  lvnames <- c(paste0("eta_x", 1:n_eta_x), paste0("eta_y", 1:n_eta_y)) # factor labels
  
  # B; matrix of regresion coefficients (structural relations)
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
                   dimnames = list(obsnames, lvnames))
  
  ## for x factors
  if (n_eta_x == 1L) {
    LAMBDA["x1", "eta_x1"] <- 1L
    LAMBDA[paste0("x", 2:n_x), "eta_x1"] <- lambda
  } else if (n_eta_x == 3L) {
    for (xx in 1:n_eta_x) { # only conditions wherein all factors have equal number of indicators supported for now
      if(xx == 1L) {
        LAMBDA["x1", paste0("eta_x", xx)] <- 1L
        LAMBDA[paste0("x", 2:(n_x/n_eta_x)), paste0("eta_x", xx)] <- lambda
      } else if (xx == 2L) {
        LAMBDA[paste0("x", (n_x/n_eta_x + 1)), paste0("eta_x", xx)] <- 1L
        LAMBDA[paste0("x", (n_x/n_eta_x + 2):(n_x-n_x/n_eta_x)), paste0("eta_x", xx)] <- lambda
      } else if (xx == 3L) {
        LAMBDA[paste0("x", (n_x-n_x/n_eta_x+1)), paste0("eta_x", xx)] <- 1L
        LAMBDA[paste0("x", (n_x-n_x/n_eta_x+2):n_x), paste0("eta_x", xx)] <- lambda
      }
    }
  }
  
  ## for y factor
  LAMBDA["y1", "eta_y1"] <- 1L
  LAMBDA[paste0("y",2:n_y), "eta_y1"] <- lambda
  
  # PHI; calculate from LAMBDA and B, factor covariances incorporating structural part
  Iden <- diag(1, nrow = n_eta_x + n_eta_y)
  PHI <- solve(Iden - B) %*% PSI %*% t(solve(Iden - B))
  #FIXME my exogenous factors need to be correlated with one another, right?
  
  # THETA; residual variances of indicators
  THETA.dash <- diag(LAMBDA%*%PHI%*%t(LAMBDA))
  THETA.star <- THETA.dash*(1/r - 1)
  THETA <- diag(THETA.star, nrow = n_x + n_y)
  dimnames(THETA) <- list(obsnames, obsnames)
  
  # introduce misspecification
  if(misspecify == T) {
    if (is.null(miss.part) | is.null(miss.strength)) {
      stop("specify which part of the model to misspecify and strength of misspecification")
    } else {
      # miss.part = c("xx:rescov", "xx:crossload", "xy:direct", 
      #               "both:cov", "both:load")
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
gendat <- function() {
}

# sigma <- genCovmat(n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1, misspecify = T, 
#                     miss.part = "both:cov", miss.strength = "strong")$SIGMA.pop
# 
# dat <- gendat(sampID = 1, nCal = 250, nPred = 1e5, covmat = sigma)

#----