## Aditi M. Bhangale
## Last updated: 8 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### data generation file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(mvtnorm)
library(lavaan)
library(portableParallelSeeds) # remotes::install_github("wjakethompson/portableParallelSeeds")

# create seed object for data generation and partition generation----
allSeeds <- seedCreator(nReps = 5e3, streamsPerRep = 2, seed = 10824)

# possible to run sampIDs 1:5000. adjust nReps if we want to go higher.
# 2 streams. stream 1 to be used for data generation and stream 2 to be used for partitioning
#----

# correctly specified and misspecified models----
data("PoliticalDemocracy")

PoliticalDemocracy$dem65_sum <- apply(PoliticalDemocracy[, paste0("y", 5:8)],
                                      1, sum) # create a sum score for y (to be predicted)

cs.mod <- ' 
  # latent variable definitions
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    
  # regressions
    dem65_sum ~ ind60 + dem60
' # correctly specified model

ms.mod <- ' 
  # latent variable definitions
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    
  # regressions
    dem65_sum ~ ind60 + dem60 + x1 + x2 + x3 + y1 + y2 + y3 + y4
' # misspecified model

# the above model will not produce SEs, possible due to (empirical) non-identification
# but we still use the estimates, because there are no Heywood cases and because
# we will only use the estimates to generate data. this model will NOT be fitted again

#----

# data generation function----
gendat <- function(sampID = NULL, nCal, nPred, misspecify, seed = NULL) {
  
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
  
  # (do not) fit model, return start values
  fit <- if(!misspecify) {
    sem(cs.mod, data = PoliticalDemocracy, meanstructure = T)
  } else {
    sem(ms.mod, data = PoliticalDemocracy, meanstructure = T)
  }
  popStats <- lavInspect(fit, "implied") # population covariance matrix and mean vector
  
  # data generation
  calibration <- as.matrix(rmvnorm(n = nCal, mean = popStats$mean, sigma = popStats$cov,
                                 pre0.9_9994 = T)) # calibration set
  prediction  <- as.matrix(rmvnorm(n = nPred, mean = popStats$mean, sigma = popStats$cov,
                                 pre0.9_9994 = T)) # prediction set
  
  datlist <- list(calibration = calibration, prediction = prediction)
  
  ## add attributes in case you need them later
  attr(datlist, "sampID")     <- ifelse(!is.null(sampID), sampID, NA)
  attr(datlist, "nCal")       <- nCal
  attr(datlist, "nPred")      <- nPred
  attr(datlist, "misspecify") <- misspecify
  attr(datlist, "seed")       <- ifelse(!is.null(seed), seed, NA)
  
  return(datlist)
}

# test function
# gendat(sampID = 1, nCal = 20, nPred = 25, misspecify = F)
# gendat(sampID = 1, nCal = 25, nPred = 20, misspecify = T)
#----

