## Aditi M. Bhangale
## Last updated: 21 April 2025

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
mod1 <- '
# factor loadings
F1 =~ 0.67*x1 + 0.71*x2 + 0.58*x3 + 0.83*x4 + 0.64*x5 + 0.55*x6 + 0.68*x7

# factor variance
F1 ~~ 1*F1

# item (co)variances
x1 ~~ 0.22*x1
x2 ~~ 0.41*x2
x3 ~~ 0.55*x3
x4 ~~ 0.49*x4
x5 ~~ 0.71*x5
x6 ~~ 0.39*x6
x7 ~~ 0.61*x7

# factor mean
F1 ~ 0*1

# item means
x1 ~ 0.91*1
x2 ~ 0.43*1
x3 ~ 0.73*1
x4 ~ 0.44*1
x5 ~ 0.83*1
x6 ~ 0.66*1
x7 ~ 0.54*1
'

mod2 <- '
# factor loadings
F1 =~ 0.67*x1 + 0.71*x2 + 0.58*x3
F2 =~ 0.83*x4 + 0.64*x5 + 0.55*x6 + 0.68*x7

# factor (co)variances
F1 ~~ 1*F1 + 0.3*F2
F2 ~~ 1*F2

# item (co)variances
x1 ~~ 0.22*x1
x2 ~~ 0.41*x2
x3 ~~ 0.55*x3
x4 ~~ 0.49*x4
x5 ~~ 0.71*x5
x6 ~~ 0.39*x6
x7 ~~ 0.61*x7

# factor means
F1 ~ 0*1
F2 ~ 0*1

# item means
x1 ~ 0.91*1
x2 ~ 0.43*1
x3 ~ 0.73*1
x4 ~ 0.44*1
x5 ~ 0.83*1
x6 ~ 0.66*1
x7 ~ 0.54*1
'

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
    lavaan(mod1, do.fit = F)
  } else {
    lavaan(mod2, do.fit = F)
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

