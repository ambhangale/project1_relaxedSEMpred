## Aditi M. Bhangale
## Last updated: 1 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### data generation file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(mvtnorm)
library(lavaan)
library(portableParallelSeeds)

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

gendat <- function(sampID = NULL, nCal, nPred, misspecify, seed) {
  
  if (!is.null(sampID)) {
    # create multiple seeds, so that same data is generated each time when a specific `sampID` is used
    allSeeds <- seedCreator(nReps = 5e3, streamsPerRep = 1, seed = seed)
    
    setSeeds(projSeeds = allSeeds, run = sampID) # set seed based on `sampID`
  } else {
    set.seed(seed)
  }
  
  # (do not) fit model, return start values
  fit <- if(!misspecify) {
    lavaan(mod1, do.fit = F)
  } else {
    lavaan(mod2, do.fit = F)
  }
  popStats <- lavInspect(fit, "implied") # population covariance matrix and mean vector
  
  # data generation
  calibration <- as.data.frame(rmvnorm(n = nCal, mean = popStats$mean, sigma = popStats$cov,
                                 pre0.9_9994 = T)) # calibration set
  prediction  <- as.data.frame(rmvnorm(n = nPred, mean = popStats$mean, sigma = popStats$cov,
                                 pre0.9_9994 = T)) # prediction set
  
  datlist <- list(calibration = calibration, prediction = prediction)
  
  ## add attributes in case you need them later
  attr(datlist, "sampID")     <- ifelse(!is.null(sampID), sampID, NA)
  attr(datlist, "nCal")       <- nCal
  attr(datlist, "nPred")      <- nPred
  attr(datlist, "misspecify") <- misspecify
  attr(datlist, "seed")       <- seed
  
  return(datlist)
}

# test function
# gendat(sampID = 1, nCal = 20, nPred = 25, misspecify = F, seed = 10824)
# gendat(sampID = 1, nCal = 25, nPred = 20, misspecify = T, seed = 10824)


