## Aditi M. Bhangale
## Last updated: 22 April 2025
## Fixes and improvements: 4 April 2025 (Julian D. Karch)

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### functions to create data partitions/folds

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "datagen_RDA_CFA.R"))

# dat <- gendat(sampID = 2, nCal = 250, nPred = 250, misspecify = F, seed = 10824) # dummy data for now
# calibration <- dat$calibration; prediction <- dat$prediction

# function to assign partition IDs per observation----
partidx <- function(ndat, sampID = NULL, K, nK, seed = NULL) {
  if (!is.null(sampID)) {
    # set seed based on `sampID`
    setSeeds(projSeeds = allSeeds, run = sampID)
    useStream(2) # use the second stream of the `sampID` rep for partition creation
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    } else {
      stop("specify `seed` argument if `sampID` is left blank") 
    }
  }
  
  partvec <- rep(1:K, nK)
  
  # randomly assign observations to one of the K parts such that there are nK 
  # observations per part
  part.id <- sample(x = partvec, size = ndat, replace = F)
  
  return(part.id)
  
}

# partidx(ndat = 100, sampID = 2, K = 10, nK = 10)

#----

# function to partition data into K parts----
partition <- function(sampID = NULL, dat, K, nK, seed = NULL) { # input here is the calibration data
  
  partid <- partidx(ndat = nrow(dat), sampID = sampID, K = K, nK = nK, seed = seed)
  
  final <- vector("list", K)
  
  for (k in seq_len(K)) {
    test  <- dat[partid == k, ]
    train <- dat[partid != k, ]
    
    final[[k]] <- list(train = train, test = test)
  }
  
  # names(final) <- paste0("part", seq_len(K)) #FIXME maybe not needed? never used, i think
  
  return(final)
}

# partition(sampID = 2, dat = calibration, K = 10, nK = 25) # test function

#----
