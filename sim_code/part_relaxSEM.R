## Aditi M. Bhangale
## Last updated: 29 May 2025
## Fixes and improvements: 4 April 2025 (Julian D. Karch)

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to create data partitions/folds

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

source("datagen_relaxSEM.R")

# dat <- gendat(sampID = 2, nCal = 250, nPred = 250, misspecify = F, seed = 10824) # dummy data for now
# calibration <- dat$calibration; prediction <- dat$prediction

# function to assign partition IDs per observation----
partidx <- function(ndat, sampID = NULL, K, seed = NULL) {
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
  
  partvec <- rep(x = 1:K, length.out = ndat)
  
  # randomly assign observations to one of the K parts such that there are nK 
  # observations per part
  part.id <- sample(x = partvec, size = ndat, replace = F)
  
  return(part.id)
  
}

# partidx(ndat = 103, sampID = 2, K = 10)

#----

# function to partition data into K parts----
partition <- function(partid, dat, K) { # input here is the calibration data
  
  final <- vector("list", K)
  
  for (k in seq_len(K)) {
    test  <- dat[partid == k, ]
    train <- dat[partid != k, ]
    
    final[[k]] <- list(train = train, test = test)
  }
  
  # names(final) <- paste0("part", seq_len(K)) #FIXME maybe not needed? never used, i think
  
  return(final)
}

# part.ids <- partidx(ndat = 250, sampID = 2, K = 10)
# partition(partid = part.ids, dat = calibration, K = 10) # test function

#----

