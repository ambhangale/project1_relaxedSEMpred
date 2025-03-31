## Aditi M. Bhangale
## Last updated: 31 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### mini simulation testing prediction rule in cross-validation setting

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "cv_RDA_CFA.R"))

# CL <- makePSOCKcluster(detectCores() - 1L)

conds <- expand.grid(sampID = 1:100, nCal = c(250, 1e3, 1e4), 
                     misspecify = c(F, T))

t0 <- Sys.time()
resList <- mcmapply(FUN = predict.y.cv, sampID = conds$sampID,
                    nCal = conds$nCal, nPred = 1e4,
                    misspecify = conds$misspecify, SIMPLIFY = F, 
                    mc.cores = detectCores() - 1L) # parallelised version of mapply (parallel package)
t1 <- Sys.time()
diff <- difftime(t1, t0, "hour")

# resList <- mapply(FUN = predict.y.cv, sampID = conds$sampID,
#                   nCal = conds$nCal, nPred = 1e4,
#                   misspecify = conds$misspecify, SIMPLIFY = F) # not parallelised

