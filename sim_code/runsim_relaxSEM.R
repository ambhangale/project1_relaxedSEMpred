## Aditi M. Bhangale
## Last updated: 30 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### runsim

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

source("wrapper_relaxSEM.R")

conds <- expand.grid(sampID = 1:100, nCal = c(100, 250, 1e3, 1e4),
                     misspecify = c(F,T))

t0 <- Sys.time()
resList <- mcmapply(FUN = wrapper.predict.y, sampID = conds$sampID,
                    nCal = conds$nCal, nPred = 1e4,
                    misspecify = conds$misspecify, SIMPLIFY = F,
                    mc.cores = detectCores() - 1L)
# parallelised version of mapply (parallel package). parallel should already be
# loaded, since it's a dependency for portableParallelSeeds used in the datagen file
t1 <- Sys.time()
diff <- difftime(t1, t0, "hour")

saveRDS(resList, file = paste0("PD_resList_itemsdirect_", Sys.Date(), ".rds"))

# resList <- readRDS()


