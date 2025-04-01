## Aditi M. Bhangale
## Last updated: 1 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### mini simulation testing prediction rule in cross-validation setting

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "cv_RDA_CFA.R"))
source(here("RDA_SEM", "customfunc_RDA_SEM.R"))

# CL <- makePSOCKcluster(detectCores() - 1L)

conds <- expand.grid(sampID = 1:100, nCal = c(250, 1e3, 1e4), 
                     misspecify = c(F, T))

# t0 <- Sys.time()
# resList <- mcmapply(FUN = predict.y.cv, sampID = conds$sampID,
#                     nCal = conds$nCal, nPred = 1e4,
#                     misspecify = conds$misspecify, SIMPLIFY = F,
#                     mc.cores = detectCores() - 1L) # parallelised version of mapply (parallel package)
# t1 <- Sys.time()
# diff <- difftime(t1, t0, "hour")
# saveRDS(resList, file = paste0("resList_", Sys.Date(),".rds"))

# resList <- mapply(FUN = predict.y.cv, sampID = conds$sampID,
#                   nCal = conds$nCal, nPred = 1e4,
#                   misspecify = conds$misspecify, SIMPLIFY = F) # not parallelised

resList <- readRDS("resList_2025-04-01.rds")

RMSEp <- sum_result(resList = resList)$RMSEp
RMSEp$alphas <- paste0(RMSEp$alpha1, ",", RMSEp$alpha2)

tab <- table(RMSEp$alphas, RMSEp$nCal, RMSEp$misspecify)

## looking at SEM prediction rule, OLS prediction rule, perfect compromise, and two in-between scenarios
tab[,,"0"][c("0,0", "0,1", "0.5,0.5", "1,0", "1,1"),] ## for correctly specified model
tab[,,"1"][c("0,0", "0,1", "0.5,0.5", "1,0", "1,1"),] ## for misspecified model

