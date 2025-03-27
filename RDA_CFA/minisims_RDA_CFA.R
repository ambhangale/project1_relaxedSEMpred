## Aditi M. Bhangale
## Last updated: 25 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### mini simulations to prediction the method

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "testfunc_RDA_CFA.R"))
source(here("RDA_SEM", "customfunc_RDA_SEM.R"))

### CORRECTLY SPECIFIED MODEL ###
conds.cs <- expand.grid(misspecify = F, alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
                        stringsAsFactors = F)

# mini simulation for cs model with small calibration set (nCal = 100) and large prediction set (nPred = 10000)----
t0.cs.1e2 <- Sys.time()
resList.cs.1e2 <- mapply(testrule, nCal = 100, nPred = 1e4, 
                         misspecify = conds.cs$misspecify, 
                         alpha1 = conds.cs$alpha1, 
                         alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e2 <- Sys.time()
diff.cs.1e2 <- difftime(t1.cs.1e2, t0.cs.1e2, "sec")

RMSE.cs.1e2 <- sum_result(resList.cs.1e2)

RMSEp.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEp, plot.stat = "RMSEp",
                                plot.title = "RMSEp for SMALL calibration set, correctly specified model")

RMSEpr.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for SMALL calibration set, correctly specified model",
                                  ynames = paste0("x", 1:3))

#----

## mini simulation for cs model with large calibration set (10000) and large prediction set (10000)----
t0.cs.1e4 <- Sys.time()
resList.cs.1e4 <- mapply(testrule, nCal = 1e4, nPred = 1e4, 
                         misspecify = conds.cs$misspecify, 
                         alpha1 = conds.cs$alpha1, 
                         alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e4 <- Sys.time()
diff.cs.1e4 <- difftime(t1.cs.1e4, t0.cs.1e4, "sec")

RMSE.cs.1e4 <- sum_result(resList.cs.1e4)

RMSEp.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for LARGE calibration set, correctly specified model")

RMSEpr.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for LARGE calibration set, correctly specified model",
                                  ynames = paste0("x", 1:3))
##----

### MISSPECIFIED MODEL ###
conds.ms <- expand.grid(misspecify = T, alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
                        stringsAsFactors = F)

# mini simulation for ms model with small calibration set (nCal = 100) and large prediction set (nPred = 10000)----
t0.ms.1e2 <- Sys.time()
resList.ms.1e2 <- mapply(testrule, nCal = 100, nPred = 1e4, 
                         misspecify = conds.ms$misspecify, 
                         alpha1 = conds.ms$alpha1, 
                         alpha2 = conds.ms$alpha2, SIMPLIFY = F)
t1.ms.1e2 <- Sys.time()
diff.ms.1e2 <- difftime(t1.ms.1e2, t0.ms.1e2, "sec")

RMSE.ms.1e2 <- sum_result(resList.ms.1e2) # none of the lv cov matrices is positive definite

RMSEp.plot.ms.1e2 <- plot_result(RMSE.ms.1e2$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for SMALL calibration set, misspecified model")

RMSEpr.plot.ms.1e2 <- plot_result(RMSE.ms.1e2$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for SMALL calibration set, misspecified model",
                                  ynames = paste0("x", 1:3))

#----

# mini simulation for ms model with large calibration set (nCal = 10000) and large prediction set (nPred = 10000)----
t0.ms.1e4 <- Sys.time()
resList.ms.1e4 <- mapply(testrule, nCal = 1e4, nPred = 1e4, 
                         misspecify = conds.ms$misspecify, 
                         alpha1 = conds.ms$alpha1, 
                         alpha2 = conds.ms$alpha2, SIMPLIFY = F)
t1.ms.1e4 <- Sys.time()
diff.ms.1e4 <- difftime(t1.ms.1e4, t0.ms.1e4, "sec")

RMSE.ms.1e4 <- sum_result(resList.ms.1e4) # all lv cov matrices are positive definite

RMSEp.plot.ms.1e4 <- plot_result(RMSE.ms.1e4$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for LARGE calibration set, misspecified model")

RMSEpr.plot.ms.1e4 <- plot_result(RMSE.ms.1e4$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for LARGE calibration set, misspecified model",
                                  ynames = paste0("x", 1:3))

#----

# mini simulation for ms model with nCal = 250 and nPred = 10000----
t0.ms.250 <- Sys.time()
resList.ms.250 <- mapply(testrule, nCal = 250, nPred = 250, 
                         misspecify = conds.ms$misspecify, 
                         alpha1 = conds.ms$alpha1, 
                         alpha2 = conds.ms$alpha2, SIMPLIFY = F)
t1.ms.250 <- Sys.time()
diff.ms.250 <- difftime(t1.ms.250, t0.ms.250, "sec")

RMSE.ms.250 <- sum_result(resList.ms.250) # all lv cov matrices are positive definite

RMSEp.plot.ms.250 <- plot_result(RMSE.ms.250$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for nCal = 250, misspecified model")

RMSEpr.plot.ms.250 <- plot_result(RMSE.ms.250$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for nCal = 250, misspecified model",
                                  ynames = paste0("x", 1:3))
#----

