## Aditi M. Bhangale
## Last updated: 11 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## mini simulations to test the method

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(here)
source(here("RDA_SEM", "testfunc_RDA_SEM.R"))
source(here("RDA_SEM", "customfunc_RDA_SEM.R"))

# experimenting only with formula 3b because 3a(i) and (ii) are special cases of this formula anyways

#TODO insert preliminary check?

# WHEN MODEL IS CORRECTLY SPECIFIED

conds.cs <- expand.grid(std.data = T, misspecify = F, 
                      regXY = TRUE, XYtype = NA, alpha1 = seq(0,1,0.1), 
                            alpha2 = seq(0,1,0.1), stringsAsFactors = F) # correctly-specified conds

## mini simulation with small train set (100) and large test set (10000)----
t0.cs.1e2 <- Sys.time()
resList.cs.1e2 <- mapply(testrule, ntrain = 100, ntest = 1e4, 
                         std.data = conds.cs$std.data, misspecify = conds.cs$misspecify, 
                         regXY = conds.cs$regXY, XYtype = conds.cs$XYtype,
                         alpha1 = conds.cs$alpha1, alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e2 <- Sys.time()
diff.cs.1e2 <- difftime(t1.cs.1e2, t0.cs.1e2, "sec")

RMSE.cs.1e2 <- sum_result(resList.cs.1e2)

RMSEp.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for SMALL train set, correctly specified model")

RMSEpr.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEpr, plot.stat = "RMSEpr",
                                 plot.title = "RMSEpr for SMALL train set, correctly specified model")

##----

## mini simulation with large train set (10000) and large test set (10000)----
t0.cs.1e4 <- Sys.time()
resList.cs.1e4 <- mapply(testrule, ntrain = 1e4, ntest = 1e4, 
                         std.data = conds.cs$std.data, misspecify = conds.cs$misspecify, 
                         regXY = conds.cs$regXY, XYtype = conds.cs$XYtype,
                         alpha1 = conds.cs$alpha1, alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e4 <- Sys.time()
diff.cs.1e4 <- difftime(t1.cs.1e4, t0.cs.1e4, "sec")

RMSE.cs.1e4 <- sum_result(resList.cs.1e4)

RMSEp.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for LARGE train set, correctly specified model")

RMSEpr.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for LARGE train set, correctly specified model")
##----

# WHEN MODEL IS MISSPECIFIED

# TODO conds.ms <- 

## TODO mini simulation with strength of misspecification manipulated?
# estimate misspecified effects
# fix misspecified effects to a particular value


