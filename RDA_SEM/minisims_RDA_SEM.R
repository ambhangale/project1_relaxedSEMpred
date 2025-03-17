## Aditi M. Bhangale
## Last updated: 11 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## SEM on political democracy dataset
### mini simulations to test the method

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(here)
source(here("RDA_SEM", "testfunc_RDA_SEM.R"))
source(here("RDA_SEM", "customfunc_RDA_SEM.R"))

# experimenting only with formula 3b because 3a(i) and (ii) are special cases of this formula anyways

#TODO insert preliminary check?

# WHEN MODEL IS CORRECTLY SPECIFIED

conds.cs <- expand.grid(std.data = T, misspecify = F, 
                      regXY = T, alpha1 = seq(0,1,0.1), 
                      alpha2 = seq(0,1,0.1), stringsAsFactors = F) # correctly-specified conds

## mini simulation with small train set (100) and large test set (10000), correctly specified model----
t0.cs.1e2 <- Sys.time()
resList.cs.1e2 <- mapply(testrule, ntrain = 100, ntest = 1e4, 
                         std.data = conds.cs$std.data, misspecify = conds.cs$misspecify, 
                         regXY = conds.cs$regXY, alpha1 = conds.cs$alpha1, 
                         alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e2 <- Sys.time()
diff.cs.1e2 <- difftime(t1.cs.1e2, t0.cs.1e2, "sec")

RMSE.cs.1e2 <- sum_result(resList.cs.1e2)
min.RMSEp.cs.1e2 <- subset(RMSE.cs.1e2$RMSEp, RMSE.cs.1e2$RMSEp$RMSEp == min(RMSE.cs.1e2$RMSEp$RMSEp)) # TODO custom function

RMSEp.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for SMALL train set, correctly specified model") +
  geom_point(data = min.RMSEp.cs.1e2, colour = "red") # TODO custom function

RMSEpr.plot.cs.1e2 <- plot_result(RMSE.cs.1e2$RMSEpr, plot.stat = "RMSEpr",
                                 plot.title = "RMSEpr for SMALL train set, correctly specified model")

##----

## mini simulation with large train set (10000) and large test set (10000), correctly specified model----
t0.cs.1e4 <- Sys.time()
resList.cs.1e4 <- mapply(testrule, ntrain = 1e4, ntest = 1e4, 
                         std.data = conds.cs$std.data, misspecify = conds.cs$misspecify, 
                         regXY = conds.cs$regXY, alpha1 = conds.cs$alpha1, 
                         alpha2 = conds.cs$alpha2, SIMPLIFY = F)
t1.cs.1e4 <- Sys.time()
diff.cs.1e4 <- difftime(t1.cs.1e4, t0.cs.1e4, "sec")

RMSE.cs.1e4 <- sum_result(resList.cs.1e4)
min.RMSEp.cs.1e4 <- subset(RMSE.cs.1e4$RMSEp, RMSE.cs.1e4$RMSEp$RMSEp == min(RMSE.cs.1e4$RMSEp$RMSEp)) # TODO custom function

RMSEp.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEp, plot.stat = "RMSEp",
                                 plot.title = "RMSEp for LARGE train set, correctly specified model") +
  geom_point(data = min.RMSEp.cs.1e4, colour = "red") # TODO custom function

RMSEpr.plot.cs.1e4 <- plot_result(RMSE.cs.1e4$RMSEpr, plot.stat = "RMSEpr",
                                  plot.title = "RMSEpr for LARGE train set, correctly specified model")
##----

# WHEN MODEL IS MISSPECIFIED

conds.ms <- rbind(expand.grid(std.data = T, misspecify = T,
                              cor.strength = NA, reg.strength = NA, 
                              regXY = T, alpha1 = seq(0,1,0.1), 
                              alpha2 = seq(0,1,0.1), stringsAsFactors = F),
                  expand.grid(std.data = T, misspecify = T,
                              cor.strength = -0.3, reg.strength = 0.4, 
                              regXY = T, alpha1 = seq(0,1,0.1), 
                              alpha2 = seq(0,1,0.1), stringsAsFactors = F),
                  expand.grid(std.data = T, misspecify = T,
                              cor.strength = 0.1, reg.strength = 0.9, 
                              regXY = T, alpha1 = seq(0,1,0.1), 
                              alpha2 = seq(0,1,0.1), stringsAsFactors = F))

## mini simulation with small train set (100) and large test set (10000), misspecified model----
t0.ms.1e2 <- Sys.time()
resList.ms.1e2 <- mapply(testrule, ntrain = 100, ntest = 1e2, 
                         std.data = conds.ms$std.data, misspecify = conds.ms$misspecify,
                         cor.strength = conds.ms$cor.strength, reg.strength = conds.ms$reg.strength,
                         regXY = conds.ms$regXY, alpha1 = conds.ms$alpha1, 
                         alpha2 = conds.ms$alpha2, SIMPLIFY = F)
t1.ms.1e2 <- Sys.time()
diff.ms.1e2 <- difftime(t1.ms.1e2, t0.ms.1e2, "sec")

RMSE.ms.1e2 <- sum_result(resList.ms.1e2)
RMSE.ms.1e2$RMSEp$ms.strength <- ifelse(is.na(RMSE.ms.1e2$RMSEp$cor.strength) & 
                                          is.na(RMSE.ms.1e2$RMSEp$reg.strength), "estimated",
                                        ifelse(RMSE.ms.1e2$RMSEp$cor.strength == -0.3 & 
                                                 RMSE.ms.1e2$RMSEp$reg.strength == 0.4, "smaller",
                                               "larger")) # just for now, create labels in columns for RMSEp
min.RMSEp.ms.1e2 <- do.call("rbind", lapply(unique(RMSE.ms.1e2$RMSEp$ms.strength), 
                                            function(x) subset(RMSE.ms.1e2$RMSEp[RMSE.ms.1e2$RMSEp$ms.strength == x,], 
                                                               RMSE.ms.1e2$RMSEp[RMSE.ms.1e2$RMSEp$ms.strength == x,]$RMSEp == 
                                                                 min(RMSE.ms.1e2$RMSEp[RMSE.ms.1e2$RMSEp$ms.strength == x,]$RMSEp)))) # TODO custom function

ggplot(data = RMSE.ms.1e2$RMSEp, mapping = aes(x = alpha1, y = RMSEp)) +
  geom_point() + 
  facet_wrap(~ ms.strength + alpha2, ncol = 11) + geom_point(data = min.RMSEp.ms.1e2, colour = "red") + # TODO custom function
  ggtitle("RMSEp for SMALL train set, misspecified model")
## TODO eventually add to custom function. also below.

##----

## mini simulation with large train set (10000) and large test set (10000), misspecified model----
t0.ms.1e4 <- Sys.time()
resList.ms.1e4 <- mapply(testrule, ntrain = 1e4, ntest = 1e4, 
                         std.data = conds.ms$std.data, misspecify = conds.ms$misspecify,
                         cor.strength = conds.ms$cor.strength, reg.strength = conds.ms$reg.strength,
                         regXY = conds.ms$regXY, alpha1 = conds.ms$alpha1, 
                         alpha2 = conds.ms$alpha2, SIMPLIFY = F)
t1.ms.1e4 <- Sys.time()
diff.ms.1e4 <- difftime(t1.ms.1e4, t0.ms.1e4, "sec")

RMSE.ms.1e4 <- sum_result(resList.ms.1e4)
RMSE.ms.1e4$RMSEp$ms.strength <- ifelse(is.na(RMSE.ms.1e4$RMSEp$cor.strength) & 
                                          is.na(RMSE.ms.1e4$RMSEp$reg.strength), "estimated",
                                        ifelse(RMSE.ms.1e4$RMSEp$cor.strength == -0.3 & 
                                                 RMSE.ms.1e4$RMSEp$reg.strength == 0.4, "smaller",
                                               "larger")) # just for now, create labels in columns for RMSEp
min.RMSEp.ms.1e4 <- do.call("rbind", lapply(unique(RMSE.ms.1e4$RMSEp$ms.strength), 
                                            function(x) subset(RMSE.ms.1e4$RMSEp[RMSE.ms.1e4$RMSEp$ms.strength == x,], 
                                                               RMSE.ms.1e4$RMSEp[RMSE.ms.1e4$RMSEp$ms.strength == x,]$RMSEp == 
                                                                 min(RMSE.ms.1e4$RMSEp[RMSE.ms.1e4$RMSEp$ms.strength == x,]$RMSEp)))) # TODO custom function

ggplot(data = RMSE.ms.1e4$RMSEp, mapping = aes(x = alpha1, y = RMSEp)) +
  geom_point() + 
  facet_wrap(~ ms.strength + alpha2, ncol = 11) + geom_point(data = min.RMSEp.ms.1e4, colour = "red") + # TODO custom function
  ggtitle("RMSEp for LARGE train set, misspecified model")
## TODO eventually add to custom function. also below.

##----


