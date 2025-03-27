## Aditi M. Bhangale
## Last updated: 24 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### testing the results from the function(s) file (match to `lavPredictY()` and `lmPred()`)

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "func_RDA_CFA.R"))
library(ggplot2)

dat <- gendat(nCal = 250, nPred = 250, misspecify = F) # data to compare results with `lavaan` and `lm()`

# compare with `lavPredictY()` to see if `alpha1/2 = 0` is equivalent to SEM prediction rule----
fit <- fitmod(dat$calibration)

Ypred.lav <- lavPredictY(fit, newdata = dat$prediction, ynames = paste0("x", 1:3), 
                         xnames = paste0("x", 4:7))
Ypred.A0 <- testrule(nCal = 250, nPred = 250, misspecify = F, alpha1 = 0, alpha2 = 0)$Ypred

stopifnot(all(round(Ypred.lav,9) == round(Ypred.A0,9))) 
## everything matches, `testrule` with `alpha1 = alpha2 = 0` is equivalent to the SEM prediction rule

#----

# compare with OLS regression to see if `alpha1/2 = 1` is equivalent to OLS prediction----
lmPred <- function(calibration, prediction, 
                   xnames = paste0("x", 4:7), 
                   ynames = paste0("x", 1:3)) {
  Ypred <- list()
  for (y in ynames) {
    lmfit <- lm(formula(paste0(y, "~", paste0(xnames, collapse = "+"))), data = calibration)
    Ypred[[match(y, ynames)]] <- predict.lm(lmfit, newdata = prediction)
  }
  Ypred.df <- do.call("cbind", Ypred)
  colnames(Ypred.df) <- ynames
  return(Ypred.df)
}

Ypred.lm <- lmPred(dat$calibration, dat$prediction)
Ypred.A1 <- testrule(nCal = 250, nPred = 250, misspecify = F, alpha1 = 1, alpha2 = 1)$Ypred

stopifnot(all(round(Ypred.lm,9) == round(Ypred.A1,9)))
## everything matches, `testrule` with `alpha1 = alpha2 = 1` is equivalent to the OLS prediction rule

#----

