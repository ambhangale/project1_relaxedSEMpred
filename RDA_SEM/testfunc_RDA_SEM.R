## Aditi M. Bhangale
## Last updated: 11 March 2025
## Fixes and improvements: 3 March 2025 (Julian D. Karch)

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## SEM on political democracy dataset
### testing the results from the function(s) file (match to `lavPredictY()` and `lmPred()`)

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(here)
source(here("RDA_SEM", "func_RDA_SEM.R"))
library(ggplot2)

dat <- gendat(ntrain = 250, ntest = 250, std.data = T, misspecify = F) # for comparison with `lavaan()` and `lm()`

## compare with `lavPredictY()` function to see if `alpha1/2 = 0` with `XYtype = "Sigma.xy"` is equivalent----

### `lavPredictY()`

fit <- fitmod(dat$train)
Ypred.lav <- lavPredictY(fit, newdata = dat$test, ynames = paste0("y", 5:8),
                         xnames = c(paste0("x", 1:3), paste0("y", 1:4)))

### my function
Ypred.A1 <- testrule(ntrain = 250, ntest = 250, std.data = T, misspecify = F, regXY = F, 
                     XYtype = "Sigma.xy", alpha1 = 0)
Ypred.A2 <- testrule(ntrain = 250, ntest = 250, std.data = T, misspecify = F, regXY = T, 
                     alpha1 = 0, alpha2 = 0)

stopifnot(all(round(Ypred.lav,9) == round(Ypred.A1$Ypred,9)))
stopifnot(all(round(Ypred.lav,9) == round(Ypred.A2$Ypred,9)))
stopifnot(all(round(Ypred.A1$Ypred,9) == round(Ypred.A2$Ypred,9)))
# the three are comparable
# also works when `std.data = F` in `gendata()` and `testrule()`
# rounding only to compare, because R returns FALSE due to rounding error at some nth term. also below.

##----

## compare with OLS regression to see if `alpha1/2 = 1` with `XYtype = "S.xy"` is equivalent----

lmPred <- function(train, test, 
                   xnames = c(paste0("x", 1:3), paste0("y", 1:4)), 
                   ynames = paste0("y", 5:8)) {
  Ypred <- list()
  for (y in ynames) {
    lmfit <- lm(formula(paste0(y, "~", paste0(xnames, collapse = "+"))), data = train)
    Ypred[[match(y, ynames)]] <- predict.lm(lmfit, newdata = test)
  }
  Ypred.df <- do.call("cbind", Ypred)
  colnames(Ypred.df) <- ynames
  return(Ypred.df)
}

Ypred.lm <- lmPred(train = dat$train, test = dat$test)

Ypred.A3 <- testrule(ntrain = 250, ntest = 250, std.data = T, misspecify = F, regXY = F, 
                     XYtype = "S.xy", alpha1 = 1)
Ypred.A4 <- testrule(ntrain = 250, ntest = 250, std.data = T, misspecify = F, regXY = T, 
                     alpha1 = 1, alpha = 1)

stopifnot(all(round(Ypred.lm,10) == round(Ypred.A3$Ypred,10)))
stopifnot(all(round(Ypred.lm,10) == round(Ypred.A4$Ypred,10)))
stopifnot(all(round(Ypred.A3$Ypred,10) == round(Ypred.A4$Ypred,10)))
# the three are comparable
# also works when `std.data = F` in `gendata()` and `testrule()`

#----


