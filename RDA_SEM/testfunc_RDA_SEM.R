## Aditi M. Bhangale
## Last updated: 21 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## testing the results from the function(s) file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

source("func_RDA_SEM.R")

## compare with `lavPredictY()` function to see if `alpha1/2 = 0` with `XYtype = "Sigma.xy"` is equivalent----

### `lavPredictY()`
dat <- gendat(ntrain = 250, ntest = 250, misspecify = F)

mod <- ' 
    # latent variable definitions
      ind60 =~ x1 + x2 + x3
      dem60 =~ y1 + y2 + y3 + y4 
      dem65 =~ y5 + y6 + y7 + y8
      
    # regressions
      dem60 ~ ind60
      dem65 ~ ind60 + dem60
      
    # residual correlations
      y1 ~~ y5
      y2 ~~ y4 + y6
      y3 ~~ y7
      y4 ~~ y8
      y6 ~~ y8
  '

fit <- sem(mod, data = dat$train, meanstructure = T)
Ypred.lav <- lavPredictY(fit, newdata = dat$test, ynames = paste0("y", 5:8),
                         xnames = c(paste0("x", 1:3), paste0("y", 1:4)))

### my function
Ypred.A1 <- testrule(ntrain = 250, ntest = 250, misspecify = F, regXY = F, 
                    XYtype = "Sigma.xy", alpha1 = 0)
Ypred.A2 <- testrule(ntrain = 250, ntest = 250, misspecify = F, regXY = T, 
                     alpha1 = 0, alpha2 = 0)

isFALSE(any(round(Ypred.lav,100) == round(Ypred.A1$Ypred,100)))
isFALSE(any(round(Ypred.lav,100) == round(Ypred.A2$Ypred,100)))
isFALSE(any(round(Ypred.A1$Ypred,100) == round(Ypred.A2$Ypred,100)))
# the three are comparable
# rounding only to compare, because R returns FALSE due to rounding error at some nth term 

##----

