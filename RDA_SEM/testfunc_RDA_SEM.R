## Aditi M. Bhangale
## Last updated: 24 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## testing the results from the function(s) file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

source("func_RDA_SEM.R")
library(ggplot2)

dat <- gendat(ntrain = 250, ntest = 250, misspecify = F) # for comparison with `lavaan()` and `lm()`

## compare with `lavPredictY()` function to see if `alpha1/2 = 0` with `XYtype = "Sigma.xy"` is equivalent----

### `lavPredictY()`

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

round(Ypred.lav,9) == round(Ypred.A1$Ypred,9)
round(Ypred.lav,9) == round(Ypred.A2$Ypred,9)
round(Ypred.A1$Ypred,9) == round(Ypred.A2$Ypred,9)
# the three are comparable
# rounding only to compare, because R returns FALSE due to rounding error at some nth term 

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

Ypred.A3 <- testrule(ntrain = 250, ntest = 250, misspecify = F, regXY = F, 
                     XYtype = "S.xy", alpha1 = 1)
Ypred.A4 <- testrule(ntrain = 250, ntest = 250, misspecify = F, regXY = T, 
                     alpha1 = 1, alpha = 1)

round(Ypred.lm,10) == round(Ypred.A3$Ypred,10)
round(Ypred.lm,10) == round(Ypred.A4$Ypred,10)
round(Ypred.A3$Ypred,10) == round(Ypred.A4$Ypred,10)
# the three are comparable

#----

## preliminary check of results from function----
conds1 <- rbind(expand.grid(regXY = FALSE, XYtype = c("S.xy", "Sigma.xy"), 
                            alpha1 = seq(0,1,0.1), alpha2 = NA, 
                            misspecify = c(FALSE,TRUE), stringAsFactors = F),
                expand.grid(regXY = TRUE, XYtype = NA, alpha1 = seq(0,1,0.1), 
                            alpha2 = seq(0,1,0.1), 
                            misspecify = c(FALSE,TRUE), stringAsFactors = F))
# dim(conds) should be (11*2*2) + (11*11*2) = 286

# results
t0 <- Sys.time()
resList1 <- mapply(testrule, ntrain = 250, ntest = 250, 
                  misspecify = conds1$misspecify, regXY = conds1$regXY, 
                  XYtype = as.character(conds1$XYtype), alpha1 = conds1$alpha1, 
                  alpha2 = conds1$alpha2, SIMPLIFY = F)
t1 <- Sys.time()
diff <- difftime(t1, t0, units = "sec")

RMSEpr <- as.data.frame(do.call("rbind", 
                                lapply(1:length(resList1), function(x) resList1[[x]]$RMSEpr.result)))
RMSEpr$meanBias <- as.numeric(RMSEpr$meanBias)
RMSEpr$RMSEpr <- as.numeric(RMSEpr$RMSEpr)
  
RMSEp <- as.data.frame(do.call("rbind", 
                               lapply(1:length(resList1), function(x) resList1[[x]]$RMSEp.result)))

RMSEp$RMSEp <- as.numeric(RMSEp$RMSEp)
# head(RMSEpr)
# head(RMSEp) # check

##----

