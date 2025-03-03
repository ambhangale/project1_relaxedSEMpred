## Aditi M. Bhangale
## Last updated: 27 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## testing the results from the function(s) file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

source("func_RDA_SEM.R")
library(ggplot2)

dat <- gendat(ntrain = 250, ntest = 250, misspecify = F) # for comparison with `lavaan()` and `lm()`

## compare with `lavPredictY()` function to see if `alpha1/2 = 0` with `XYtype = "Sigma.xy"` is equivalent----

### `lavPredictY()`

fit <- fitmod(dat$train)
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
diff1 <- difftime(t1, t0, units = "sec")

RMSEpr1 <- as.data.frame(do.call("rbind", 
                                lapply(1:length(resList1), function(x) resList1[[x]]$RMSEpr.result)))
RMSEpr1$meanBias <- as.numeric(RMSEpr1$meanBias)
RMSEpr1$RMSEpr <- as.numeric(RMSEpr1$RMSEpr)
  
RMSEp1 <- as.data.frame(do.call("rbind", 
                               lapply(1:length(resList1), function(x) resList1[[x]]$RMSEp.result)))

RMSEp1$RMSEp <- as.numeric(RMSEp1$RMSEp)
# head(RMSEpr1)
# head(RMSEp1) # check

## plot results

### RMSEp1 (all outcomes combined)
# RMSEp1, regXY = F
ggplot(data = RMSEp1[RMSEp1$regXY == F & RMSEp1$XYtype == "S.xy" &  is.na(RMSEp1$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT1: RMSEp1 values for regXY = F and XYtype = 'S.xy'")
ggplot(data = RMSEp1[RMSEp1$regXY == F & RMSEp1$XYtype == "Sigma.xy" &  is.na(RMSEp1$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT2: RMSEp1 values for regXY = F and XYtype = 'Sigma.xy'")

# RMSEp1, regXY = T
# misspecify = F
ggplot(data = RMSEp1[RMSEp1$regXY == T & RMSEp1$misspecify == F &  !is.na(RMSEp1$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT3: RMSEp1 values for regXY = T and misspecify = F")
# misspecify = T
ggplot(data = RMSEp1[RMSEp1$regXY == T & RMSEp1$misspecify == T &  !is.na(RMSEp1$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT4: RMSEp1 values for regXY = T and misspecify = T")


### RMSEpr1
# RMSEpr1, regXY = F
ggplot(data = RMSEpr1[RMSEpr1$regXY == F & RMSEpr1$XYtype == "S.xy" &  is.na(RMSEpr1$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEpr)) + geom_point() + facet_wrap(~ yname + misspecify, ncol = 2) +
  ggtitle("PLOT5: RMSEpr1 values for regXY = F, and XYtype = 'S.xy'")
ggplot(data = RMSEpr1[RMSEpr1$regXY == F & RMSEpr1$XYtype == "Sigma.xy" &  is.na(RMSEpr1$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEpr)) + geom_point() + facet_wrap(~ yname + misspecify, ncol = 2) +
  ggtitle("PLOT6: RMSEpr1 values for regXY = F, and XYtype = 'Sigma.xy'")

# RMSEpr1, regXY = T
# misspecify = F
ggplot(data = RMSEpr1[RMSEpr1$regXY == T & RMSEpr1$misspecify == F &  !is.na(RMSEpr1$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEpr)) + geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) +
  ggtitle("PLOT7: RMSEpr1 values for regXY = T and misspecify = F")
# misspecify = T
ggplot(data = RMSEpr1[RMSEpr1$regXY == T & RMSEpr1$misspecify == T &  !is.na(RMSEpr1$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEpr)) + geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) +
  ggtitle("PLOT8: RMSEpr1 values for regXY = T and misspecify = T")

# ggplot(data = RMSEpr1[RMSEpr1$regXY == T & RMSEpr1$misspecify == F &  !is.na(RMSEpr1$alpha2) & RMSEpr1$yname == "y5",], 
#        mapping = aes(x = alpha2, y = RMSEpr1)) + geom_point() + facet_wrap(~ alpha1)

# FIXME eventually create a custom function to make all these plots (maybe even for extracting RMSEp1(r))

##----

## experimenting with large training set (1000) and small test set (100)----

# use conds1

# results
t0 <- Sys.time()
resList2 <- mapply(testrule, ntrain = 1000, ntest = 100, 
                   misspecify = conds1$misspecify, regXY = conds1$regXY, 
                   XYtype = as.character(conds1$XYtype), alpha1 = conds1$alpha1, 
                   alpha2 = conds1$alpha2, SIMPLIFY = F)
t1 <- Sys.time()
diff2 <- difftime(t1, t0, units = "sec")

RMSEpr2 <- as.data.frame(do.call("rbind", 
                                 lapply(1:length(resList2), function(x) resList2[[x]]$RMSEpr.result)))
RMSEpr2$meanBias <- as.numeric(RMSEpr2$meanBias)
RMSEpr2$RMSEpr <- as.numeric(RMSEpr2$RMSEpr)

RMSEp2 <- as.data.frame(do.call("rbind", 
                                lapply(1:length(resList2), function(x) resList2[[x]]$RMSEp.result)))

RMSEp2$RMSEp <- as.numeric(RMSEp2$RMSEp)

## plot results

### RMSEp2 (all outcomes combined)
# RMSEp2, regXY = F
ggplot(data = RMSEp2[RMSEp2$regXY == F & RMSEp2$XYtype == "S.xy" &  is.na(RMSEp2$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT1: RMSEp2 values for regXY = F and XYtype = 'S.xy'")
ggplot(data = RMSEp2[RMSEp2$regXY == F & RMSEp2$XYtype == "Sigma.xy" &  is.na(RMSEp2$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT2: RMSEp2 values for regXY = F and XYtype = 'Sigma.xy'")

# RMSEp2, regXY = T
# misspecify = F
ggplot(data = RMSEp2[RMSEp2$regXY == T & RMSEp2$misspecify == F &  !is.na(RMSEp2$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT3: RMSEp2 values for regXY = T and misspecify = F")
# misspecify = T
ggplot(data = RMSEp2[RMSEp2$regXY == T & RMSEp2$misspecify == T &  !is.na(RMSEp2$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT4: RMSEp2 values for regXY = T and misspecify = T")

#----

## experimenting with small training set (100) and small test set (100)----

# use conds1

# results
t0 <- Sys.time()
resList3 <- mapply(testrule, ntrain = 100, ntest = 100, 
                   misspecify = conds1$misspecify, regXY = conds1$regXY, 
                   XYtype = as.character(conds1$XYtype), alpha1 = conds1$alpha1, 
                   alpha2 = conds1$alpha2, SIMPLIFY = F)
t1 <- Sys.time()
diff3 <- difftime(t1, t0, units = "sec")

RMSEpr3 <- as.data.frame(do.call("rbind", 
                                 lapply(1:length(resList3), function(x) resList3[[x]]$RMSEpr.result)))
RMSEpr3$meanBias <- as.numeric(RMSEpr3$meanBias)
RMSEpr3$RMSEpr <- as.numeric(RMSEpr3$RMSEpr)

RMSEp3 <- as.data.frame(do.call("rbind", 
                                lapply(1:length(resList3), function(x) resList3[[x]]$RMSEp.result)))

RMSEp3$RMSEp <- as.numeric(RMSEp3$RMSEp)

## plot results

### RMSEp3 (all outcomes combined)
# RMSEp3, regXY = F
ggplot(data = RMSEp3[RMSEp3$regXY == F & RMSEp3$XYtype == "S.xy" &  is.na(RMSEp3$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT1: RMSEp3 values for regXY = F and XYtype = 'S.xy'")
ggplot(data = RMSEp3[RMSEp3$regXY == F & RMSEp3$XYtype == "Sigma.xy" &  is.na(RMSEp3$alpha2),], 
       mapping = aes(x = alpha1, y = RMSEp)) + geom_point() + facet_wrap(~ misspecify) +
  ggtitle("PLOT2: RMSEp3 values for regXY = F and XYtype = 'Sigma.xy'")

# RMSEp3, regXY = T
# misspecify = F
ggplot(data = RMSEp3[RMSEp3$regXY == T & RMSEp3$misspecify == F &  !is.na(RMSEp3$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT3: RMSEp3 values for regXY = T and misspecify = F")
# misspecify = T
ggplot(data = RMSEp3[RMSEp3$regXY == T & RMSEp3$misspecify == T &  !is.na(RMSEp3$alpha2),], 
       mapping = aes(x = alpha2, y = RMSEp)) + geom_point() + facet_wrap(~ alpha1) +
  ggtitle("PLOT4: RMSEp3 values for regXY = T and misspecify = T")

##----


