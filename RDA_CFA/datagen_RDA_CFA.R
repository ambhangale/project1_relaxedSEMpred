## Aditi M. Bhangale
## Last updated: 20 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### data generation file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(mvtnorm)
library(lavaan)

mod1 <- '
# factor loadings
F1 =~ 0.67*x1 + 0.71*x2 + 0.58*x3
F2 =~ 0.83*x4 + 0.64*x5 + 0.55*x6 + 0.68*x7

# factor (co)variances
F1 ~~ 1*F1 + 0.3*F2
F2 ~~ 1*F2

# item (co)variances
x1 ~~ 0.22*x1
x2 ~~ 0.41*x2
x3 ~~ 0.55*x3
x4 ~~ 0.49*x4
x5 ~~ 0.71*x5
x6 ~~ 0.39*x6
x7 ~~ 0.61*x7

# factor means
F1 ~ 0*1
F2 ~ 0*1

# item means
x1 ~ 0.91*1
x2 ~ 0.43*1
x3 ~ 0.73*1
x4 ~ 0.44*1
x5 ~ 0.83*1
x6 ~ 0.66*1
x7 ~ 0.54*1
'

mod2 <- '
# factor loadings
F1 =~ 0.67*x1 + 0.71*x2 + 0.58*x3 + 0.83*x4 + 0.64*x5 + 0.55*x6 + 0.68*x7

# factor variance
F1 ~~ 1*F1

# item (co)variances
x1 ~~ 0.22*x1
x2 ~~ 0.41*x2
x3 ~~ 0.55*x3
x4 ~~ 0.49*x4
x5 ~~ 0.71*x5
x6 ~~ 0.39*x6
x7 ~~ 0.61*x7

# factor mean
F1 ~ 0*1

# item means
x1 ~ 0.91*1
x2 ~ 0.43*1
x3 ~ 0.73*1
x4 ~ 0.44*1
x5 ~ 0.83*1
x6 ~ 0.66*1
x7 ~ 0.54*1
'

gendat <- function(ntrain, ntest, misspecify, seed = 10824) {
  
  # (do not) fit model, return start values
  fit <- if(!misspecify) {lavaan(mod1, do.fit = F)
  } else {
    lavaan(mod2, do.fit = F)
  }
  pop.cov <- lavInspect(fit, "cov.ov") # population covariance matrix
  pop.mean <- lavInspect(fit, "mean.ov") # population mean vector
  
  # data generation
  set.seed(seed)
  train <- as.data.frame(rmvnorm(n = ntrain, mean = pop.mean, sigma = pop.cov,
                                 pre0.9_9994 = T)) # training set
  test  <- as.data.frame(rmvnorm(n = ntest, mean = pop.mean, sigma = pop.cov,
                                 pre0.9_9994 = T)) # test set
  
  datlist <- list(train = train, test = test)
  
  ## add attributes in case you need them later
  attr(datlist, "ntrain") <- ntrain
  attr(datlist, "ntest") <- ntest
  attr(datlist, "misspecify") <- misspecify
  attr(datlist, "seed") <- seed
  
  return(datlist)
}

# test function
# gendat(ntrain = 20, ntest = 25, misspecify = F)
# gendat(ntrain = 25, ntest = 20, misspecify = T)


