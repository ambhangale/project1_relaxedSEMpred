## Aditi M. Bhangale
## Last updated: 26 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## data generation file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(lavaan)
library(MASS) # for mvrnorm()

pol_model <- ' 
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

## function for data generation
gendat <- function(ntrain, ntest, misspecify, mod = pol_model, seed = 10824) {
  
  data("PoliticalDemocracy")
  
  if (misspecify) mod <- paste(mod, 'dem60 ~ x1 \n', collapse = "\n") # add additional direct effect for misspecification
  
  fit <- sem(mod, data = PoliticalDemocracy, meanstructure = T) # mean structure is saturated
  popStats <- lavInspect(fit, "implied") # use to generate new data
  
  set.seed(seed)
  train <- as.data.frame(mvrnorm(n = ntrain, mu = popStats$mean, Sigma = popStats$cov))
  test <-  as.data.frame(mvrnorm(n = ntest , mu = popStats$mean, Sigma = popStats$cov))
  
  datlist <- list(train = train, test = test)
  attr(datlist, "misspecify") <- misspecify
  return(datlist)
}

# test function
# gendat(ntrain = 15, ntest = 15, misspecify = T)
# gendat(ntrain = 15, ntest = 15, misspecify = F)


