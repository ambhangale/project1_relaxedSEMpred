## Aditi M. Bhangale
## Last updated: 10 March 2025
## Fixes and improvements: 3 March 2025 (Julian D. Karch)

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
' # population model for political democracy data

## function for data generation
gendat <- function(ntrain, ntest, std.data, misspecify, mod = pol_model, 
                   cor.strength = NA, reg.strength = NA, seed = 10824) {
  
  data("PoliticalDemocracy")
  
  if (misspecify) mod <- paste0(c(mod, paste0('x2 ~~ ', 
                                       ifelse(is.na(cor.strength), NA, cor.strength), 
                                       '*x3'), 
                           paste0('dem60 ~ ', 
                                  ifelse(is.na(reg.strength), NA, reg.strength), 
                                  '*x1')), sep = "\n")
  # add additional direct effect and covariance/correlation for misspecification
  
  if(!misspecify & (!is.na(cor.strength) | !is.na(reg.strength))) message("`cor.strength` and `reg.strength` not used unless `misspecify = T`")
  
  if (!std.data) {
    fit <- sem(mod, data = PoliticalDemocracy, std.lv = T, meanstructure = T) 
    # mean structure is saturated, UVI constraint
    
    popStats <- lavInspect(fit, "implied") # use to generate new data
  } else {
    fit <- sem(mod, data = as.data.frame(scale(PoliticalDemocracy)), std.lv = T, meanstructure = T) 
    # above + set obs V means to 0 and variance to 1
    
    popStats <- lapply(c("cor.ov", "mean.ov"), function(x) lavInspect(fit, x)) # when dat stdised, use correlation matrix instead
    names(popStats) <- list("cov", "mean")
  }
  
  # generate training and test set
  set.seed(seed)
  train <- as.data.frame(mvrnorm(n = ntrain, mu = popStats$mean, Sigma = popStats$cov))
  test <-  as.data.frame(mvrnorm(n = ntest , mu = popStats$mean, Sigma = popStats$cov))
  
  datlist <- list(train = train, test = test)
  
  # add attributes in case required later
  attr(datlist, "misspecify") <- misspecify
  attr(datlist, "std.data") <- std.data
  attr(datlist, "cor.strength") <- cor.strength
  attr(datlist, "reg.strength") <- reg.strength
  
  return(datlist)
}

# test function
# gendat(ntrain = 15, ntest = 15, std.data = F,  misspecify = T)
# gendat(ntrain = 15, ntest = 15, std.data = T, misspecify = F)
# gendat(ntrain = 15, ntest = 15, std.data = T, reg.strength = 0.4, misspecify = T)
# gendat(ntrain = 15, ntest = 15, std.data = T, cor.strength = 0.3, reg.strength = 0.4, misspecify = T)

