## Aditi M. Bhangale
## Last updated: 15 July 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### function to predict using elastic net regression

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

library(glmnet) # to perform elastic net regression
source("part_relaxSEM.R")

# prediction with an elastic net regresison model with cross-validation----
en.predict.y.cv <- function(calidat, preddat, alphas, partid, n_y, xnames, ynames) {

  cv.errors <- rep(NA, length(alphas)) 
  
  ## in the `glmnet()` and `cv.glmnet()` calls below, i explicitly specify
  ## `standardize = T, standardize.response = F` even though they are the defaults
  ## just in case the defaults change in the future and so i can be sure of the 
  ## procedure in my simulation
  for (a in 1:length(alphas)) { # returns the optimal elastic net mixing parameter
    # save the minimum cross-validated error for each value in `alphas`
    cv.errors[a] <- min(cv.glmnet(x = as.matrix(calidat[, xnames]),
                                  y = as.matrix(calidat[, ynames]), 
                                  foldid = partid,
                                  family = ifelse(n_y == 1L, "gaussian", "mgaussian"), 
                                  alpha = alphas[a], 
                                  standardize = T, standardize.response = F)$cvm) 
    # above, use default lambda sequence because that's what Mark did in the original simulation 
  }
  
  min.alpha <- alphas[which.min(cv.errors)] # alpha value with minimum cross-validated error
  
  lambda <- cv.glmnet(x = as.matrix(calidat[, xnames]),
                      y = as.matrix(calidat[, ynames]), 
                      foldid = partid,
                      family = ifelse(n_y == 1L, "gaussian", "mgaussian"), 
                      alpha = min.alpha, 
                      standardize = T, standardize.response = F)$lambda.min # minimum/optimal tuning parameter (lambda) value
  
  out <- glmnet(x = as.matrix(calidat[, xnames]),
                y = as.matrix(calidat[, ynames]),
                family = ifelse(n_y == 1L, "gaussian", "mgaussian"), 
                alpha = min.alpha, 
                standardize = T, standardize.response = F) # final model to use for predicting new values
  
  Ypred <- if (n_y == 1L) {
    as.matrix(predict(out, newx = as.matrix(preddat[,xnames]), s = lambda)[,1])
  } else {
    as.matrix(predict(out, newx = as.matrix(preddat[,xnames]), s = lambda)[,,1])
  }
  
  
  attr(Ypred, "alpha")  <- min.alpha # en mixing parameter
  attr(Ypred, "lambda") <- lambda # tuning parameter
  
  return(Ypred)
}

# en.predict.y.cv(calidat = calibration, preddat = prediction, alphas = seq(0,1,0.1),
#                 partid = part.ids, xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----
