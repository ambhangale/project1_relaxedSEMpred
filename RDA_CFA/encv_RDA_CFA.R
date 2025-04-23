## Aditi M. Bhangale
## Last updated: 23 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### function to predict using elastic net regression

library(here)
library(glmnet) # to perform elastic net regression
source(here("RDA_CFA", "part_RDA_CFA.R")) 

# prediction with an elastic net regresison model with cross-validation----
en.predict.y.cv <- function(calidat, preddat, alphas, partid, xnames, ynames) {

  cv.errors <- rep(NA, length(alphas)) 
  
  for (a in 1:length(alphas)) { # returns the optimal elastic net mixing parameter
    # save the minimum cross-validated error for each value in `alphas`
    cv.errors[a] <- min(cv.glmnet(x = calidat[, xnames],
                                  y = calidat[, ynames], 
                                  foldid = partid,
                                  family = "mgaussian", alpha = alphas[a])$cvm) 
    # above, use default lambda sequence because that's what Mark did in the original simulation 
  }
  
  min.alpha <- alphas[which.min(cv.errors)] # alpha value with minimum cross-validated error
  
  lambda <- cv.glmnet(x = calidat[, xnames],
                      y = calidat[, ynames], 
                      foldid = partid,
                      family = "mgaussian", 
                      alpha = min.alpha)$lambda.min # minimum/optimal tuning parameter (lambda) value
  
  out <- glmnet(x = calidat[, xnames],
                y = calidat[, ynames],
                family = "mgaussian", alpha = min.alpha) # final model to use for predicting new values
  
  Ypred <- as.matrix(predict(out, newx = preddat[,xnames], s = lambda)[,,1])
  
  attr(Ypred, "alpha")  <- min.alpha # en mixing parameter
  attr(Ypred, "lambda") <- lambda # tuning parameter
  
  return(Ypred)
}

# en.predict.y.cv(calidat = calibration, preddat = prediction, alphas = seq(0,1,0.1),
#                 partid = part.ids, xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----
