## Aditi M. Bhangale
## Last updated: 23 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### wrapper function for De Rooij, OLS, lavcv, and encv rules

library(here)
source(here("RDA_CFA", "lavcv_RDA_CFA.R"))
source(here("RDA_CFA", "encv_RDA_CFA.R")) 

wrapper.predict.y <- function(sampID, nCal, nPred, misspecify, lav.CV = TRUE,
                              lav.alpha1 = seq(0,1,0.1), lav.alpha2 = seq(0,1,0.1), 
                              en.alphas = seq(0,1,0.1), K = 10, nK = NULL, 
                              xnames = paste0("x", 4:7), ynames = paste0("x", 1:3),
                              seed = NULL) {
  
  dat <- gendat(sampID = sampID, nCal = nCal, nPred = nPred, 
                misspecify = misspecify) # always require `sampID` when called in this wrapper function
  calibration <- dat$calibration # calibration set
  prediction <- dat$prediction # prediction set
  
  nK <- ifelse(is.null(nK), nCal/K, nK) # compute nK manually if left blank
  
  # partition IDs to be used for lav.predict.y.cv and en.predict.y.cv
  partIDx <- partidx(ndat = nCal, sampID = sampID, K = K, nK = nK)
  
  lav.fit <- fitmod(dat = calibration) # lavaan model fitted on complete calibration set
  
  # check if cov.lv and cov.ov are positive definite
  PD.lv <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.lv"))$values > 0), F, T)
  PD.ov <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.ov"))$values > 0), F, T)
  
  # Predictions using De Rooij et al. (2022) rule
  DeRooij.Ypred <- lavPredictY(object = lav.fit, newdata = prediction, 
                               ynames = ynames, xnames = xnames)
  
  # Predictions using ordinary least squares regression
  OLS.Ypred.list <- list()
  for (y in ynames) {
    lmfit <- lm(formula(paste0(y, "~", paste0(xnames, collapse = "+"))), 
                data = as.data.frame(calibration))
    OLS.Ypred.list[[match(y, ynames)]] <- predict.lm(lmfit, newdata = as.data.frame(prediction))
  }
  OLS.Ypred <- do.call("cbind", OLS.Ypred.list)
  colnames(OLS.Ypred) <- ynames
  
  # Predictions using regularised SEM rule (with cross-validation)
  lavcv.Ypred <- lav.predict.y.cv(calidat = calibration, preddat = prediction, 
                                  califit = lav.fit, CV = lav.CV, 
                                  alpha1 = lav.alpha1, alpha2 = lav.alpha2,
                                  K = K, nK = nK, partid = partIDx,
                                  xnames = xnames, ynames = ynames)
  
  # Predictions using elastic net regression (with cross-validation)
  encv.Ypred <- en.predict.y.cv(calidat = calibration, preddat = prediction,
                                alphas = en.alphas, partid = partIDx,
                                xnames = xnames, ynames = ynames)
  
  return() #TODO
  
  # TODO separate times saved for each method? & overall time saved as an attribute?
}

# foo <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = F)
# bar <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = T)
