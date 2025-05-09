## Aditi M. Bhangale
## Last updated: 9 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### wrapper function for De Rooij, OLS, lavcv, and encv rules

library(here)
source(here("sim_code", "lavcv_relaxSEM.R"))
source(here("sim_code", "encv_relaxSEM.R")) 

# sampID = 1; nCal = 250; nPred = 250; misspecify = F; lav.CV = T;
# lav.alpha1 = seq(0,1,0.1); lav.alpha2 = seq(0,1,0.1);
# en.alphas = seq(0,1,0.1); K = 10; nK = NULL;
# xnames = c(paste0("x",1:3), paste0("y",1:4)); ynames = "dem65_sum"; seed = NULL

wrapper.predict.y <- function(sampID, nCal, nPred, misspecify, lav.CV = TRUE,
                              lav.alpha1 = seq(0,1,0.1), lav.alpha2 = seq(0,1,0.1), 
                              en.alphas = seq(0,1,0.1), K = 10, nK = NULL, 
                              xnames = c(paste0("x",1:3), paste0("y",1:4)), 
                              ynames = "dem65_sum",
                              seed = NULL) {
  t0 <- Sys.time()
  dat <- gendat(sampID = sampID, nCal = nCal, nPred = nPred, 
                misspecify = misspecify) # always require `sampID` when called in this wrapper function
  calibration <- dat$calibration # calibration set
  prediction <- dat$prediction # prediction set
  
  Ytrue <- prediction[,ynames] # true values of outcome variable(s)
  
  nK <- ifelse(is.null(nK), nCal/K, nK) # compute nK manually if left blank
  
  # partition IDs to be used for lav.predict.y.cv and en.predict.y.cv
  partIDx <- partidx(ndat = nCal, sampID = sampID, K = K, nK = nK)
  
  lav.fit <- fitmod(dat = calibration) # lavaan model fitted on complete calibration set
  
  # check if cov.lv and cov.ov are positive definite
  PD.lv <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.lv"))$values > 0), F, T)
  PD.ov <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.ov"))$values > 0), F, T)
  
  # check whether test for exact fit was rejected or not
  fit.measures <- fitMeasures(lav.fit)
  exact.fit    <- ifelse(fit.measures["pvalue"] < .05, "reject", "fail_to_reject")
  RMSEA        <- fit.measures["rmsea"]
  RMSEA.lowCI  <- fit.measures["rmsea.ci.lower"]
  RMSEA.upCI   <- fit.measures["rmsea.ci.upper"]
  
  # Predictions using De Rooij et al. (2022) rule
  DeRooij.t0 <- Sys.time()
  DeRooij.Ypred <- lavPredictY(object = lav.fit, newdata = prediction, 
                               ynames = ynames, xnames = xnames)
  DeRooij.t1 <- Sys.time()
  DeRooij.bias <- DeRooij.Ypred - Ytrue
  DeRooij.meanBias <- colMeans(DeRooij.bias)
  DeRooij.diff <- difftime(DeRooij.t1, DeRooij.t0, "sec")
  DeRooij.RMSEp <- cbind(method = "DeRooij", PD.lv = PD.lv, PD.ov = PD.ov, 
                         exact.fit = exact.fit, RMSEA = RMSEA,
                         RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                         RMSEp = sqrt(sum((DeRooij.bias)^2)/(length(ynames)*nPred)),
                         runTime = DeRooij.diff)
  DeRooij.RMSEpr <- cbind(method = "DeRooij", PD.lv = PD.lv, PD.ov = PD.ov,
                          exact.fit = exact.fit, RMSEA = RMSEA,
                          RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                          yname = ynames, 
                          meanBias = DeRooij.meanBias,
                          RMSEpr = sqrt(DeRooij.meanBias^2),
                          runTime = DeRooij.diff)
  
  
  # Predictions using ordinary least squares regression
  OLS.t0 <- Sys.time()
  OLS.Ypred.list <- list()
  for (y in ynames) {
    lmfit <- lm(formula(paste0(y, "~", paste0(xnames, collapse = "+"))), 
                data = as.data.frame(calibration))
    OLS.Ypred.list[[match(y, ynames)]] <- predict.lm(lmfit, newdata = as.data.frame(prediction))
  }
  OLS.Ypred <- do.call("cbind", OLS.Ypred.list)
  colnames(OLS.Ypred) <- ynames
  OLS.t1 <- Sys.time()
  OLS.bias <- OLS.Ypred - Ytrue
  OLS.meanBias <- colMeans(OLS.bias)
  OLS.diff <- difftime(OLS.t1, OLS.t0, "sec")
  OLS.RMSEp <- cbind(method = "OLS", 
                     RMSEp = sqrt(sum((OLS.bias)^2)/(length(ynames)*nPred)),
                     runTime = OLS.diff)
  OLS.RMSEpr <- cbind(method = "OLS", yname = ynames, meanBias = OLS.meanBias, 
                      RMSEpr = sqrt(OLS.meanBias^2),
                      runTime = OLS.diff)
  
  # Predictions using regularised SEM rule (with cross-validation)
  lavcv.t0 <- Sys.time()
  lavcv.Ypred <- lav.predict.y.cv(calidat = calibration, preddat = prediction, 
                                  califit = lav.fit, CV = lav.CV, 
                                  alpha1 = lav.alpha1, alpha2 = lav.alpha2,
                                  K = K, nK = nK, partid = partIDx,
                                  xnames = xnames, ynames = ynames)
  lavcv.t1 <- Sys.time()
  lavcv.bias <- lavcv.Ypred - Ytrue
  lavcv.meanBias <- colMeans(lavcv.bias)
  lavcv.diff <- difftime(lavcv.t1, lavcv.t0, "sec")
  lavcv.RMSEp <- cbind(method = "lavcv", PD.lv = PD.lv, PD.ov = PD.ov, 
                       exact.fit = exact.fit, RMSEA = RMSEA,
                       RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                       lav.alpha1 = attr(lavcv.Ypred, "alpha1"), 
                       lav.alpha2 = attr(lavcv.Ypred, "alpha2"),
                       RMSEp = sqrt(sum((lavcv.bias)^2)/(length(ynames)*nPred)),
                       runTime = lavcv.diff)
  lavcv.RMSEpr <- cbind(method = "lavcv", PD.lv = PD.lv, PD.ov = PD.ov, 
                        exact.fit = exact.fit, RMSEA = RMSEA,
                        RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                        lav.alpha1 = attr(lavcv.Ypred, "alpha1"), 
                        lav.alpha2 = attr(lavcv.Ypred, "alpha2"),
                        yname = ynames,
                        meanBias = lavcv.meanBias,
                        RMSEpr = sqrt(lavcv.meanBias^2),
                        runTime = lavcv.diff)
  
  # Predictions using elastic net regression (with cross-validation)
  encv.t0 <- Sys.time()
  encv.Ypred <- en.predict.y.cv(calidat = calibration, preddat = prediction,
                                alphas = en.alphas, partid = partIDx,
                                xnames = xnames, ynames = ynames)
  encv.t1 <- Sys.time()
  encv.bias <- encv.Ypred - Ytrue
  encv.meanBias <- colMeans(encv.bias)
  encv.diff <- difftime(encv.t1, encv.t0, "sec")
  encv.RMSEp <- cbind(method = "encv", en.alpha = attr(encv.Ypred, "alpha"), 
                      en.lambda = attr(encv.Ypred, "lambda"),
                      RMSEp = sqrt(sum((encv.bias)^2)/(length(ynames)*nPred)),
                      runTime = encv.diff)
  encv.RMSEpr <- cbind(method = "encv", en.alpha = attr(encv.Ypred, "alpha"), 
                       en.lambda = attr(encv.Ypred, "lambda"),
                       yname = ynames,
                       meanBias = encv.meanBias,
                       RMSEpr = sqrt(encv.meanBias^2),
                       runTime = encv.diff)
  
  RMSEp <- cbind(sampID = sampID, nCal = nCal, nPred = nPred, misspecify = misspecify,
                 Reduce(function(x,y) merge(x, y, all = T), 
                        list (DeRooij.RMSEp, OLS.RMSEp, lavcv.RMSEp, encv.RMSEp)))
  
  RMSEpr <- cbind(sampID = sampID, nCal = nCal, nPred = nPred, misspecify = misspecify,
                 Reduce(function(x,y) merge(x, y, all = T), 
                        list (DeRooij.RMSEpr, OLS.RMSEpr, lavcv.RMSEpr, encv.RMSEpr)))
  
  t1 <- Sys.time()
  diff <- difftime(t1, t0, "sec")
  
  final <- list(RMSEp = RMSEp, RMSEpr = RMSEpr)
  
  attr(final, "sampID")       <- sampID
  attr(final, "nCal")         <- nCal
  attr(final, "nPred")        <- nPred
  attr(final, "misspecify")   <- misspecify
  attr(final, "lav.CV")       <- lav.CV
  attr(final, "K")            <- K
  attr(final, "nK")           <- nK
  attr(final, "PD.lv")        <- PD.lv
  attr(final, "PD.ov")        <- PD.ov
  attr(final, "exact.fit")    <- exact.fit
  attr(final, "RMSEA")        <- RMSEA
  attr(final, "RMSEA.lowCI")  <- RMSEA.lowCI
  attr(final, "RMSEA.upCI")   <- RMSEA.upCI
  attr(final, "xnames")       <- xnames
  attr(final, "ynames")       <- ynames
  attr(final, "lav.alpha1")   <- attr(lavcv.Ypred, "alpha1")
  attr(final, "lav.alpha2")   <- attr(lavcv.Ypred, "alpha2")
  attr(final, "en.alpha")     <- attr(encv.Ypred, "alpha")
  attr(final, "en.lambda")    <- attr(encv.Ypred, "lambda")
  attr(final, "seed")         <- ifelse(!is.null(seed), seed, NA)
  attr(final, "runtime")      <- diff
  
  return(final) 
}

# foo <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = F)
# bar <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = T)
