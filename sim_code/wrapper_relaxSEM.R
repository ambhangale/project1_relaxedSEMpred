## Aditi M. Bhangale
## Last updated: 31 October 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### wrapper function for De Rooij, OLS, lavcv, and encv rules

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

source("lavcv_relaxSEM.R")
source("encv_relaxSEM.R") 

# sampID = 1; nCal = 250; nPred = 250; misspecify = F; lav.CV = T;
# lav.alpha1 = seq(0,1,0.1); lav.alpha2 = seq(0,1,0.1);
# en.alphas = seq(0,1,0.1); K = 10
# xnames = c(paste0("x",1:3), paste0("y",1:4)); ynames = "dem65_sum"; seed = NULL

wrapper.predict.y <- function(sampID, nCal, nPred = 1e4, covmat, lav.CV = TRUE,
                              lav.alpha1 = seq(0,1,0.1), lav.alpha2 = seq(0,1,0.1), 
                              lav.equal.alphas = F,
                              en.alphas = seq(0,1,0.1), K = 10, seed = NULL,
                              save.out = F) {
  
  covmat.attr <- attributes(covmat) # attributes of covmat used to generate data
  
  # xnames and ynames
  allnames <- covmat.attr$dimnames[[1]]
  xnames <- allnames[grep("x", allnames)]
  ynames <- allnames[grep("y", allnames)]
  
  # details about model
  n_x <- covmat.attr$n_x
  n_eta_x <- covmat.attr$n_eta_x
  n_y <- covmat.attr$n_y
  n_eta_y <- covmat.attr$n_eta_y
  
  # details about misspecification
  misspecify <- covmat.attr$misspecify
  miss.part <- covmat.attr$miss.part
  miss.strength <- covmat.attr$miss.strength
  
  t0 <- Sys.time()
  dat <- gendat(sampID = sampID, nCal = nCal, nPred = nPred, 
                covmat = covmat) # always require `sampID` when called in this wrapper function
  calibration <- dat$calibration # calibration set
  prediction <- dat$prediction # prediction set
  
  Ytrue <- prediction[,ynames] # true values of outcome variable(s)
  
  # partition IDs to be used for lav.predict.y.cv and en.predict.y.cv
  partIDx <- partidx(ndat = nCal, sampID = sampID, K = K)
  
  lav.fit <- fitmod(dat = calibration, 
                    n_x = n_x, n_eta_x = n_eta_x,
                    n_y = n_y, n_eta_y = n_eta_y) # lavaan model fitted on complete calibration set
  
  # check if cov.lv and cov.ov are positive definite
  PD.lv <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.lv"))$values > 0), F, T)
  PD.ov <- ifelse(!all(eigen(lavInspect(lav.fit, "cov.ov"))$values > 0), F, T)
  
  # save number of estimated parameters 
  npar <- lav.fit@Fit@npar
  
  # check whether test for exact fit was rejected or not
  if (lav.fit@Fit@converged) {
    fit.measures <- fitMeasures(lav.fit)
    exact.fit    <- ifelse(fit.measures["pvalue"] < .05, "reject", "fail_to_reject")
    df           <- fit.measures["df"]
    CFI          <- fit.measures["cfi"]
    RMSEA        <- fit.measures["rmsea"]
    RMSEA.lowCI  <- fit.measures["rmsea.ci.lower"]
    RMSEA.upCI   <- fit.measures["rmsea.ci.upper"] 
  } else {
    exact.fit <- df <- CFI <- RMSEA <- RMSEA.lowCI <- RMSEA.upCI <- NA
  }
  
  # Predictions using De Rooij et al. (2022) rule
  DeRooij.t0 <- Sys.time()
  DeRooij.Ypred <- lavPredictY(object = lav.fit, newdata = prediction, 
                               ynames = ynames, xnames = xnames)
  DeRooij.t1 <- Sys.time()
  DeRooij.bias <- DeRooij.Ypred - Ytrue
  DeRooij.diff <- difftime(DeRooij.t1, DeRooij.t0, "sec")
  DeRooij.RMSEp <- cbind(method = "DeRooij", PD.lv = PD.lv, PD.ov = PD.ov, npar = npar, 
                         df = df, exact.fit = exact.fit, CFI = CFI, RMSEA = RMSEA,
                         RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                         RMSEp = sqrt(sum((DeRooij.bias)^2)/(length(ynames)*nPred)),
                         runTime = DeRooij.diff)
  DeRooij.RMSEpr <- cbind(method = "DeRooij", PD.lv = PD.lv, PD.ov = PD.ov, npar = npar, 
                          df = df, exact.fit = exact.fit, CFI = CFI, RMSEA = RMSEA,
                          RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                          yname = ynames,
                          RMSEpr = sqrt(colSums(DeRooij.bias^2)/nPred),
                          runTime = DeRooij.diff)
  
  # Predictions using the structural after measurement approach
  SAM.t0 <- Sys.time()
  lav.fit.sam <- fitmod(dat = calibration, 
                        n_x = n_x, n_eta_x = n_eta_x,
                        n_y = n_y, n_eta_y = n_eta_y, SAM = T)
  PD.lv.sam <- ifelse(!all(eigen(lavInspect(lav.fit.sam, "cov.lv"))$values > 0), F, T)
  PD.ov.sam <- ifelse(!all(eigen(lavInspect(lav.fit.sam, "cov.ov"))$values > 0), F, T)
  
  SAM.Ypred <- lavPredictY(object = lav.fit.sam, newdata = prediction, 
                           ynames = ynames, xnames = xnames)
  SAM.t1 <- Sys.time()
  SAM.bias <- SAM.Ypred - Ytrue
  SAM.diff <- difftime(SAM.t1, SAM.t0, "sec")
  SAM.RMSEp <- cbind(method = "SAM", PD.lv = PD.lv.sam, PD.ov = PD.ov.sam,
                     RMSEp = sqrt(sum((SAM.bias)^2)/(length(ynames)*nPred)),
                     runTime = SAM.diff)
  SAM.RMSEpr <- cbind(method = "SAM", PD.lv = PD.lv.sam, PD.ov = PD.ov.sam,
                      yname = ynames,
                      RMSEpr = sqrt(colSums(SAM.bias^2)/nPred),
                      runTime = SAM.diff)
  
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
  OLS.diff <- difftime(OLS.t1, OLS.t0, "sec")
  OLS.RMSEp <- cbind(method = "OLS", 
                     RMSEp = sqrt(sum((OLS.bias)^2)/(length(ynames)*nPred)),
                     runTime = OLS.diff)
  OLS.RMSEpr <- cbind(method = "OLS", yname = ynames, 
                      RMSEpr = sqrt(colSums(OLS.bias^2)/nPred),
                      runTime = OLS.diff)
  
  # Predictions using regularised SEM rule (with cross-validation)
  lavcv.t0 <- Sys.time()
  lavcv.Ypred <- lav.predict.y.cv(calidat = calibration, preddat = prediction, 
                                  califit = lav.fit, CV = lav.CV, 
                                  alpha1 = lav.alpha1, alpha2 = lav.alpha2,
                                  equal.alphas = lav.equal.alphas, 
                                  n_x = n_x, n_eta_x = n_eta_x,
                                  n_y = n_y, n_eta_y = n_eta_y,
                                  K = K, partid = partIDx,
                                  xnames = xnames, ynames = ynames)
  lavcv.t1 <- Sys.time()
  lavcv.bias <- lavcv.Ypred - Ytrue
  lavcv.diff <- difftime(lavcv.t1, lavcv.t0, "sec")
  lavcv.RMSEp <- cbind(method = "lavcv", PD.lv = PD.lv, PD.ov = PD.ov, npar = npar,
                       df = df, exact.fit = exact.fit, CFI = CFI, RMSEA = RMSEA,
                       RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                       lav.fullSRMR = attr(lavcv.Ypred, "fullSRMR"),
                       lav.xxSRMR = attr(lavcv.Ypred, "xxSRMR"), 
                       lav.yySRMR = attr(lavcv.Ypred, "yySRMR"), 
                       lav.yxSRMR = attr(lavcv.Ypred, "yxSRMR"), 
                       lav.alpha1 = attr(lavcv.Ypred, "alpha1"), 
                       lav.alpha2 = attr(lavcv.Ypred, "alpha2"),
                       RMSEp = sqrt(sum((lavcv.bias)^2)/(length(ynames)*nPred)),
                       runTime = lavcv.diff)
  lavcv.RMSEpr <- cbind(method = "lavcv", PD.lv = PD.lv, PD.ov = PD.ov, npar = npar,
                        df = df, exact.fit = exact.fit, CFI = CFI, RMSEA = RMSEA,
                        RMSEA.lowCI = RMSEA.lowCI, RMSEA.upCI = RMSEA.upCI,
                        lav.fullSRMR = attr(lavcv.Ypred, "fullSRMR"),
                        lav.xxSRMR = attr(lavcv.Ypred, "xxSRMR"), 
                        lav.yySRMR = attr(lavcv.Ypred, "yySRMR"), 
                        lav.yxSRMR = attr(lavcv.Ypred, "yxSRMR"), 
                        lav.alpha1 = attr(lavcv.Ypred, "alpha1"), 
                        lav.alpha2 = attr(lavcv.Ypred, "alpha2"),
                        lav.alpha1 = attr(lavcv.Ypred, "alpha1"), 
                        lav.alpha2 = attr(lavcv.Ypred, "alpha2"),
                        yname = ynames,
                        RMSEpr = sqrt(colSums(lavcv.bias^2)/nPred),
                        runTime = lavcv.diff)
  
  # Predictions using elastic net regression (with cross-validation)
  encv.t0 <- Sys.time()
  encv.Ypred <- en.predict.y.cv(calidat = calibration, preddat = prediction,
                                alphas = en.alphas, partid = partIDx, 
                                n_y = n_y,
                                xnames = xnames, ynames = ynames)
  encv.t1 <- Sys.time()
  encv.bias <- encv.Ypred - Ytrue
  encv.diff <- difftime(encv.t1, encv.t0, "sec")
  encv.RMSEp <- cbind(method = "encv", en.alpha = attr(encv.Ypred, "alpha"), 
                      en.lambda = attr(encv.Ypred, "lambda"),
                      RMSEp = sqrt(sum((encv.bias)^2)/(length(ynames)*nPred)),
                      runTime = encv.diff)
  encv.RMSEpr <- cbind(method = "encv", en.alpha = attr(encv.Ypred, "alpha"), 
                       en.lambda = attr(encv.Ypred, "lambda"),
                       yname = ynames,
                       RMSEpr = sqrt(colSums(encv.bias^2)/nPred),
                       runTime = encv.diff)
  
  # save Ytrue and Ypred for all methods
  if (n_y == 1L){
    Ytrue <- as.data.frame(Ytrue)
    colnames(Ytrue) <- ynames
  } 
  colnames(Ytrue) <- paste0("Ytrue.", colnames(Ytrue))
  colnames(DeRooij.Ypred) <- paste0("DeRooij.", colnames(DeRooij.Ypred))
  colnames(SAM.Ypred) <- paste0("SAM.", colnames(SAM.Ypred))
  colnames(OLS.Ypred) <- paste0("OLS.", colnames(OLS.Ypred))
  colnames(lavcv.Ypred) <- paste0("lavcv.", colnames(lavcv.Ypred))
  colnames(encv.Ypred) <- paste0("encv.", colnames(encv.Ypred)) 
  Y <- as.data.frame(cbind(sampID = sampID, nCal = nCal, nPred = nPred, 
             misspecify = misspecify, miss.part = miss.part, miss.strength = miss.strength,
             Ytrue, DeRooij.Ypred, SAM.Ypred, OLS.Ypred, lavcv.Ypred, encv.Ypred))
  
  # save alpha values for lavcv
  lavcv.alphas <- cbind(sampID = sampID, nCal = nCal, nPred = nPred, 
                        misspecify = misspecify, miss.part = miss.part, 
                        miss.strength = miss.strength,
                        alpha1 = attr(lavcv.Ypred, "alpha1"), 
                        alpha2 = attr(lavcv.Ypred, "alpha2"))
  
  # save RMSEp and RMSEpr
  RMSEp <- cbind(sampID = sampID, nCal = nCal, nPred = nPred, 
                 misspecify = misspecify, miss.part = miss.part, miss.strength = miss.strength,
                 Reduce(function(x,y) merge(x, y, all = T), 
                        list (DeRooij.RMSEp, SAM.RMSEp, OLS.RMSEp, lavcv.RMSEp, encv.RMSEp)))
  
  RMSEpr <- cbind(sampID = sampID, nCal = nCal, nPred = nPred, 
                  misspecify = misspecify, miss.part = miss.part, miss.strength = miss.strength,
                  Reduce(function(x,y) merge(x, y, all = T), 
                         list (DeRooij.RMSEpr, SAM.RMSEpr, OLS.RMSEpr, lavcv.RMSEpr, encv.RMSEpr)))
  
  t1 <- Sys.time()
  diff <- difftime(t1, t0, "sec")
  
  final <- list(Y = Y, RMSEp = RMSEp, RMSEpr = RMSEpr, lavcv.alphas = lavcv.alphas)
  
  attr(final, "sampID")        <- sampID
  attr(final, "nCal")          <- nCal
  attr(final, "nPred")         <- nPred
  attr(final, "misspecify")    <- misspecify
  attr(final, "miss.part")     <- miss.part
  attr(final, "miss.strength") <- miss.strength
  attr(final, "lav.CV")        <- lav.CV
  attr(final, "K")             <- K
  attr(final, "PD.lv")         <- PD.lv
  attr(final, "PD.ov")         <- PD.ov
  attr(final, "PD.lv.sam")     <- PD.lv.sam
  attr(final, "PD.ov.sam")     <- PD.ov.sam
  attr(final, "exact.fit")     <- exact.fit
  attr(final, "df")            <- df
  attr(final, "CFI")           <- CFI
  attr(final, "RMSEA")         <- RMSEA
  attr(final, "RMSEA.lowCI")   <- RMSEA.lowCI
  attr(final, "RMSEA.upCI")    <- RMSEA.upCI
  attr(final, "xnames")        <- xnames
  attr(final, "ynames")        <- ynames
  attr(final, "lav.alpha1")    <- attr(lavcv.Ypred, "alpha1")
  attr(final, "lav.alpha2")    <- attr(lavcv.Ypred, "alpha2")
  attr(final, "en.alpha")      <- attr(encv.Ypred, "alpha")
  attr(final, "en.lambda")     <- attr(encv.Ypred, "lambda")
  attr(final, "seed")          <- ifelse(!is.null(seed), seed, NA)
  attr(final, "runtime")       <- diff
  
  if(save.out) saveRDS(final, paste0("ID", sampID, ".rds")) # FIXME just placeholder for now
  
  return(final) 
}

# foo <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = F)
# bar <- wrapper.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = T)
