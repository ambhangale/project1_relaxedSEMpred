## Aditi M. Bhangale
## Last updated: 27 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### function(s) file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "datagen_RDA_CFA.R"))

fitmod <- function(dat) {
  mod <- '
# factor loadings
F1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7

# factor variance
F1 ~~ 1*F1

# item (co)variances
x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
x6 ~~ x6
x7 ~~ x7

# factor mean
F1 ~ 0*1

# item means
x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
x6 ~ 1
x7 ~ 1
'
  
  fit <- lavaan(mod, data = dat, meanstructure = T)
  
  return(fit)
}

testrule <- function(nCal, nPred, misspecify,
                     alpha1, alpha2, 
                     xnames = paste0("x", 4:7), ynames = paste0("x", 1:3), 
                     seed = 10824) {
  dat <- gendat(nCal = nCal, nPred = nPred, 
                misspecify = misspecify, seed = seed)
  
  calibration <- dat$calibration
  prediction  <- dat$prediction
  
  # rescale covariance matrix of training set to make consistent with `lavaan`
  S    <- (cov(calibration)*(nrow(calibration)-1)) / nrow(calibration) 
  S_xx <- S[xnames, xnames]
  S_yx <- S[ynames, xnames] ## using _yx to avoid using t()
  
  # values from prediction dataset
  X0    <- prediction[,xnames] # to be inputted into formulae
  Ytest <- prediction[,ynames] # original Y values from prediction dataset
  
  fit <- fitmod(dat = calibration)
  ImpliedStats <- lavInspect(fit, "implied")
  Sigma_xx     <- ImpliedStats$cov[xnames, xnames]
  Sigma_yx     <- ImpliedStats$cov[ynames, xnames] ## using _yx to avoid using t()
  Mu_x         <- ImpliedStats$mean[xnames]
  Mu_y         <- ImpliedStats$mean[ynames]
  
  PD.lv <- ifelse(!all(eigen(lavInspect(fit, "cov.lv"))$values > 0), F, T)
  PD.ov <- ifelse(!all(eigen(lavInspect(fit, "cov.ov"))$values > 0), F, T)
  
  if (0L <= alpha1 && alpha1 <= 1L && 0L <= alpha2 && alpha2 <= 1L) {
    Ypred <- t(Mu_y + ((1-alpha2)*Sigma_yx + alpha2*S_yx) %*% 
                 solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                 (t(X0) - Mu_x))
  } else {
    stop("specify values between 0 and 1 for `alpha1` and `alpha2`")
  }
  
  bias <- Ypred - Ytest
  
  RMSEpr.result <- as.data.frame(cbind(nCal     = nCal, 
                                       nPred    = nPred, 
                                       alpha1   = alpha1, 
                                       alpha2   = ifelse(!is.null(alpha2), alpha2, NA),
                                       PD.lv    = PD.lv,
                                       PD.ov    = PD.ov,
                                       meanBias = colMeans(bias),
                                       RMSEpr   = sqrt(colMeans((bias)^2)),
                                       yname    = ynames))
  
  RMSEp.result <- as.data.frame(cbind(nCal       = nCal, 
                                      nPred      = nPred, 
                                      alpha1     = alpha1, 
                                      alpha2     = ifelse(!is.null(alpha2), alpha2, NA),
                                      PD.lv      = PD.lv,
                                      PD.ov      = PD.ov,
                                      misspecify = misspecify,
                                      RMSEp      = sqrt(sum((bias)^2)/(length(ynames)*nPred))))
  
  final <- list(Ypred = Ypred, Ytest = Ytest, bias = bias,
                RMSEpr.result = RMSEpr.result, RMSEp.result = RMSEp.result)
  # save all arguments as attributes, just in case we need them later
  attr(final, "nCal")         <- nCal
  attr(final, "nPred")        <- nPred
  attr(final, "misspecify")   <- misspecify
  attr(final, "alpha1")       <- alpha1
  attr(final, "alpha2")       <- alpha2
  attr(final, "PD.lv")        <- PD.lv
  attr(final, "PD.ov")        <- PD.ov
  attr(final, "xnames")       <- xnames
  attr(final, "ynames")       <- ynames
  
  return(final)
}

# testrule(nCal = 100, nPred = 100, misspecify = F, alpha1 = 0, alpha2 = 0)
# testrule(nCal = 100, nPred = 100, misspecify = F, alpha1 = 1, alpha2 = 0)
# testrule(nCal = 100, nPred = 100, misspecify = T, alpha1 = 0, alpha2 = 1)
# testrule(nCal = 100, nPred = 100, misspecify = T, alpha1 = 1, alpha2 = 1)
