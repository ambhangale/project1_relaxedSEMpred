## Aditi M. Bhangale
## Last updated: 1 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### functions to test the prediction rule in a k-fold cross-validation setting

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_CFA")

library(here)
source(here("RDA_CFA", "datagen_RDA_CFA.R"))

# dat <- gendat(sampID = 1, nCal = 250, nPred = 250, misspecify = F, seed = 10824) # dummy data for now
# calibration <- dat$calibration; prediction <- dat$prediction

# function to partition data into K parts----
partition <- function(dat, K, nK, seed) { # input here is the calibration data
  
  set.seed(seed)
  
  partvec <- rep(1:K, nK)
  
  # randomly assign observations to one of the K parts such that there are nK 
  # observations per part
  partidx <- sample(partvec, nrow(dat), replace = F)
  
  final <- list()
  
  for (k in 1:K) {
    assign(paste0("test", k), dat[partidx == k,])
    assign(paste0("train", k), dat[partidx != k,])
    
    final[[paste0("part",k)]] <- list(train = get(paste0("train", k)),
                                      test = get(paste0("test", k)))
  }
  
  return(final)
}

# partition(calibration, K = 10, nK = 25) # test function

#----

# fit model in lavaan (from 'func_RDA_CFA.R' file)----
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

#----

# prediction rule----
predict.y <- function(calidat, preddat, califit, 
                      alpha1, alpha2, xnames, ynames) {
  S <- (cov(calidat)*(nrow(calidat)-1)) / nrow(calidat) 
  S_xx <- S[xnames, xnames]
  S_yx <- S[ynames, xnames] ## using _yx to avoid using t()
  
  # values from prediction dataset
  X0    <- preddat[,xnames] # to be inputted into formulae
  Ytrue <- preddat[,ynames] # true Y values from prediction dataset
  
  ImpliedStats <- lavInspect(califit, "implied")
  Sigma_xx     <- ImpliedStats$cov[xnames, xnames]
  Sigma_yx     <- ImpliedStats$cov[ynames, xnames] ## using _yx to avoid using t()
  Mu_x         <- ImpliedStats$mean[xnames]
  Mu_y         <- ImpliedStats$mean[ynames]
  
  PD.lv <- ifelse(!all(eigen(lavInspect(califit, "cov.lv"))$values > 0), F, T)
  PD.ov <- ifelse(!all(eigen(lavInspect(califit, "cov.ov"))$values > 0), F, T)
  
  if (0L <= alpha1 && alpha1 <= 1L && 0L <= alpha2 && alpha2 <= 1L) {
    Ypred <- t(Mu_y + ((1-alpha2)*Sigma_yx + alpha2*S_yx) %*% 
                 solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                 (t(X0) - Mu_x))
  } else {
    stop("specify values between 0 and 1 for `alpha1` and `alpha2`")
  }
  
  final <- list(Ypred = Ypred, Ytrue = Ytrue)
  
  attr(final, "PD.lv")        <- PD.lv
  attr(final, "PD.ov")        <- PD.ov
  
  return(final) 
  
}

# fit <- fitmod(dat = calibration)
# predict.y(calibration, prediction, fit, alpha1 = 0.5, alpha2 = 0.3)

#----

# prediction for the K partitions----
predict.y.part <- function(dat, K, nK, 
                           alpha1, alpha2, xnames, ynames, seed) { #TODO is seed argument necessary?
  partdat <- partition(dat = dat, K = K, nK = nK, seed = seed) # partitioned data
  
  mat.rows <- do.call("c", lapply(1:K, 
                                  function(k) apply(expand.grid(k, 1:nK), 1, 
                                                    paste0, collapse = ".")))
  mat.cols <- apply(expand.grid(ynames, as.character(alpha1), as.character(alpha2)), 
                    1, paste0, collapse = ",") 
  # use as.charcater() above to paste only 0/1 instead of 0.0 and 1.0
  # because in the for loops, a1/a2 are 0/1 not 0.0/0.1
  biasmat <- matrix(NA, K*nK, length(alpha1)*length(alpha2)*length(ynames),
                    dimnames = list(mat.rows, mat.cols)) # matrix with predicted values of partitioned data
 
 for (k in 1:K) {
   for (a1 in alpha1) {
     for (a2 in alpha2) {
       fitpart <- fitmod(dat = partdat[[k]]$train) # fit to only the training data
       
       predpart <- predict.y(calidat = partdat[[k]]$train,
                             preddat = partdat[[k]]$test,
                             califit = fitpart,
                             alpha1 = a1, alpha2 = a2, xnames = xnames,
                             ynames = ynames)
       
       biasmat[paste0(k, ".", 1:nK), paste0(ynames, ",", a1,",", a2)] <- 
         as.matrix(predpart$Ypred - predpart$Ytrue)
     }
   }
 } 
  return(biasmat)
}

# t0 <- Sys.time()
# foo <- predict.y.part(dat = calibration)
# t1 <- Sys.time()
# diff <- difftime(t1, t0, "sec")

#----

# compute RMSEp(r) for each alpha1,alpha2 combination and return alpha1,2 values with min(RMSEp)----
predict.y.alpha <- function(dat, K, nK, 
                            alpha1, alpha2, xnames, ynames, seed) {
  biasmat <- predict.y.part(dat = dat, K = K, nK = nK, alpha1 = alpha1,
                            alpha2 = alpha2, xnames = xnames, ynames = ynames, 
                            seed = seed)
  
  RMSEp  <- expand.grid(alpha1 = alpha1, alpha2 = alpha2, RMSEp = NA)
  # RMSEp <- matrix(NA, length(alpha1)*length(alpha2), 3,
  #                 dimnames = list(NULL, c("alpha1", "alpha2", "RMSEp")))
  
  # TODO only doing RMSEp for now, but if necessary, can add RMSEpr later
  
  for (k in 1:K) {
    for (a1 in alpha1) {
      for (a2 in alpha2) {
        
        RMSEp.val <- sum((biasmat[paste0(k, ".", 1:nK), 
                              paste0(ynames, ",", a1,",", a2)])^2)/(nK*length(ynames))
        
        RMSEp[RMSEp$alpha1 == a1 & RMSEp$alpha2 == a2, "RMSEp"] <- RMSEp.val
      }
    }
  }
  
  min.RMSE.val <- RMSEp[which(RMSEp$RMSEp == min(RMSEp$RMSEp)),] # minimum RMSE value and associated alpha1/2
 
 return(list(alpha1 = min.RMSE.val$alpha1, alpha2 = min.RMSE.val$alpha2))  
}

# predict.y.alpha(dat = calibration)

#----

# prediction rule with cross-validation----
predict.y.cv <- function(sampID, nCal, nPred, misspecify,
                         alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
                         K = 10, nK = NULL, 
                         xnames = paste0("x", 4:7), ynames = paste0("x", 1:3),
                         seed = 10824) {
  dat <- gendat(sampID = sampID, nCal = nCal, nPred = nPred, misspecify = misspecify, seed = seed)
  calibration <- dat$calibration # calibration set
  prediction  <- dat$prediction # prediction set
  
  nK <- ifelse(is.null(nK), nCal/K, nK) # compute nK manually if left blank
  
  alpha.vals <- predict.y.alpha(dat = calibration, K = K, nK = nK,
                                alpha1 = alpha1, alpha2 = alpha2,
                                xnames = xnames, ynames = ynames, 
                                seed = seed)
  
  fit <- fitmod(dat = calibration) # model fitted on complete calibration set
  
  predVals <- predict.y(calidat = calibration, preddat = prediction,
                        califit = fit, alpha1 = alpha.vals$alpha1,
                        alpha2 = alpha.vals$alpha2,
                        xnames = xnames, ynames = ynames)
  
  bias <- predVals$Ypred - predVals$Ytrue
  
  PD.lv <- attributes(predVals)$PD.lv
  PD.ov <- attributes(predVals)$PD.ov
  
  RMSEpr.result <- as.data.frame(cbind(sampID   = sampID, 
                                       nCal     = nCal, 
                                       nPred    = nPred, 
                                       alpha1   = alpha.vals$alpha1, 
                                       alpha2   = alpha.vals$alpha2,
                                       PD.lv    = PD.lv,
                                       PD.ov    = PD.ov,
                                       meanBias = colMeans(bias),
                                       RMSEpr   = sqrt(colMeans((bias)^2)),
                                       yname    = ynames))
  
  RMSEp.result <- as.data.frame(cbind(sampID     = sampID, 
                                      nCal       = nCal, 
                                      nPred      = nPred, 
                                      alpha1     = alpha.vals$alpha1, 
                                      alpha2     = alpha.vals$alpha2,
                                      PD.lv      = PD.lv,
                                      PD.ov      = PD.ov,
                                      misspecify = misspecify,
                                      RMSEp      = sqrt(sum((bias)^2)/(length(ynames)*nPred))))
  
  final <- list(Ypred = predVals$Ypred, Ytrue = predVals$Ytrue, bias = bias,
                RMSEpr.result = RMSEpr.result, RMSEp.result = RMSEp.result)
  
  # save all arguments as attributes, just in case we need them later
  attr(final, "sampID")       <- sampID
  attr(final, "nCal")         <- nCal
  attr(final, "nPred")        <- nPred
  attr(final, "misspecify")   <- misspecify
  attr(final, "alpha1")       <- alpha.vals$alpha1
  attr(final, "alpha2")       <- alpha.vals$alpha2
  attr(final, "K")            <- K
  attr(final, "nK")           <- nK
  attr(final, "PD.lv")        <- PD.lv
  attr(final, "PD.ov")        <- PD.ov
  attr(final, "xnames")       <- xnames
  attr(final, "ynames")       <- ynames
  attr(final, "seed")         <- seed
  
  return(final)
  
}

# t0 <- Sys.time()
# bar <- predict.y.cv(sampID = 1, nCal = 250, nPred = 250, misspecify = F)
# t1 <- Sys.time()
# diff <- difftime(t1,t0,"sec")

#----


