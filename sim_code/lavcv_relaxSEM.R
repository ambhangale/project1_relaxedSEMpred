## Aditi M. Bhangale
## Last updated: 22 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to test the prediction rule in a k-fold cross-validation setting

library(here)
source(here("sim_code", "part_relaxSEM.R"))

# fit model in lavaan----
fitmod <- function(dat) {
  mod <- cs.mod
  
  fit <- sem(mod, data = dat, meanstructure = T)
  
  return(fit)
}

# fitmod(calibration)

#----

# prediction rule----
lav.predict.y <- function(calidat, preddat, califit, 
                      alpha1, alpha2, xnames, ynames) {
  # rescale covariance matrix to scale as `lavaan` does
  S <- (cov(calidat)*(nrow(calidat)-1)) / nrow(calidat)
  S_xx <- S[xnames, xnames]
  S_yx <- S[ynames, xnames] ## using _yx to avoid using t()
  
  # values from prediction dataset
  X0    <- preddat[,xnames] # to be inputted into formulae
  Ytrue <- preddat[,ynames] # true Y values from prediction dataset
  # Ytrue will be used to compute squared deviations in `lav.predict.y.part()`
  
  ImpliedStats <- lavInspect(califit, "implied")
  Sigma_xx     <- ImpliedStats$cov[xnames, xnames]
  Sigma_yx     <- ImpliedStats$cov[ynames, xnames] ## using _yx to avoid using t()
  Mu_x         <- ImpliedStats$mean[xnames]
  Mu_y         <- ImpliedStats$mean[ynames]
  
  if (0L <= alpha1 && alpha1 <= 1L && 0L <= alpha2 && alpha2 <= 1L) {
    Ypred <- t(Mu_y + ((1-alpha2)*Sigma_yx + alpha2*S_yx) %*% 
                 solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                 (t(X0) - Mu_x))
  } else {
    stop("specify values between 0 and 1 for `alpha1` and `alpha2`")
  }
  
  final <- list(Ypred = Ypred, Ytrue = Ytrue)
  
  return(final) 
  
}

# fit <- fitmod(dat = calibration)
# lav.predict.y(calibration, prediction, fit, alpha1 = 0.5, alpha2 = 0.3,
#           xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----

# prediction for the K partitions----
lav.predict.y.part <- function(dat, K, nK, partid, 
                           alpha1, alpha2, xnames, ynames) {
  partdat <- partition(partid = partid, dat = dat, K = K) # partitioned data
  
  # row and column names
  # FIXME JDK: would abandon this in favour of array, see previous comment 
  mat.rows <- do.call("c", lapply(1:K, 
                                  function(k) apply(expand.grid(k, 1:nK), 1, 
                                                    paste0, collapse = ".")))
  mat.cols <- apply(expand.grid(ynames, as.character(alpha1), as.character(alpha2)), 
                    1, paste0, collapse = ",") 
  # use as.character() above to paste only 0/1 instead of 0.0 and 1.0
  # because in the for loops, a1/a2 are 0/1 not 0.0/0.1
  # FIXME JDK: below, K*nK will not always work (when n is not divisible by K) 
  # FIXME JDK: I would recommend using an array instead of dim (n, length(alpha1), length(alpha2), length(ynames)) 
  sqdevmat <- matrix(NA, K*nK, length(alpha1)*length(alpha2)*length(ynames),
                     dimnames = list(mat.rows, mat.cols)) # matrix with squared deviations
  
  for (k in 1:K) {
    for (a1 in alpha1) {
      for (a2 in alpha2) {
        fitpart <- fitmod(dat = partdat[[k]]$train) # fit to only the training data
        
        predpart <- lav.predict.y(calidat = partdat[[k]]$train,
                              preddat = partdat[[k]]$test,
                              califit = fitpart,
                              alpha1 = a1, alpha2 = a2, xnames = xnames,
                              ynames = ynames)
        
        sqdevmat[paste0(k, ".", 1:nK), paste0(ynames, ",", a1,",", a2)] <- 
          as.matrix((predpart$Ypred - predpart$Ytrue)^2)
      }
    }
  } 
  return(sqdevmat)
}

# t0 <- Sys.time()
# foo <- lav.predict.y.part(dat = calibration, K = 10, nK = 25, partid = part.ids,
#                       alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
#                       xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1, t0, "sec")

#----

# compute RMSEp(r) for each alpha1,alpha2 combination and return alpha1,2 values with min(RMSEp)----
lav.predict.y.alpha <- function(dat, K, nK, partid, 
                            alpha1, alpha2, xnames, ynames) {
  sqdevmat <- lav.predict.y.part(dat = dat, K = K, nK = nK, 
                                 partid = partid, alpha1 = alpha1, alpha2 = alpha2, 
                                 xnames = xnames, ynames = ynames)
  
  RMSEp  <- expand.grid(alpha1 = alpha1, alpha2 = alpha2, RMSEp = NA)
  # RMSEp <- matrix(NA, length(alpha1)*length(alpha2), 3,
  #                 dimnames = list(NULL, c("alpha1", "alpha2", "RMSEp")))
  
  # TODO only doing RMSEp for now, but if necessary, can add RMSEpr later
  
  for (k in 1:K) {
    for (a1 in alpha1) {
      for (a2 in alpha2) {
        
        # compute single RMSE value per alpha1-alpha2 combination
        RMSEp.val <- sum(sqdevmat[paste0(k, ".", 1:nK), 
                                  paste0(ynames, ",", a1,",", a2)])/(nK*length(ynames))
        
        RMSEp[RMSEp$alpha1 == a1 & RMSEp$alpha2 == a2, "RMSEp"] <- RMSEp.val
      }
    }
  }
  
  min.RMSE.val <- RMSEp[which(RMSEp$RMSEp == min(RMSEp$RMSEp)),] # minimum RMSE value and associated alpha1/2
  
  return(list(alpha1 = min.RMSE.val$alpha1, alpha2 = min.RMSE.val$alpha2))  
}

# lav.predict.y.alpha(dat = calibration, K = 10, nK = 25, partid = part.ids,
#                 alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
#                 xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----

# prediction rule with cross-validation----
lav.predict.y.cv <- function(calidat, preddat, califit, CV,
                             alpha1, alpha2,
                             K, nK, partid, 
                             xnames, ynames) {
  
  if (CV) {
    if (is.numeric(alpha1) && is.numeric(alpha2) && 
        length(alpha1) > 1L && length(alpha2) > 1L && 
        all(0L <= alpha1) && all(alpha1 <= 1L) && 
        all(0L <= alpha2) && all(alpha2 <= 1L)) {
      alpha.vals <- lav.predict.y.alpha(dat = calidat, K = K, nK = nK,
                                        partid = partid, alpha1 = alpha1, alpha2 = alpha2,
                                        xnames = xnames, ynames = ynames)
    } else {
      stop("specify numeric vectors with length > 1L and only containing values between 
         0 and 1 for `alpha1` and `alpha2` when `CV = TRUE`")
    }
  } else {
    if (is.numeric(alpha1) && is.numeric(alpha2) && 
        length(alpha1) == 1L && length(alpha2) == 1L &&
        0L <= alpha1 && alpha1 <= 1L && 0L <= alpha2 && alpha2 <= 1L) {
      alpha.vals <- list(alpha1 = alpha1, alpha2 = alpha2)
    } else {
      stop("specify numeric values between 0 and 1 for `alpha1` and `alpha2` when `CV = FALSE`")
    }
  }
  
  Ypred <- lav.predict.y(calidat = calidat, preddat = preddat,
                         califit = califit, alpha1 = alpha.vals$alpha1,
                         alpha2 = alpha.vals$alpha2,
                         xnames = xnames, ynames = ynames)$Ypred
  
  attr(Ypred, "CV")     <- CV # was cross-validation performed?
  attr(Ypred, "alpha1") <- alpha.vals$alpha1 # for xx part
  attr(Ypred, "alpha2") <- alpha.vals$alpha2 # for xy part
  
  return(Ypred)
  
}

# t0 <- Sys.time()
# bar <- lav.predict.y.cv(calidat = calibration, preddat = prediction,
#                         califit = fit, CV = T,
#                     alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1), K = 10, nK = 25,
#                     partid = part.ids, xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1,t0,"sec")

#----
