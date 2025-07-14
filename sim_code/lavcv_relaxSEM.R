## Aditi M. Bhangale
## Last updated: 14 July 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to test the prediction rule in a k-fold cross-validation setting

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

source("part_relaxSEM.R")

# fit model in lavaan----
fitmod <- function(dat, n_x, n_eta_x, n_y,  n_eta_y) {
  obsxnames <- paste0("x", 1:n_x)
  obsynames <- paste0("y", 1:n_y)
  lvxnames <- paste0("eta_x", 1:n_eta_x)
  lvynames <- paste0("eta_y", 1:n_eta_y)
  
  # factor loadings
  if (n_eta_x == 1L) {
    FLx <- paste0(lvxnames, "=~", paste0(obsxnames, collapse = "+"), " \n")
  } else if (n_eta_x == 3L) {
    for (xx in 1:n_eta_x) { # only conditions wherein all factors have equal number of indicators supported for now
      if(xx == 1) {
        FLx1 <- paste0(lvxnames[xx], "=~", paste0(obsxnames[1:(n_x/n_eta_x)], collapse = "+"), " \n")
      } else if (xx == 2) {
        FLx2 <- paste0(lvxnames[xx], "=~", paste0(obsxnames[(n_x/n_eta_x + 1):(n_x-n_x/n_eta_x)], 
                                                  collapse = "+"), " \n")
      } else if (xx == 3) {
        FLx3 <- paste0(lvxnames[xx], "=~", paste0(obsxnames[(n_x-n_x/n_eta_x+1):n_x], 
                                                  collapse = "+"), " \n")
      }
    }
    
    FLx <- paste0(FLx1, FLx2, FLx3, collapse = "\n")
    
  } else stop("only `n_eta_x = 1 or 3` supported for now")
  
  if (n_eta_y == 1L) {
    FLy <- paste0(lvynames, "=~", paste0(obsynames, collapse = "+"), " \n")
  } else stop("only `n_eta_y = 1` supported for now")
  
  FL <- paste0(FLx, FLy, collapse = " \n")
  
  # regression slopes
  if (n_eta_x == 1L & n_eta_y == 1L) {
    B <- paste0(lvynames, "~", lvxnames, " \n")
  } else if (n_eta_x == 3L & n_eta_y == 1L) {
    B <- paste0(lvynames, "~", paste(lvxnames, collapse = "+"), " \n")
  } else stop("only `n_eta_x = 1 or 3` with `n_eta_y = 1` supported for now")
  
  mod <- paste(FL, B, collapse = " \n") # complete model
  
  fit <- sem(mod = mod, dat = dat, meanstructure = T) # factor and indicator variances will be added by default
  
  return(fit)
}

# fitmod(calibration, n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1)

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
lav.predict.y.part <- function(dat, K, partid, 
                               alpha1, alpha2, n_x, n_eta_x, n_y,  n_eta_y,
                               xnames, ynames) {
  partdat <- partition(partid = partid, dat = dat, K = K) # partitioned data
  
  # row and column names
  mat.rows <- do.call("c", 
                      lapply(1:K, function(k) 
                        paste0(k, ".", 1:nrow(partdat[[k]]$test))))
  mat.cols <- apply(expand.grid(ynames, as.character(alpha1), as.character(alpha2)), 
                    1, paste0, collapse = ",") 
  # use as.character() above to paste only 0/1 instead of 0.0 and 1.0
  # because in the for loops, a1/a2 are 0/1 not 0.0/0.1
  sqdevmat <- matrix(NA, nrow(dat), length(alpha1)*length(alpha2)*length(ynames),
                     dimnames = list(mat.rows, mat.cols)) # matrix with squared deviations
  
  for (k in 1:K) {
    for (a1 in alpha1) {
      for (a2 in alpha2) {
        fitpart <- fitmod(dat = partdat[[k]]$train,
                          n_x = n_x, n_eta_x = n_eta_x, 
                          n_y = n_y,  n_eta_y = n_eta_y) # fit to only the training data
        
        predpart <- lav.predict.y(calidat = partdat[[k]]$train,
                                  preddat = partdat[[k]]$test,
                                  califit = fitpart,
                                  alpha1 = a1, alpha2 = a2, xnames = xnames,
                                  ynames = ynames)
        
        sqdevmat[rownames(sqdevmat)[grep(pattern = paste0("^", k, "\\."), 
                                         rownames(sqdevmat))], 
                 paste0(ynames, ",", a1,",", a2)] <- 
          as.matrix((predpart$Ypred - predpart$Ytrue)^2)
      }
    }
  } 
  return(sqdevmat)
}

# t0 <- Sys.time()
# foo <- lav.predict.y.part(dat = calibration, K = 10, partid = part.ids,
#                       alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
#                       xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1, t0, "sec")

#----

# compute RMSEp(r) for each alpha1,alpha2 combination and return alpha1,2 values with min(RMSEp)----
lav.predict.y.alpha <- function(dat, K, partid, 
                                alpha1, alpha2, xnames, ynames) {
  sqdevmat <- lav.predict.y.part(dat = dat, K = K, 
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
        RMSEp.val <- 
          sum(sqdevmat[rownames(sqdevmat)[grep(pattern = paste0(k, "\\."), 
                                               rownames(sqdevmat))],
                       paste0(ynames, ",", a1,",", a2)])/(length(grep(pattern = paste0(k, "\\."), 
                                                                      rownames(sqdevmat)))*length(ynames))
        
        RMSEp[RMSEp$alpha1 == a1 & RMSEp$alpha2 == a2, "RMSEp"] <- RMSEp.val
      }
    }
  }
  
  min.RMSE.val <- RMSEp[which(RMSEp$RMSEp == min(RMSEp$RMSEp)),] # minimum RMSE value and associated alpha1/2
  
  return(list(alpha1 = min.RMSE.val$alpha1, alpha2 = min.RMSE.val$alpha2))  
}

# lav.predict.y.alpha(dat = calibration, K = 10, partid = part.ids,
#                 alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1),
#                 xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----

# prediction rule with cross-validation----
lav.predict.y.cv <- function(calidat, preddat, califit, CV,
                             alpha1, alpha2, K, partid, 
                             xnames, ynames) {
  
  if (CV) {
    if (is.numeric(alpha1) && is.numeric(alpha2) && 
        length(alpha1) > 1L && length(alpha2) > 1L && 
        all(0L <= alpha1) && all(alpha1 <= 1L) && 
        all(0L <= alpha2) && all(alpha2 <= 1L)) {
      alpha.vals <- lav.predict.y.alpha(dat = calidat, K = K, partid = partid, 
                                        alpha1 = alpha1, alpha2 = alpha2,
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
#                     alpha1 = seq(0,1,0.1), alpha2 = seq(0,1,0.1), K = 10,
#                     partid = part.ids, xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1,t0,"sec")

#----
