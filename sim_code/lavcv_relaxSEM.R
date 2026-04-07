## Aditi M. Bhangale
## Last updated: 7 April 2026

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### functions to test the prediction rule in a k-fold cross-validation setting

# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/sim_code")

source("part_relaxSEM.R")
source("SRMR_relaxSEM.R")

library(lavaan)

# fit model in lavaan----
fitmod <- function(dat, n_x, n_eta_x, n_y, n_eta_y, SAM = FALSE) {
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
  
  fit <- if (!SAM) {
    sem(mod = mod, dat = dat, meanstructure = T, 
        se = "none", baseline = F) # use ULI constraint
  } else {
    sam(mod = mod, dat = dat, meanstructure = T, 
        se = "none", baseline = F) # SAM, use ULI constraint
  }
  
  return(fit)
}

# fitmod(calibration, n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1)

#----

# prediction rule----
lav.predict.y <- function(preddat, califit, 
                      gamma1, gamma2, xnames, ynames, srmr = F) {
  # extract sample statistics from fitted model object
  SampStats <- lavInspect(califit, "sampstat")
  S         <- SampStats$cov
  S_xx      <- S[xnames, xnames]
  S_xy      <- S[xnames, ynames]
  S_yy      <- S[ynames, ynames]
  # M         <- SampStats$mean
  # M_x       <- SampStats$mean[xnames]
  # M_y       <- SampStats$mean[ynames]
  
  # values from prediction dataset
  X0    <- preddat[,xnames] # to be inputted into formulae
  Ytrue <- preddat[,ynames] # true Y values from prediction dataset
  # Ytrue will be used to compute squared deviations in `lav.predict.y.part()`
  
  # model implied mean vector and covariance matrix
  ImpliedStats <- lavInspect(califit, "implied")
  Sigma        <- ImpliedStats$cov
  Sigma_xx     <- Sigma[xnames, xnames]
  Sigma_xy     <- Sigma[xnames, ynames]
  Sigma_yy     <- Sigma[ynames, ynames]
  Mu           <- ImpliedStats$mean
  Mu_x         <- Mu[xnames]
  Mu_y         <- Mu[ynames]
  
  if (0L <= gamma1 && gamma1 <= 1L && 0L <= gamma2 && gamma2 <= 1L) {
    vec_one <- rep(1, nrow(preddat))
    Mu_x_mat <- vec_one %*% t(Mu_x)
    Mu_y_mat <- vec_one %*% t(Mu_y)
    
    Ypred <- Mu_y_mat + ((X0 - Mu_x_mat) %*%
                           chol2inv(chol((1-gamma1)*Sigma_xx + gamma1*S_xx)) %*%
                           ((1-gamma2)*Sigma_xy + gamma2*S_xy))
  } else {
    stop("specify values between 0 and 1 for `gamma1` and `gamma2`")
  }
  
  if (srmr) { # currently only implemented when `include.mean = F` as mean structure is saturated 
    fullSRMR  <- SRMR(S = S, sigma.hat = Sigma)
    xxSRMR    <- SRMR(S = S_xx, sigma.hat = Sigma_xx)
    yySRMR    <- SRMR(S = S_yy, sigma.hat = Sigma_yy)
    yxSRMR    <- SRMR.xy(S = S, Sxy = S_xy, sigmaxy.hat = Sigma_xy, 
                         xnames = xnames, ynames = ynames)
  } else {
    fullSRMR <- xxSRMR <- yySRMR <- yxSRMR <- NULL
  }
  
  attr(Ypred, "fullSRMR") <- ifelse(!is.null(fullSRMR), fullSRMR, NA)
  attr(Ypred, "xxSRMR")   <- ifelse(!is.null(xxSRMR), xxSRMR, NA)
  attr(Ypred, "yySRMR")   <- ifelse(!is.null(yySRMR), yySRMR, NA)
  attr(Ypred, "yxSRMR")   <- ifelse(!is.null(yxSRMR), yxSRMR, NA)
  
  final <- list(Ypred = Ypred, Ytrue = Ytrue)
  
  return(final) 
  
}

# fit <- fitmod(dat = calibration)
# lav.predict.y(prediction, fit, gamma1 = 0.5, gamma2 = 0.3,
#           xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----

# prediction for the K partitions----
lav.predict.y.part <- function(dat, K, partid, 
                               gamma1, gamma2, equal.gammas, 
                               n_x, n_eta_x, n_y,  n_eta_y,
                               xnames, ynames) {
  partdat <- partition(partid = partid, dat = dat, K = K) # partitioned data
  
  # row and column names
  mat.rows <- do.call("c", 
                      lapply(1:K, function(k) 
                        paste0(k, ".", 1:nrow(partdat[[k]]$test))))
  if (!equal.gammas) {
    mat.cols <- apply(expand.grid(ynames, as.character(gamma1), as.character(gamma2)), 
                      1, paste0, collapse = ",") 
    # use as.character() above to paste only 0/1 instead of 0.0 and 1.0
    # because in the for loops, g1/g2 are 0/1 not 0.0/0.1
    sqdevmat <- matrix(NA, nrow(dat), length(gamma1)*length(gamma2)*length(ynames),
                       dimnames = list(mat.rows, mat.cols)) # matrix with squared deviations  
  } else {
    mat.cols <- apply(expand.grid(ynames, 
                                  paste0(as.character(gamma1), ",", as.character(gamma2))), 
                      1, paste0, collapse = ",")
    sqdevmat <- matrix(NA, nrow(dat), length(gamma1)*length(ynames),
                       dimnames = list(mat.rows, mat.cols)) # `ncol` computed this way because only equal gammas are considered
  }
  
  for (k in 1:K) {
    for (g1 in gamma1) {
      for (g2 in gamma2) {
        if (equal.gammas && g1 != g2) {
          next
        }
        fitpart <- fitmod(dat = partdat[[k]]$train,
                          n_x = n_x, n_eta_x = n_eta_x, 
                          n_y = n_y,  n_eta_y = n_eta_y) # fit to only the training data
        
        predpart <- lav.predict.y(preddat = partdat[[k]]$test,
                                  califit = fitpart,
                                  gamma1 = g1, gamma2 = g2, xnames = xnames,
                                  ynames = ynames)
        
        sqdevmat[rownames(sqdevmat)[grep(pattern = paste0("^", k, "\\."), 
                                         rownames(sqdevmat))], 
                 paste0(ynames, ",", g1,",", g2)] <- 
          as.matrix((predpart$Ypred - predpart$Ytrue)^2)
      }
    }
  } 
  return(sqdevmat)
}

# t0 <- Sys.time()
# foo <- lav.predict.y.part(dat = calibration, K = 10, partid = part.ids,
#                       gamma1 = seq(0,1,0.1), gamma2 = seq(0,1,0.1),
#                       equal.gammas = F, 
#                       n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1,
#                       xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1, t0, "sec")

#----

# compute RMSEp(r) for each gamma1,gamma2 combination and return gamma1,2 values with min(RMSEp)----
lav.predict.y.gamma <- function(dat, K, partid, 
                                gamma1, gamma2, equal.gammas,
                                n_x, n_eta_x, n_y,  n_eta_y,
                                xnames, ynames) {
  sqdevmat <- lav.predict.y.part(dat = dat, K = K, 
                                 partid = partid, gamma1 = gamma1, gamma2 = gamma2,
                                 equal.gammas = equal.gammas,
                                 n_x = n_x, n_eta_x = n_eta_x, 
                                 n_y = n_y,  n_eta_y = n_eta_y,
                                 xnames = xnames, ynames = ynames)
  
  RMSEp  <- if (!equal.gammas) {
    expand.grid(gamma1 = gamma1, gamma2 = gamma2, RMSEp = NA)
  } else {
    as.data.frame(cbind(gamma1 = gamma1, gamma2 = gamma2, RMSEp = NA))
  }
  
  # RMSEp <- matrix(NA, length(gamma1)*length(gamma2), 3,
  #                 dimnames = list(NULL, c("gamma1", "gamma2", "RMSEp")))
  
  # TODO only doing RMSEp for now, but if necessary, can add RMSEpr later
  
  for (k in 1:K) {
    for (g1 in gamma1) {
      for (g2 in gamma2) {
        if (equal.gammas && g1 != g2) {
          next
        }
        # compute single RMSE value per gamma1-gamma2 combination
        RMSEp.val <- 
          sum(sqdevmat[rownames(sqdevmat)[grep(pattern = paste0(k, "\\."), 
                                               rownames(sqdevmat))],
                       paste0(ynames, ",", g1,",", g2)])/(length(grep(pattern = paste0(k, "\\."), 
                                                                      rownames(sqdevmat)))*length(ynames))
        
        RMSEp[RMSEp$gamma1 == g1 & RMSEp$gamma2 == g2, "RMSEp"] <- RMSEp.val
      }
    }
  }
  
  min.RMSE.val <- RMSEp[which(RMSEp$RMSEp == min(RMSEp$RMSEp)),] # minimum RMSE value and associated gamma1/2
  
  return(list(gamma1 = min.RMSE.val$gamma1, gamma2 = min.RMSE.val$gamma2))  
}

# lav.predict.y.gamma(dat = calibration, K = 10, partid = part.ids,
#                 gamma1 = seq(0,1,0.1), gamma2 = seq(0,1,0.1),
#                 equal.gammas = F,
#                 n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1,
#                 xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))

#----

# prediction rule with cross-validation----
lav.predict.y.cv <- function(calidat, preddat, califit, CV,
                             gamma1, gamma2, equal.gammas, 
                             n_x, n_eta_x, n_y,  n_eta_y,
                             K, partid, 
                             xnames, ynames) {
  
  if (CV) {
    if (is.numeric(gamma1) && is.numeric(gamma2) && 
        length(gamma1) > 1L && length(gamma2) > 1L && 
        all(0L <= gamma1) && all(gamma1 <= 1L) && 
        all(0L <= gamma2) && all(gamma2 <= 1L)) {
      gamma.vals <- lav.predict.y.gamma(dat = calidat, K = K, partid = partid, 
                                        gamma1 = gamma1, gamma2 = gamma2,
                                        equal.gammas = equal.gammas,
                                        n_x = n_x, n_eta_x = n_eta_x,
                                        n_y = n_y, n_eta_y = n_eta_y,
                                        xnames = xnames, ynames = ynames)
    } else {
      stop("specify numeric vectors with length > 1L and only containing values between 
         0 and 1 for `gamma1` and `gamma2` when `CV = TRUE`")
    }
  } else {
    if (is.numeric(gamma1) && is.numeric(gamma2) && 
        length(gamma1) == 1L && length(gamma2) == 1L &&
        0L <= gamma1 && gamma1 <= 1L && 0L <= gamma2 && gamma2 <= 1L) {
      gamma.vals <- list(gamma1 = gamma1, gamma2 = gamma2)
    } else {
      stop("specify numeric values between 0 and 1 for `gamma1` and `gamma2` when `CV = FALSE`")
    }
  }
  
  Ypred <- lav.predict.y(preddat = preddat,
                         califit = califit, gamma1 = gamma.vals$gamma1,
                         gamma2 = gamma.vals$gamma2,
                         xnames = xnames, ynames = ynames, srmr = T)$Ypred
  colnames(Ypred) <- ynames
  
  attr(Ypred, "CV")     <- CV # was cross-validation performed?
  attr(Ypred, "gamma1") <- gamma.vals$gamma1 # for xx part
  attr(Ypred, "gamma2") <- gamma.vals$gamma2 # for xy part
  
  return(Ypred)
  
}

# t0 <- Sys.time()
# bar <- lav.predict.y.cv(calidat = calibration, preddat = prediction,
#                         califit = fit, CV = T,
#                     gamma1 = seq(0,1,0.1), gamma2 = seq(0,1,0.1), 
#                     equal.gammas = F,
#                     n_x = 27, n_eta_x = 3, n_y = 9, n_eta_y = 1,
#                     K = 10,
#                     partid = part.ids, xnames = paste0("x", 4:7), ynames = paste0("x", 1:3))
# t1 <- Sys.time()
# diff <- difftime(t1,t0,"sec")

#----
