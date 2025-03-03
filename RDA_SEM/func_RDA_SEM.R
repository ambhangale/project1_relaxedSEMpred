## Aditi M. Bhangale
## Last updated: 26 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## function(s) file

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

source("datagen_RDA_SEM.R")


# function to fit model
fitmod <- function(train, mod = pol_model) {
  fit <- sem(mod, data = train, meanstructure = T) # mean structure saturated for now
  return(fit)
}


# function to apply formulae
testrule <- function(ntrain, ntest, misspecify, 
                     regXY, XYtype = NULL, alpha1, alpha2 = NULL,
                     xnames = c(paste0("x", 1:3), paste0("y", 1:4)), 
                     ynames = paste0("y", 5:8)) {
  
  dat <- gendat(ntrain = ntrain, ntest = ntest, misspecify = misspecify)
  train <- dat$train
  test <- dat$test
  
  S <- (cov(train)*(nrow(train)-1)) / nrow(train) 
  # `cov()` computes covmat using denominator 'n-1', but `lavaan` uses 'n' to scale
  # covmats and SEs. thus, rescaled the covmat produced by `cov()` so that it scales by 'n'
  S_xx <- S[xnames, xnames]
  S_xy <- S[xnames, ynames] #FIXME eventually use _yx to avoid the t()?
  
  # values from test dataset
  X0 <- test[,xnames] # to be inputted into formulae
  Ytest <- test[,ynames] # original Y values from test dataset
  
  fit <- fitmod(train = train)
  ImpliedStats <- lavInspect(fit, "implied")
  Sigma_xx     <- ImpliedStats$cov[xnames, xnames]
  Sigma_xy     <- ImpliedStats$cov[xnames, ynames] #FIXME eventually use _yx to avoid the t()?
  Mu_x         <- ImpliedStats$mean[xnames]
  Mu_y         <- ImpliedStats$mean[ynames]
  
  if(!regXY) {
    if (XYtype == "S.xy") {
      if (0L <= alpha1 && alpha1 <= 1L) {
      Ypred <- t(Mu_y + t(S_xy) %*% 
                   solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                   (t(X0) - Mu_x)) 
      # t(X0) to make compatible with Mu_x and t(Mu_y + ...) to change to long format. also below.
      } else {
        stop("specify value between 0 and 1 for `alpha1`")
      }
    } else if (XYtype == "Sigma.xy") { # if alpha = 0, result will match De Rooij et al. (2022) prediction rule
      if (0L <= alpha1 && alpha1 <= 1L) {
      Ypred <- t(Mu_y + t(Sigma_xy) %*% 
                   solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                   (t(X0) - Mu_x)) 
      } else {
        stop("specify value between 0 and 1 for `alpha1`")
      }
    } else {
      stop("specify valid type for XY covariance matrix")
    }
  } else {
    if(!is.null(alpha2)) {
      if (0L <= alpha1 && alpha1 <= 1L && 0L <= alpha2 && alpha2 <= 1L) {
      Ypred <- t(Mu_y + ((1-alpha2)*t(Sigma_xy) + alpha2*t(S_xy)) %*% 
                   solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% 
                   (t(X0) - Mu_x))
      } else {
        stop("specify values between 0 and 1 for `alpha1` and `alpha2`")
      }
    } else {
      stop("specify value for `alpha2`")
    }
  }
  
  bias <- Ypred - Ytest
  
  RMSEpr.result <- as.data.frame(cbind(regXY = regXY, 
                                       XYtype = ifelse(!is.null(XYtype), XYtype, NA),
                                       alpha1 = alpha1, 
                                       alpha2 = ifelse(!is.null(alpha2), alpha2, NA),
                                       misspecify = misspecify,
                                       meanBias = colMeans(bias),
                                       RMSEpr = sqrt(colMeans((bias)^2)),
                                       yname = ynames))
  
  RMSEp.result <- as.data.frame(cbind(regXY = regXY, 
                                      XYtype = ifelse(!is.null(XYtype), XYtype, NA),
                                      alpha1 = alpha1, 
                                      alpha2 = ifelse(!is.null(alpha2), alpha2, NA),
                                      misspecify = misspecify,
                                      RMSEp = sqrt(sum((bias)^2)/(length(ynames)*ntest))))
  
  final <- list(Ypred = Ypred, Ytest = Ytest, bias = bias,
                RMSEpr.result = RMSEpr.result, RMSEp.result = RMSEp.result)
  # save all arguments as attributes, just in case we need them later
  attr(final, "ntrain")     <- ntrain
  attr(final, "ntest")      <- ntest
  attr(final, "misspecify") <- misspecify
  attr(final, "regXY")      <- regXY
  attr(final, "XYtype")     <- ifelse(!is.null(XYtype), XYtype, "NA")
  attr(final, "alpha1")     <- alpha1
  attr(final, "alpha2")     <- ifelse(!is.null(alpha2), alpha2, "NA")
  attr(final, "xnames")     <- xnames
  attr(final, "ynames")     <- ynames
  
  return(final)
  
}

# test `testrule()`
# foo <- testrule(ntrain = 100, ntest = 250, misspecify = F, regXY = F,
#                 XYtype = "Sigma.xy", alpha1 = 0)
# bar <- testrule(ntrain = 100, ntest = 250, misspecify = T, regXY = T, 
#                 alpha1 = 0.4, alpha2 = 0.3)



#TODO adapt function so that you can run `lavPredictY()` on the `fitmod()` output?
# or an alternative is to just use the `gendat()` function and fit the model yourself
# since it only needs to happen once

