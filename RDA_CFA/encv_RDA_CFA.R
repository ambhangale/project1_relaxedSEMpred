## Aditi M. Bhangale
## Last updated: 22 April 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### function to predict using elastic net regression

library(here)
library(glmnet) # to perform elastic net regression
source(here("RDA_CFA", "part_RDA_CFA.R")) 

# prediction with an elastic net regresison model----
en.predict.y <- function(sampID, nCal, nPred, misspecify, 
                         alphas = seq(0,1,0.1), 
                         xnames = paste0("x", 4:7), ynames = paste0("x", 1:3),
                         seed = 10824) {
  t0 <- Sys.time()
  dat <- gendat(sampID = sampID, nCal = nCal, nPred = nPred, 
                misspecify = misspecify, seed = seed)
  calibration <- dat$calibration # calibration set
  prediction  <- dat$prediction # prediction set
  
  cv.errors <- rep(NA, length(alphas)) 
  
  # FIXME if i want to use same partitions for CV in my rule + en 
  # TODO can create a separate `partidx()` function that only comes up with partition assignments
  # TODO can replace redundant code in `partition()` function if i do the above
  # TODO can just include `foldid = (output of partidx() function)` in `cv.glmnet()`
  
  for (a in 1:length(alphas)) { # returns the optimal elastic net mixing parameter
    # save the minimum cross-validated error for each value in `alphas`
    cv.errors[a] <- min(cv.glmnet(x = dat$calibration[, xnames],
                                  y = dat$calibration[, ynames], 
                                  family = "mgaussian", alpha = alphas[a])$cvm) 
    # above, use default lambda sequence because that's what Mark did in the original simulation 
  }
  
  min.alpha <- alphas[which.min(cv.errors)] # alpha value with minimum cross-validated error
  
  lambda <- cv.glmnet(x = dat$calibration[, xnames],
                      y = dat$calibration[, ynames], 
                      family = "mgaussian", 
                      alpha = min.alpha)$lambda.min # minimum/optimal tuning parameter (lambda) value
  
  # FIXME: something to figure out: are two different fold patterns used for the two `cv.glmnet()` calls here?
  
  out <- glmnet(x = dat$calibration[, xnames],
                y = dat$calibration[, ynames],
                family = "mgaussian", alpha = min.alpha) # final model to use for predicting new values
  
  Ypred <- as.matrix(predict(out, newx = dat$prediction[,xnames], s = lambda)[,,1])
  
  Ytrue <- dat$prediction[,ynames] # true Y values
  
  bias <- Ypred - Ytrue # bias
  
  RMSEpr.result <- as.data.frame(cbind(sampID   = sampID, 
                                       nCal     = nCal, 
                                       nPred    = nPred, 
                                       alpha    = min.alpha, 
                                       lambda   = lambda,
                                       meanBias = colMeans(bias),
                                       RMSEpr   = sqrt(colMeans((bias)^2)),
                                       yname    = ynames))
  
  RMSEp.result <- as.data.frame(cbind(sampID     = sampID, 
                                      nCal       = nCal, 
                                      nPred      = nPred, 
                                      alpha      = min.alpha, 
                                      lambda     = lambda,
                                      misspecify = misspecify,
                                      RMSEp      = sqrt(sum((bias)^2)/(length(ynames)*nPred))))
  
  final <- list(Ypred = Ypred, Ytrue = Ytrue, bias = bias,
                RMSEpr.result = RMSEpr.result, RMSEp.result = RMSEp.result)
  
  t1   <- Sys.time()
  diff <- difftime(t1, t0, "sec")
  
  # save all arguments as attributes, just in case we need them later
  attr(final, "sampID")       <- sampID
  attr(final, "nCal")         <- nCal
  attr(final, "nPred")        <- nPred
  attr(final, "misspecify")   <- misspecify
  attr(final, "alpha")        <- min.alpha
  attr(final, "lambda")       <- lambda
  attr(final, "xnames")       <- xnames
  attr(final, "ynames")       <- ynames
  attr(final, "seed")         <- seed
  attr(final, "runtime")      <- diff
  
  return(final)
}

# en.predict.y(sampID = 1, nCal = 250, nPred = 250, misspecify = F)

#----
