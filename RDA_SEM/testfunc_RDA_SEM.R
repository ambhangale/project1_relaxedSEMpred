## Aditi M. Bhangale
## Last updated: 19 February 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(lavaan) 
library(MASS) # for mvrnorm()

# n = 100; regSxy = F; alpha1 = 0.5; misspecify = F
# xnames = c(paste0("x", 1:3), paste0("y", 1:4))
# ynames = paste0("y", 5:8)

# 3 (a)
# Y <- mu.Y + t(S_xy)*solve((1-alpha)*Sigma_xx + alpha*S_xx)*(x0 - mu.X)

# 3 (b)
# Y <- mu.Y + ((1-alpha)*t(Sigma_xy) + alpha*t(S_xy))*solve((1-alpha)*Sigma_xx + alpha*S_xx)*(x0 - mu.X)


data("PoliticalDemocracy")

mod <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4 
     dem65 =~ y5 + y6 + y7 + y8
    
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
' # empirical example in De Rooij et al. (2023)

fit <- sem(mod, data = PoliticalDemocracy, meanstructure = T) # mean structure is saturated
# summary(fit, fit.measures = T, std = T)
popStats <- lavInspect(fit, "implied") # use to generate new data

# generate training and test datasets
# FIXME, can also put in `testrule()` function with `n` argument if necessary. but for now, this works
set.seed(123)
traindat <- mvrnorm(n = 10000, mu = popStats$mean, Sigma = popStats$cov)
testdat <- mvrnorm(n = 10000, mu = popStats$mean, Sigma = popStats$cov) 

fit <- sem(mod, traindat, meanstructure = T) # for `lavpredictY()`

testrule <- function(train, test, regSxy, alpha1, alpha2 = NULL, 
                     xnames = c(paste0("x", 1:3), paste0("y", 1:4)), 
                     ynames = paste0("y", 5:8), misspecify = F) {
  
  # sample statistics of the training dataset
  S <- cov(train) #FIXME nrow(train)
  S_xx <- S[xnames, xnames]
  S_xy <- S[xnames, ynames]
  
  # values from test dataset
  X0 <- test[, xnames] # to be inputted in formulae
  Ytest <- test[, ynames] # original Y values from the test dataset
  
  # model
  mod <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4 
     dem65 =~ y5 + y6 + y7 + y8
    
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
    
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
' 

  if (misspecify) mod <- paste(mod, 'dem60 ~ x1', collapse = "\n") # add additional direct effect
    
  fit <- sem(mod, data = train, meanstructure = T)
  ImpliedStats <- lavInspect(fit, "implied")
  Sigma_xx <- ImpliedStats$cov[xnames, xnames]
  Sigma_xy <- ImpliedStats$cov[xnames, ynames]
  Mu_x <- ImpliedStats$mean[xnames]
  Mu_y <- ImpliedStats$mean[ynames]
  
  if (!regSxy) {
    Ypred <- t(Mu_y + t(S_xy) %*% solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% (t(X0) - Mu_x)) # t(X0) to make compatible with Mu_x and t(Mu_y + ...) to change to long format. also below.
  } else {
    if(!is.null(alpha2)) {
    Ypred <- t(Mu_y + ((1-alpha2)*t(Sigma_xy) + alpha2*t(S_xy)) %*% solve((1-alpha1)*Sigma_xx + alpha1*S_xx) %*% (t(X0) - Mu_x))
    } else {
    stop("specify value for `alpha2`")
    }
  }
  
  return(list(Ypred = Ypred, Ytest = Ytest, bias = Ypred - Ytest, 
              final = as.data.frame(cbind(regSxy = regSxy, alpha1 = alpha1, 
                                          alpha2 = ifelse(!is.null(alpha2), alpha2, NA),
                                          misspecify = misspecify, 
                                          meanBias = colMeans(Ypred - Ytest),
                                          RMSE = sqrt(mean((Ypred - Ytest)^2)),
                                          yname = ynames)))) 
  # save `final` in long format for `ggplot()`
}


# test the function
# predY.A <- testrule(train = traindat, test = testdat, regSxy = F, 
#                     alpha1 = 0.3)
# predY.A <- testrule(train = traindat, test = testdat, regSxy = T, 
#                     alpha1 = 0.3, alpha2 = 0.2)

predY.og <- lavPredictY(fit, newdata = testdat, 
                   xnames = c(paste0("x", 1:3), paste0("y", 1:4)), 
                   ynames = paste0("y", 5:8)) 
## predY.A with alpha1 = 0 and predY.og will not match exactly as `lavPredictY()` uses 
## Sigma_xy in the formula, whereas `testrule()` uses S_xy in the formula

conds <- rbind(expand.grid(regSxy = F, alpha1 = seq(0,1,0.1), alpha2 = NA, misspecify = c(F,T)),
                expand.grid(regSxy = T, alpha1 = seq(0,1,0.1), 
                            alpha2 = seq(0,1,0.1), misspecify = c(F,T)))
# dim(conds) should be (11*2) + (11*11*2) = 264

# results
t0 <- Sys.time()
res_list <- mapply(testrule, regSxy = conds$regSxy, alpha1 = conds$alpha1, 
              alpha2 = conds$alpha2, misspecify = conds$misspecify, 
              MoreArgs = list(train = traindat, test = testdat), SIMPLIFY = F)
t1 <- Sys.time()
diff <- difftime(t1, t0, units = "sec")

# save only final
res_final <- as.data.frame(do.call("rbind", 
                     lapply(1:length(res_list), function(i) res_list[[i]]$final)))
res_final$meanBias <- as.numeric(res_final$meanBias) # convert to numeric
res_final$RMSE <- as.numeric(res_final$RMSE) # convert to numeric


library(ggplot2)


## mean 
# for regSxy == F
ggplot(data = res_final[res_final$regSxy == F & is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha1, y = meanBias)) +
  geom_point() + facet_wrap(~ yname + misspecify, ncol = 2) + geom_hline(yintercept = 0) 

# for regSxy == T & misspecify == F
ggplot(data = res_final[res_final$regSxy == T & res_final$misspecify == F & !is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha2, y = meanBias)) +
  geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) + geom_hline(yintercept = 0) + ggtitle("regSxy = T and misspecify = F")

# for regSxy == T & misspecify == T
ggplot(data = res_final[res_final$regSxy == T & res_final$misspecify == T & !is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha2, y = meanBias)) +
  geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) + geom_hline(yintercept = 0) + ggtitle("regSxy = T and misspecify = F") 


## RMSE 
# for regSxy == F
ggplot(data = res_final[res_final$regSxy == F & is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha1, y = RMSE)) +
  geom_point() + facet_wrap(~ yname + misspecify, ncol = 2) 

# for regSxy == T & misspecify == F
ggplot(data = res_final[res_final$regSxy == T & res_final$misspecify == F & !is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha2, y = RMSE)) +
  geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) + ggtitle("regSxy = T and misspecify = F")

# for regSxy == T & misspecify == T
ggplot(data = res_final[res_final$regSxy == T & res_final$misspecify == T & !is.na(res_final$alpha2), ], 
       mapping = aes(x = alpha2, y = RMSE)) +
  geom_point() + facet_wrap(~ yname + alpha1, ncol = 11) + ggtitle("regSxy = T and misspecify = T") 


testrule_CV <- function() {
  # TODO maybe try (nested) CV for the above formulae
  
  # TODO add MSE as an element in the returned object
}




