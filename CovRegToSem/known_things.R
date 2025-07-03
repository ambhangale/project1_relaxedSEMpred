source("shared_functions.R")
repsLambda <- 10^2
repsMSE <- 10^2 
estimate <- TRUE
p <- 5
n <- 10^2
Sigma <- simulate_one_factor_cov(p, add_resid_cov = FALSE)
# Sigma <- random_cov_matrix(p)

clip01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

if(estimate){
  Sigma <- Sigma
}else{
  data <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  Sigma <- cov(data)
}


# estimate lambda via a simulation study that evaluates Equation 6 of https://www.sciencedirect.com/science/article/pii/S0167947314003107#s000015
dif_S_Sig <- dif_T_S <- cov_S_T <- dif_S_T_est <-  numeric(repsLambda)
for(i in 1:repsLambda){
  data <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  TMat <- estimate_one_factor_cov(data)
  SMat <- cov(data)
  
  dif_S_Sig[i] <- squared_frobenius_norm(SMat - Sigma) #
  dif_T_S[i] <- squared_frobenius_norm(SMat - TMat) #
  cov_S_T[i] <- sum(diag((SMat - Sigma) %*% (Sigma - TMat))) #
  if(i %% 100 == 0) {
    cat("Iteration:", i, "\n")
  }
}

lambda <- clip01((mean(dif_S_Sig) + mean(cov_S_T) / p) / mean(dif_T_S))

# compare the three estimators with found lambda
dif_S_Sig <- dif_SStar_Sig <- dif_T_Sig <- numeric(repsMSE)
for(i in 1:repsMSE){
  data <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  SMat <- cov(data)
  TMat <- estimate_one_factor_cov(data)
  SStar <- lambda * TMat + (1 - lambda) * SMat
  
  
  dif_S_Sig[i] <- squared_frobenius_norm(SMat - Sigma) 
  dif_T_Sig[i] <- squared_frobenius_norm(Sigma - TMat)
  dif_SStar_Sig[i] <- squared_frobenius_norm(SStar - Sigma) 
}

cat("Estimated?:" , estimate, "\n")
cat("lambda:", round(lambda, 4), "\n")
cat("MSE of the Sample Covariance Matrix:", round(mean(dif_S_Sig), 4), "\n")
cat("MSE of the SEM Covariance Matrix:", round(mean(dif_T_Sig), 4), "\n")
cat("MSE of the Compromise Covariance Matrix:", round(mean(dif_SStar_Sig), 4), "\n")


# for(i in 1:reps){
#   data <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
#   TMat <- estimate_one_factor_cov(data)
#   SMat <- cov(data)
#   
#   dif_frob_est[i] <- frob_risk_estimator(data)
#   dif_S_Sig[i] <- squared_frobenius_norm(SMat - Sigma)
#   dif_S_T_est[i] <- squared_frobenius_norm(SMat-TMat)
#   dif_SStar_Sig[i] <- squared_frobenius_norm(SStar - Sigma)
#   dif_T_S[i] <- squared_frobenius_norm(SMat - TMat)
#   dif_T_Sig[i] <- squared_frobenius_norm(Sigma - TMat)
#   cov_S_T[i] <- sum((diag((SMat - Sigma) %*% (Sigma - TMat))))
#   if(i %% 100 == 0) {
#     cat("Iteration:", i, "\n")
#   }
# }
# 
# 
# MSE_S <- mean(dif_S_Sig)
# print(MSE_S)
# formula_MSE_S <- 1/(n-1) * (m_trace(Sigma %*% Sigma) + m_trace(Sigma)^2) #seems correct
# print(formula_MSE_S)
# MSE_T <- mean(dif_T_Sig)
# MSE_SStar <- mean(dif_SStar_Sig)
# 
# 
# 
# shrinkage_lambda_boot <- function(X, T_fun, B = 1000, seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#   X   <- as.matrix(X); n <- nrow(X); p <- ncol(X)
#   S   <- cov(X)
#   
# 
#   dif_S_Sig <- dif_T_S <- cov_S_T <- numeric(B)
#   for (b in 1:B) {
#     data <- mvrnorm(n = n, mu = rep(0, p), Sigma = S)
#     TMat  <- T_fun(data)
#     SMat  <- cov(data)
#     
#     dif_S_Sig[b] <- squared_frobenius_norm(SMat - S)
#     dif_T_S[b] <- squared_frobenius_norm(SMat - TMat)
#     cov_S_T[b] <- sum((diag((SMat - S) %*% (S - TMat))))
#     
#     if(b %% 100 == 0) {
#       cat("Iteration:", b, "\n")
#     }
#     
#   }
#   lambda <- clip01((MSE_T_est + mean(cov_S_T) / p) / mean(dif_T_S))
# }
# MSE_T_est <- frob_risk_estimator(data)
# lambda_boot <- shrinkage_lambda_boot(data, estimate_one_factor_cov)
#   