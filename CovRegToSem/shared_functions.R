library(MASS)
library(lavaan)
# Function to simulate data from a one-factor model with fixed loadings
simulate_one_factor_cov <- function(n_vars, loading = 0.7, error_var = 1, add_resid_cov = FALSE, resid_pairs = list(c(1, 2), c(3, 4)), resid_cov_value = 0.45) {
  lambda <- matrix(runif(n = n_vars, min = .3, max = .7), ncol = 1)
  implied_cov <- lambda %*% t(lambda) + diag(runif(n_vars, min = 0.5, max = 1))
  
  if (add_resid_cov) {
    for (pair in resid_pairs) {
      i <- pair[1]
      j <- pair[2]
      implied_cov[i, j] <- implied_cov[j, i] <- resid_cov_value
    }
  }
  
  colnames(implied_cov) <- rownames(implied_cov) <- paste0("x", 1:n_vars)
  implied_cov
  # browser()
}



# Function to estimate model-implied covariance matrix using lavaan with equal loadings
estimate_one_factor_cov <- function(data) {
  var_names <- colnames(data)
  model_string <- paste0(
    "F1 =~ ", paste0("l*", var_names, collapse = " + "), "\n",
    paste0(var_names, " ~~ v*", var_names, collapse = "\n")
  )
  fit <- cfa(model_string, data = data, fixed.x = FALSE)
  fitted(fit)$cov
}

squared_frobenius_norm <- function(A) {
  sum(A^2)
}

random_cov_matrix <- function(p) {
  set.seed(12356)
  A <- matrix(rnorm(p^2), nrow = p)
  cov_mat <- crossprod(A)
  colnames(cov_mat) <- rownames(cov_mat) <- paste0("x", 1:p)
  cov_mat
}


