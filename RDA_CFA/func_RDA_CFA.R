## Aditi M. Bhangale
## Last updated: 20 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## CFA example
### function(s) file

library(here)
source(here("RDA_CFA", "datagen_RDA_CFA.R"))

fitmod <- function(train) {
  mod <- '
# factor loadings
F1 =~ x1 + x2 + x3
F2 =~ x4 + x5 + x6 + x7

# factor (co)variances
F1 ~~ 1*F1 + F2
F2 ~~ 1*F2

# item (co)variances
x1 ~~ x1
x2 ~~ x2
x3 ~~ x3
x4 ~~ x4
x5 ~~ x5
x6 ~~ x6
x7 ~~ x7

# factor means
F1 ~ 0*1
F2 ~ 0*1

# item means
x1 ~ 1
x2 ~ 1
x3 ~ 1
x4 ~ 1
x5 ~ 1
x6 ~ 1
x7 ~ 1
'
  
  fit <- lavaan(mod, data = train, meanstructure = T)
  
  return(fit)
}

testrule <- function() {
  
}
