## Aditi M. Bhangale
## Last updated: 12 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### postprocessing file

library(here)
library(ggplot2)
source(here("sim_code/results", "postfunc_relaxSEM.R"))

result <- readRDS("PD_resList_2025-05-08.rds")
## the RMSEpr values in this list are incorrect because you found a bug in your code
## right after you ran this. but it doesn't matter, because you are only interested
## in RMSEp values in this case anyways

RMSEp <- sum_result(resList = result)$RMSEp

## boxplots
ggplot(RMSEp, aes(x = method, y = RMSEp, fill = factor(method))) + 
  geom_boxplot(aes(group = factor(method))) + 
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(nCal ~ misspecify) +
  xlab("method") + 
  ylab("RMSEp") 

## violin plots
ggplot(RMSEp, aes(x = method, y = RMSEp, fill = factor(method))) + 
  geom_violin(aes(group = factor(method))) + 
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(nCal ~ misspecify) +
  xlab("method") + 
  ylab("RMSEp") 

## error bars
aggregate2 <- function(formula, data, FUN1 = "median", FUN2 = "min", FUN3 = "max") {
  A <- aggregate(x = as.formula(formula), data = data, FUN = FUN1)
  names(A)[!names(A) %in% c("method", "nCal", "misspecify")] <- 
    paste0(FUN1, "_",  names(A)[!names(A) %in% c("method", "nCal", "misspecify")])
  B <- aggregate(x = as.formula(formula), data = data, FUN = FUN2)
  names(B)[!names(B) %in% c("method", "nCal", "misspecify")] <- 
    paste0(FUN2, "_",  names(B)[!names(B) %in% c("method", "nCal", "misspecify")])
  C <- aggregate(x = as.formula(formula), data = data, FUN = FUN3)
  names(C)[!names(C) %in% c("method", "nCal", "misspecify")] <- 
    paste0(FUN3, "_",  names(C)[!names(C) %in% c("method", "nCal", "misspecify")])
  
  range <- Reduce("merge", list(A,B,C))
  
  return(range)
}

range.RMSEp <- aggregate2(formula = "RMSEp ~ method + nCal + misspecify", 
                          data = RMSEp) # min, median, and max RMSEp

ggplot(range.RMSEp, mapping = aes(x = method, y = median_RMSEp, fill = factor(method))) +
  geom_errorbar(aes(group = factor(method), ymin = min_RMSEp, ymax = max_RMSEp)) +
  geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) + 
  facet_grid(nCal ~ misspecify) +
  xlab("method") + 
  ylab("RMSEp") 


