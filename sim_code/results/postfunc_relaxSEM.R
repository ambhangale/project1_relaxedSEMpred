## Aditi M. Bhangale
## Last updated: 9 May 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### custom functions for postprocessing

sum_result <- function(resList) {
  RMSEp <- data.frame(do.call("rbind", lapply(1:length(resList),
                                   function(x) resList[[x]]$RMSEp)))
  
  makeNumeric <- c("RMSEp", "runTime", "RMSEA", "RMSEA.lowCI", "RMSEA.upCI",
                   "lav.alpha1", "lav.alpha2", "en.alpha", "en.lambda")
  RMSEp[, makeNumeric] <- apply(RMSEp[, makeNumeric], 2, 
                                function(x) as.numeric(x))
  
  RMSEpr <- as.data.frame(do.call("rbind", 
                                  lapply(1:length(resList),
                                         function(x) resList[[x]]$RMSEpr)))
  makeNumeric <- gsub("RMSEp", "RMSEpr", makeNumeric)
  RMSEpr[, makeNumeric] <- apply(RMSEpr[, makeNumeric], 2, 
                                function(x) as.numeric(x))
  
  return(list(RMSEp = RMSEp, RMSEpr = RMSEpr))
}
