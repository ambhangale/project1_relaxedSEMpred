## Aditi M. Bhangale
## Last updated: 7 January 2026

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
# relaxed SEM
### custom functions for postprocessing

## to run on ALICE
allRes <- do.call("c", lapply(dir(pattern = "res"), readRDS))

sum_result <- function(resList, save.files = T) {
  ## RMSEp
  RMSEp <- data.frame(do.call("rbind", lapply(1:length(resList),
                                   function(x) resList[[x]]$RMSEp)))
  
  makeNumeric <- c("sampID", "nCal", "nPred", "n_x", "n_eta_x", "n_y", "n_eta_y", 
                   "beta", "r" , "RMSEp", "runTime", "npar", "df", "CFI", 
                   "RMSEA", "RMSEA.lowCI", "RMSEA.upCI",
                   "lav.fullSRMR", "lav.xxSRMR", "lav.yySRMR", "lav.yxSRMR",
                   "lav.alpha1", "lav.alpha2", "en.alpha", "en.lambda")
  RMSEp[, makeNumeric] <- apply(RMSEp[, makeNumeric], 2, 
                                function(x) as.numeric(x))
  
  RMSEp$nCal <- factor(RMSEp$nCal)
  RMSEp$beta <- factor(RMSEp$beta, levels = c(0.3, 0.5), labels = c("weak", "strong"))
  RMSEp$r <- factor(RMSEp$r, levels = c(0.3, 0.7), labels = c("weak", "strong"))
  RMSEp$n_x <- factor(RMSEp$n_x)
  RMSEp$n_y <- factor(RMSEp$n_y)
  RMSEp$n_eta_x <- factor(RMSEp$n_eta_x)
  RMSEp$misspecify <- factor(RMSEp$misspecify, levels = c(FALSE, TRUE),
                             labels = c("no", "yes"))
  RMSEp$miss.part <- factor(RMSEp$miss.part, 
                            levels = c(NA, "xx:rescov", "xx:crossload", 
                                       "xy:direct", "both:cov", "both:load"), 
                            labels = c("none", "meas:rescov", "meas:crossload", 
                                       "struc:direct", "both:cov", "both:load"), 
                            exclude = NULL)
  RMSEp$miss.strength <- factor(RMSEp$miss.strength, 
                                 levels = c(NA, "weak", "strong"),
                                 labels = c("none", "weak", "strong"),
                                 exclude = NULL)
  RMSEp$mod.size <- paste0("n_eta_x=",RMSEp$n_eta_x, ",n_x=",RMSEp$n_x,
                           ",n_eta_y=",RMSEp$n_eta_y,",n_y=",RMSEp$n_y)
  RMSEp$mod.size <- factor(RMSEp$mod.size,
                           levels = c("n_eta_x=1,n_x=4,n_eta_y=1,n_y=1",
                                      "n_eta_x=1,n_x=4,n_eta_y=1,n_y=4",
                                      "n_eta_x=1,n_x=4,n_eta_y=1,n_y=8",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=1",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=4",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=8",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=1",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=4",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=8",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=1",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=4",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=8"))
  
  ## RMSEpr
  RMSEpr <- as.data.frame(do.call("rbind", 
                                  lapply(1:length(resList),
                                         function(x) resList[[x]]$RMSEpr)))
  makeNumeric <- gsub("RMSEp", "RMSEpr", makeNumeric)
  RMSEpr[, makeNumeric] <- apply(RMSEpr[, makeNumeric], 2, 
                                function(x) as.numeric(x))
  
  RMSEpr$nCal <- factor(RMSEpr$nCal)
  RMSEpr$beta <- factor(RMSEpr$beta, levels = c(0.3, 0.5), labels = c("weak", "strong"))
  RMSEpr$r <- factor(RMSEpr$r, levels = c(0.3, 0.7), labels = c("weak", "strong"))
  RMSEpr$n_x <- factor(RMSEpr$n_x)
  RMSEpr$n_y <- factor(RMSEpr$n_y)
  RMSEpr$n_eta_x <- factor(RMSEpr$n_eta_x)
  RMSEpr$misspecify <- factor(RMSEpr$misspecify, levels = c(FALSE, TRUE),
                             labels = c("no", "yes"))
  RMSEpr$miss.part <- factor(RMSEpr$miss.part, 
                            levels = c(NA, "xx:rescov", "xx:crossload", 
                                       "xy:direct", "both:cov", "both:load"), 
                            labels = c("none", "meas:rescov", "meas:crossload", 
                                       "struc:direct", "both:cov", "both:load"), 
                            exclude = NULL)
  RMSEpr$miss.strength <- factor(RMSEpr$miss.strength, 
                                 levels = c(NA, "weak", "strong"),
                                 labels = c("none", "weak", "strong"),
                                 exclude = NULL)
  RMSEpr$mod.size <- paste0("n_eta_x=",RMSEpr$n_eta_x, ",n_x=",RMSEpr$n_x,
                           ",n_eta_y=",RMSEpr$n_eta_y,",n_y=",RMSEpr$n_y)
  RMSEpr$mod.size <- factor(RMSEpr$mod.size,
                           levels = c("n_eta_x=1,n_x=4,n_eta_y=1,n_y=1",
                                      "n_eta_x=1,n_x=4,n_eta_y=1,n_y=4",
                                      "n_eta_x=1,n_x=4,n_eta_y=1,n_y=8",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=1",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=4",
                                      "n_eta_x=1,n_x=8,n_eta_y=1,n_y=8",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=1",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=4",
                                      "n_eta_x=3,n_x=12,n_eta_y=1,n_y=8",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=1",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=4",
                                      "n_eta_x=3,n_x=24,n_eta_y=1,n_y=8"))
  
  ## lavcv.alphas
  lavcv.alphas <- as.data.frame(do.call("rbind",
                                        lapply(1:length(resList),
                                               function(x) resList[[x]]$lavcv.alphas)))
  lavcv.alphas[, c(makeNumeric[1:9], paste0("alpha", 1:2))] <- apply(lavcv.alphas[, c(makeNumeric[1:9], 
                                                                                    paste0("alpha", 1:2))],
                                                                    2, function(x) as.numeric(x))
  
  if (save.files) {
    saveRDS(RMSEp, file = paste0("RMSEp_relaxSEM", Sys.Date(), ".rds"))
    saveRDS(RMSEpr, file = paste0("RMSEpr_relaxSEM", Sys.Date(), ".rds"))
    saveRDS(lavcv.alphas, file = paste0("lavcv_alphas_relaxSEM", Sys.Date(), ".rds"))
  }
  
  return(list(RMSEp = RMSEp, RMSEpr = RMSEpr, lavcv.alphas = lavcv.alphas))
}

sum_result(allRes) # source file on ALICE to run all at once
