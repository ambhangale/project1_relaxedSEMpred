## Aditi M. Bhangale
## Last updated: 24 March 2025

# Creating a function that applies the RDA-like constraints on the SEM prediction rule
## SEM on political democracy dataset
### custom functions for result summarisation and plots

# getwd()
# setwd("/Users/Aditi_2/Desktop/Universiteit Leiden/Projects/project_1_relaxedSEMpred/RDA_SEM")

library(ggplot2)
# library(gghighlight)

# summarise results
sum_result <- function(resList) {
  RMSEpr <- as.data.frame(do.call("rbind",
                                  lapply(1:length(resList), function(x) resList[[x]]$RMSEpr.result)))
  RMSEpr$meanBias <- as.numeric(RMSEpr$meanBias)
  RMSEpr$RMSEpr   <- as.numeric(RMSEpr$RMSEpr)
  
  RMSEp <- as.data.frame(do.call("rbind",
                                 lapply(1:length(resList), function(x) resList[[x]]$RMSEp.result)))
  RMSEp$RMSEp <- as.numeric(RMSEp$RMSEp)
  
  return(list(RMSEpr = RMSEpr, RMSEp = RMSEp))
}

# plot results
plot_result <- function(data, plot.stat, plot.title = NULL, facet_cols = 11) { 
  plot <- ggplot(data = data, mapping = aes(x = alpha1, y = .data[[plot.stat]])) +
    geom_point() + 
    facet_wrap(as.formula(paste("~", ifelse(plot.stat == "RMSEpr", "yname +", ""), "alpha2")), ncol = facet_cols)
    ggtitle(plot.title)
    
    if (plot.stat == "RMSEp") {
      min.val <- subset(data, data$RMSEp == min(data$RMSEp))
      plot <- plot + geom_point(data = min.val, color = "red") # highlight minimum value
    }
  
  return(plot)
}

# plot_result(RMSE.cs.1e2$RMSEpr, plot.stat = "RMSEpr", plot.title = "test") # test
