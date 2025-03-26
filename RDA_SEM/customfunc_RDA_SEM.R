## Aditi M. Bhangale
## Last updated: 26 March 2025

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
plot_result <- function(data, plot.stat, plot.title = NULL, facet_cols = 11, ynames = NULL) { 
  plot <- ggplot(data = data, mapping = aes(x = alpha1, y = .data[[plot.stat]])) +
    geom_point() + 
    facet_wrap(as.formula(paste("~", ifelse(plot.stat == "RMSEpr", "yname +", ""), "alpha2")), 
               ncol = facet_cols) + theme(axis.text.x = element_text(angle = 45, size = 7)) + 
    ggtitle(plot.title)  # base plot
    
    if (plot.stat == "RMSEp") {
      min.val   <- subset(data, data$RMSEp == min(data$RMSEp)) # minimum value
      alpha.val <- do.call("rbind", mapply(function(x,y) subset(data, data$alpha1 == x & 
                                                                  data$alpha2 == y),
                                           x = c(0,0.5,1), y = c(0,0.5,1), SIMPLIFY = F)) # SEM & OLS pred, perfect compromise
      
      plot <- plot + geom_point(data = alpha.val, aes(color = c("SEM rule", "Compromise", "OLS rule"))) +
        geom_point(data = min.val, aes(color = "Minimum")) +
        guides(color = guide_legend(title = "RMSE value"))
      
    } else if (plot.stat == "RMSEpr") {
      if (!is.null(ynames)) {
        min.val   <- do.call("rbind", lapply(ynames, function(x) 
          subset(data[data$yname == x, ], 
                 data[data$yname == x, ]$RMSEpr == min(data[data$yname == x, ]$RMSEpr)))) # minimum value
        alpha.val <- do.call("rbind", lapply(ynames, 
                                             function(i)  do.call("rbind", mapply(function(x,y) 
                                               subset(data[data$yname == i, ], 
                                                      data[data$yname == i, ]$alpha1 == x & 
                                                        data[data$yname == i, ]$alpha2 == y),
                                               x = c(0, 0.5, 1), y = c(0, 0.5, 1), SIMPLIFY = F)))) # SEM & OLS pred, perfect compromise
        
        plot <- plot + geom_point(data = alpha.val, aes(color = rep(c("SEM rule", "Compromise", "OLS rule"), length(ynames)))) +
          geom_point(data = min.val, aes(color = "Minimum")) +
          guides(color = guide_legend(title = "RMSE value"))
        
        }
    }
  attr(plot, "min.val")   <- min.val
  attr(plot, "alpha.val") <- alpha.val 
  
  return(plot)
}

# plot_result(RMSE.cs.1e2$RMSEpr, plot.stat = "RMSEpr", plot.title = "test", ynames = paste0("x",1:3)) # test
