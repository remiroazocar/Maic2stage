## This file processes the results of the simulation study and computes
## and plots the relevant performance measures.

rm(list=ls())

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/Maic2stage")
load(file="sim_settings.RData") # load simulation settings
source("functions.R") # load functions to compute performance measures

# packages required for plotting simulation study results
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(ggridges)) {install.packages("ggridges"); library(ggridges)}
if(!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}

scenario.settings <- pc # parameter combinations
scenarios <- nrow(pc) # number of scenarios
replicates <- N_sim # number of simulated data replicates per scenario
Delta.AB <- 0 # true value of marginal A vs. B treatment effect in S=2 is zero

# data frame to store performance metrics (4 methods, each row repeated 4 times)
simulation.metrics <- scenario.settings[rep(seq_len(scenarios), each=4), ]
metrics.names <- c("Method", "Bias", "Bias.MCSE", "LCI", "LCI.MCSE", "UCI", "UCI.MCSE",
                   "Cov", "Cov.MCSE", "ESE", "ESE.MCSE", "MSE", "MSE.MCSE", "number.NAs")
simulation.metrics[metrics.names] <- NA

# data frame to store all A vs. B marginal treatment effect point estimates
ate.table <- scenario.settings[rep(seq_len(scenarios), each=4*replicates),]
ate.table["Method"] <- NA
ate.table["ATE"] <- NA

# function that computes performance metrics for a given method
process.metrics <- function(means, variances, truth) {
  # remove NAs (correspond to no feasible weighting solution, empirical non-overlap)
  NAs <- is.na(means)
  means <- means[!NAs]  
  variances <- variances[!NAs]
  number.NAs <- sum(NAs) # number that do not converge (no feasible weighting solution)
  replicates <- length(means)
  bias.metric <- bias(means, truth)
  bias.metric.mcse <- bias.mcse(means)
  mse.metric <- mse(means, truth) 
  mse.metric.mcse <- mse.mcse(means, truth) 
  # construct Wald-type interval estimates using normal distribution
  lci <- means + qnorm(0.025)*sqrt(variances)
  uci <- means + qnorm(0.975)*sqrt(variances)
  lci.mean <- mean(lci)
  lci.mcse <- mcse.estimate(lci)
  uci.mean <- mean(uci)
  uci.mcse <- mcse.estimate(uci)
  cov <- coverage(lci, uci, truth)
  cov.mcse <- coverage.mcse(cov, replicates)
  empse.metric <- empse(means)
  empse.metric.mcse <- empse.mcse(empse.metric, replicates)
  list(bias.metric, bias.metric.mcse, lci.mean, lci.mcse, uci.mean, uci.mcse,
       cov, cov.mcse, empse.metric, empse.metric.mcse,
       mse.metric, mse.metric.mcse, number.NAs)
} 

j <- 1 # row counter for simulation metrics
k <- 1 # row counter for ATEs

for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])   
  ### Matching-adjusted indirect comparison (MAIC)
  load(paste0("Results/MAIC/means_", file.id, ".RData"))
  load(paste0("Results/MAIC/variances_", file.id, ".RData"))
  simulation.metrics[j,3] <- "MAIC"
  maic.metrics <- process.metrics(means, variances, Delta.AB) 
  simulation.metrics[j,4:16] <- unlist(maic.metrics)
  ate.table[k:(k+replicates-1),3] <- "MAIC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates
  ### Two-stage matching-adjusted indirect comparison (2SMAIC)
  load(paste0("Results/2SMAIC/means_", file.id, ".RData"))
  load(paste0("Results/2SMAIC/variances_", file.id, ".RData"))  
  simulation.metrics[j,3] <- "2SMAIC"
  maic.2s.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:16] <- unlist(maic.2s.metrics)
  ate.table[k:(k+replicates-1),3] <- "2SMAIC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates  
  ### Truncated matching-adjusted indirect comparison (T-MAIC)
  load(paste0("Results/T-MAIC/means_", file.id, ".RData"))
  load(paste0("Results/T-MAIC/variances_", file.id, ".RData"))    
  simulation.metrics[j,3] <- "T-MAIC"
  t.maic.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:16] <- unlist(t.maic.metrics)
  ate.table[k:(k+replicates-1),3] <- "T-MAIC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates    
  ### Truncated two-stage matching-adjusted indirect comparison (T-2SMAIC)
  load(paste0("Results/T-2SMAIC/means_", file.id, ".RData"))
  load(paste0("Results/T-2SMAIC/variances_", file.id, ".RData"))
  simulation.metrics[j,3] <- "T-2SMAIC"
  t.maic.2s.metrics <- process.metrics(means, variances, Delta.AB)
  simulation.metrics[j,4:16] <- unlist(t.maic.2s.metrics)
  ate.table[k:(k+replicates-1),3] <- "T-2SMAIC"
  ate.table[k:(k+replicates-1),4] <- means
  j <- j+1
  k <- k+replicates      
}

# Save simulation study performance metrics
write.csv(simulation.metrics, "Analysis/scenarios.csv", row.names = FALSE)

## Function generates ridgeline plot for a specific scenario
plot.results <- function(scenario) {
  i <- scenario
  scenario.ates <- subset(ate.table, meanX_AC==pc$meanX_AC[i]&N_AC==pc$N_AC[i])
  ridge.plot <- ggplot(scenario.ates, aes(x=ATE, y=Method, fill=Method)) +
    geom_density_ridges(alpha=0.65) +
    geom_vline(xintercept=0, linetype="dashed", color ="red") +
    scale_x_continuous(limits=c(-2, 2)) +
    scale_y_discrete(limits=c("T-2SMAIC", "T-MAIC", "2SMAIC", "MAIC")) +
    theme_classic() + 
    theme(legend.position = "none", 
          axis.text.y = element_text(color="grey20", size=9, 
                                     face ="plain"),
          axis.text.x = element_text(color="grey20", size=7, 
                                     face="plain"),
          plot.title = element_text(color="grey20", size=10, face ="plain"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(color="grey20", size=9, face="plain")) +
    scale_fill_brewer(palette="Dark2") + xlab("Point estimates")
  # a ridgeline plot with the spread of point estimates is returned
  return(ridge.plot)
}  

## Function generates results table for a specific scenario
table.results <- function(scenario) {
  i <- scenario
  scenario.metrics <- subset(simulation.metrics,
                             meanX_AC==pc$meanX_AC[i]&N_AC==pc$N_AC[i])
  display.table <- cbind(Method=scenario.metrics$Method,
                         Bias=paste0(format(round(scenario.metrics$Bias,digits=3),nsmall=3)," (",
                                     format(round(scenario.metrics$Bias.MCSE,digits=3),nsmall=3),")"),
                         Cov=paste0(format(round(scenario.metrics$Cov,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$Cov.MCSE,digits=3),nsmall=3),")"),
                         ESE=paste0(format(round(scenario.metrics$ESE,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$ESE.MCSE,digits=3),nsmall=3),")"),
                         MSE=paste0(format(round(scenario.metrics$MSE,digits=3),nsmall=3)," (",
                                    format(round(scenario.metrics$MSE.MCSE,digits=3),nsmall=3),")"))
  table.grob <- tableGrob(display.table, theme=ttheme_minimal(base_size=7))
  # a table with the performance measures and corresponding MCSEs is returned
  return(table.grob)
}

# ridgeline plot for each scenario
ridge.plot.s1 <- plot.results(scenario=1) + ggtitle(expression(paste(italic(n),
                                                                     "=140, strong overlap")))  
ridge.plot.s2 <- plot.results(scenario=2) + ggtitle(expression(paste(italic(n),
                                                                     "=200, strong overlap")))
ridge.plot.s3 <- plot.results(scenario=3) + ggtitle(expression(paste(italic(n),
                                                                     "=140, moderate overlap")))
ridge.plot.s4 <- plot.results(scenario=4) + ggtitle(expression(paste(italic(n),
                                                                     "=200, moderate overlap")))
ridge.plot.s5 <- plot.results(scenario=5) + ggtitle(expression(paste(italic(n),
                                                                     "=140, poor overlap")))
ridge.plot.s6 <- plot.results(scenario=6) + ggtitle(expression(paste(italic(n),
                                                                     "=200, poor overlap"))) 

# table of results for each scenario
table.grob.s1 <- table.results(scenario=1)  
table.grob.s2 <- table.results(scenario=2)  
table.grob.s3 <- table.results(scenario=3)  
table.grob.s4 <- table.results(scenario=4)  
table.grob.s5 <- table.results(scenario=5)  
table.grob.s6 <- table.results(scenario=6)  

ridge.grid <- arrangeGrob(ridge.plot.s1, table.grob.s1, ridge.plot.s2, table.grob.s2,
                          ridge.plot.s3, table.grob.s3, ridge.plot.s4, table.grob.s4,
                          ridge.plot.s5, table.grob.s5, ridge.plot.s6, table.grob.s6, 
                          ncol=2, widths=c(0.8,1.2)) 
ggsave(file="Analysis/Figure1.pdf", plot=ridge.grid, width=170, height=225, 
       units="mm", dpi = 300)