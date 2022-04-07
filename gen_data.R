# This file specifies the simulation setup and generates the simulation study data

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/Maic2stage")

rm(list=ls())

# Load packages
# package for data manipulation
if(!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
# package to sample/simulate the covariates from a multivariate normal
if(!require("MASS")) {install.packages("MASS"); library(MASS)}

set.seed(555) # set random seed for reproducibility

# Define simulation study parameters

N_sim <- 5000 # number of Monte Carlo replicates per scenario/method
allocation <- 1/2 # active treatment vs. placebo allocation ratio (1:1)

N_AC <- c(140,200) # number of subjects in the index (S=1) trial
N_BC <- 300 # number of subjects in the competitor (S=2) trial
K <- 3 # number of baseline covariates (all have main and interaction terms)
b_0 <- 5 # intercept of outcome-generating linear regression
b_trt <- -2 # conditional effect of active treatment vs. common comparator at baseline (x=0)
b_X <- 2 # main (prognostic) effect of covariates
b_int <- 1 # interaction effect of covariates
meanX_AC <- c(0.5,0.4,0.3) # mean of each normally-distributed covariate in S=1
meanX_BC <- 0.6 # mean of each normally-distributed covariate in S=2
sdX <- 0.4 # standard deviation of each covariate (same for both studies)
corX <- 0.2 # covariate pairwise correlation coefficient  
mean_eps <- 0 # mean of linear regression error terms
sd_eps <- 1 # std. deviation of linear regression error terms

# parameter combinations for each scenario
param.combinations <- expand.grid(N_AC=N_AC, meanX_AC=meanX_AC)
pc <- param.combinations

scenarios <- nrow(pc) # number of simulation scenarios
save(pc, N_sim, allocation, file="sim_settings.RData")

# generate simulated datasets for index and competitor trials
gen.data <- function(N, K, b_0, b_trt, b_X, b_int, meanX, sdX, corX, allocation,
                     mean_eps, sd_eps) {
  rho <- matrix(corX, nrow=K, ncol=K) # set correlation matrix
  diag(rho) <- rep(1, K)
  N_active <- round(N*allocation) # number of patients under active treatment
  N_control <- N - N_active # number of patients under control
  sd.vec <-rep(sdX, K) # vector of standard deviations
  cor2cov <- function(R, S) {
    # function to compute covariance matrix from correlation matrix R and vector
    # of standard deviations S. covariance matrix required as input for mvrnorm
    sweep(sweep(R, 1, S, "*"),2,S,"*")
  }
  cov.mat <- cor2cov(rho, sd.vec) # covariance matrix
  # simulate correlated continuous covariate values from multivariate normal
  # patients under active treatment
  X_active <- as.data.frame(MASS::mvrnorm(n=N_active,mu=rep(meanX,K), Sigma=cov.mat))
  # patients under control treatment
  X_control <- as.data.frame(MASS::mvrnorm(n=N_control,mu=rep(meanX,K), Sigma=cov.mat))  
  # all patients
  X <- rbind(X_active, X_control)
  colnames(X) <- c("X1","X2","X3") 
  # treatment assignment (1: active; 0: control)
  trt <- c(rep(1,N_active),rep(0,N_control)) 
  # generate outcomes using linear regression
  LP <- b_0 + b_X*X$X1 + b_X*X$X2 + b_X*X$X3 + b_trt*trt + 
    b_int*X$X1*trt + b_int*X$X2*trt + b_int*X$X3*trt # linear predictor
  # normally-distributed error terms
  eps <- rnorm(N, mean_eps, sd_eps)
  y <- LP + eps # continuous outcomes
  return(as.data.frame(cbind(X, trt, y)))
}

for (i in 1:scenarios) {
  print(i)
  # simulate IPD covariates and outcome for index trial (S=1)
  IPD.AC <- replicate(n=N_sim, expr=gen.data(pc$N_AC[i], K, b_0, b_trt, b_X, b_int, 
                                             pc$meanX_AC[i], sdX, corX, allocation,
                                             mean_eps, sd_eps),
                      simplify=FALSE)
  # simulate IPD covariates and outcome for competitor trial (S=2)
  IPD.BC <- replicate(n=N_sim, expr=gen.data(N_BC, K, b_0, b_trt, b_X, b_int, meanX_BC, 
                                             sdX, corX, allocation, mean_eps, sd_eps),
                      simplify=FALSE)
  # Summarize competitor study IPD as ALD
  ALD.BC <- lapply(1:N_sim, function(j) {
    as.data.frame(cbind(
      # aggregate the data for the competitor trial 
      summarise(IPD.BC[[j]], mean.X1=mean(X1), mean.X2=mean(X2), mean.X3=mean(X3),
                sd.X1=sd(X1), sd.X2=sd(X2), sd.X3=sd(X3),
                # estimated marginal mean difference for B vs. C with corresponding variance
                hat.Delta=summary(lm(y~trt, data=IPD.BC))$coef[2],
                var.hat.Delta=vcov(lm(y~trt, data=IPD.BC))[2,2])))
  } )
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  save(IPD.AC, file=paste0("Data/IPD_AC_", file.id, ".RData"))
  save(IPD.BC, file=paste0("Data/IPD_BC_", file.id, ".RData"))
  save(ALD.BC, file=paste0("Data/ALD_BC_", file.id, ".RData"))  
}                       
