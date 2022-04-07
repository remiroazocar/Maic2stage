# This file performs the MAIC methods on the simulated data

rm(list=ls())

# setwd("C:/Users/gndtl/OneDrive - Bayer/Personal Data/New Papers/Maic2stage")
source('functions.R') # load MAIC function
load("sim_settings.RData") # load sim. study parameter combinations

# for detect cores
if(!require(parallel)) {install.packages("parallel"); library(parallel)}
# for parallel cluster
if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
# for non-parametric bootstrap in MAIC
if(!require(boot)) {install.packages("boot"); library(boot)}

set.seed(444) # set seed for reproducibility
scenarios <- nrow(pc) # number of simulation scenarios

# settings for the methods
resamples <- 2000 # Number of bootstrap resamples

# simulated patient-level (index) and aggregate-level (competitor) datasets for all scenarios
IPD.AC.all <- vector(mode="list", scenarios)
ALD.BC.all <- vector(mode="list", scenarios)
# load data
for (i in 1:scenarios) {
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i]) 
  load(paste0("Data/IPD_AC_", file.id, ".RData")) # load index patient-level data
  load(paste0("Data/ALD_BC_", file.id, ".RData")) # load competitor aggregate-level data
  IPD.AC.all[[i]] <- IPD.AC
  ALD.BC.all[[i]] <- ALD.BC
}

### One-stage matching-adjusted indirect comparison
maic.wrapper <- function(data.AC, data.BC, resamples, truncation=FALSE, trunc.cutoff=1) { 
  # Inputs: data.AC - index trial individual patient-level data; 
  # data.BC - comparator trial aggregate-level data;  
  # resamples - number of resamples for non-parametric bootstrap;
  # truncation - whether weights are truncated or not;
  # trunc.cutoff - percentile threshold for truncation
  maic.boot <- function(data, indices) {
    dat <- data[indices,]
    x.EM <- dat[,c("X1","X2","X3")] # IPD effect modifiers 
    theta <- data.BC[c("mean.X1", "mean.X2", "mean.X3")] # ALD effect modifier means
    # center the IPD effect modifiers on the ALD means
    x.EM$X1 <- x.EM$X1 - theta$mean.X1
    x.EM$X2 <- x.EM$X2 - theta$mean.X2
    x.EM$X3 <- x.EM$X3 - theta$mean.X3
    # trial assignment weights estimated using standard method of moments
    hat.w <- maic(X.EM=x.EM) # estimated weights
    # weight truncation
    if (truncation==TRUE) {
      # weights above specified percentile set to the value at the percentile
      w.cutoff <- quantile(hat.w, trunc.cutoff)
      hat.w <- ifelse(hat.w > w.cutoff, w.cutoff, hat.w)
    }
    # fit weighted linear regression
    outcome.fit <- lm(y~trt, weights=hat.w, data=dat)
    # fitted treatment coefficient is marginal treatment effect for A vs. C
    hat.Delta.AC <- coef(outcome.fit)["trt"] 
    return(hat.Delta.AC)
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.AC, statistic=maic.boot, R=resamples)
  # bootstrap mean of marginal A vs. C treatment effect estimate
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of A vs. C treatment effect estimate   
  hat.var.Delta.AC <- var(boot.object$t)
  # published marginal B vs. C treatment effect estimate
  hat.Delta.BC <- with(data.BC, hat.Delta)
  # published variance estimate for marginal B vs. C treatment effect
  hat.var.Delta.BC <- with(data.BC, var.hat.Delta)
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC  
  list(hat.Delta.AB, hat.var.Delta.AB)
}  

### Two-stage matching-adjusted indirect comparison
maic2s.wrapper <- function(data.AC, data.BC, resamples, truncation=FALSE, trunc.cutoff=1) { 
  # Inputs: data.AC - index trial individual patient-level data; 
  # data.BC - comparator trial aggregate-level data;  
  # resamples - number of resamples for non-parametric bootstrap;
  # truncation - whether weights are truncated or not;
  # trunc.cutoff - percentile threshold for truncation
  maic2s.boot <- function(data, indices) {
    dat <- data[indices,]
    ## TREATMENT ASSIGNMENT MODEL 
    trt.model <- glm(trt~X1+X2+X3, data = dat, family = binomial(link="logit")) 
    # predict propensity scores (probability of assignment to treatment A)
    dat$ps <- predict(trt.model, data = dat, type = "response")
    # inverse probability weights
    dat$w.trt <- ifelse(dat$trt==1, 1/dat$ps, 1/(1-dat$ps)) 
    ## TRIAL ASSIGNMENT MODEL
    x.EM <- dat[,c("X1","X2","X3")] # IPD effect modifiers 
    theta <- data.BC[c("mean.X1", "mean.X2", "mean.X3")] # ALD effect modifier means
    # center the IPD effect modifiers on the ALD means
    x.EM$X1 <- x.EM$X1 - theta$mean.X1
    x.EM$X2 <- x.EM$X2 - theta$mean.X2
    x.EM$X3 <- x.EM$X3 - theta$mean.X3
    # trial assignment weights estimated using standard method of moments
    dat$w.trial <- maic(X.EM=x.EM) 
    # combine weights
    hat.w <- dat$w.trt*dat$w.trial
    # weight truncation
    if (truncation==TRUE) {
      # weights above specified percentile set to the value at the percentile
      w.cutoff <- quantile(hat.w, trunc.cutoff)
      hat.w <- ifelse(hat.w > w.cutoff, w.cutoff, hat.w)
    }
    # fit weighted linear regression
    outcome.fit <- lm(y~trt, weights=hat.w, data=dat)
    # fitted treatment coefficient is marginal treatment effect for A vs. C
    hat.Delta.AC <- coef(outcome.fit)["trt"] 
    return(hat.Delta.AC)
  }
  # non-parametric bootstrap
  boot.object <- boot::boot(data=data.AC, statistic=maic2s.boot, R=resamples)
  # bootstrap mean of marginal A vs. C treatment effect estimate
  hat.Delta.AC <- mean(boot.object$t)
  # bootstrap variance of A vs. C treatment effect estimate   
  hat.var.Delta.AC <- var(boot.object$t)
  # published marginal B vs. C treatment effect estimate
  hat.Delta.BC <- with(data.BC, hat.Delta)
  # published variance estimate for marginal B vs. C treatment effect
  hat.var.Delta.BC <- with(data.BC, var.hat.Delta)
  hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC # A vs. B
  hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC  
  list(hat.Delta.AB, hat.var.Delta.AB)
}  

# set up cluster for parallel computing
num.cores <- detectCores()-1 # leave one available
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
# progress bar
pb <- txtProgressBar(max=N_sim, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

# combine lists in parallelisation
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Run MAIC methods for all replicates/scenarios in parallel
for(i in 1:scenarios) {
  IPD.AC <- IPD.AC.all[[i]]
  ALD.BC <- ALD.BC.all[[i]]
  file.id <- paste0("N_AC", pc$N_AC[i], "meanX_AC", pc$meanX_AC[i])  
  ### Standard one-stage matching-adjusted indirect comparison
  maic.results <- foreach(j=1:N_sim, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list()), .options.snow=opts,
                          .packages=c("boot")) %dopar% {
                          results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                  resamples)
                          return(results)
                        }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  save(means, file=paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file.id, ".RData"))
  ### Two-stage matching-adjusted indirect comparison
  maic2s.results <- foreach(j=1:N_sim, .combine='comb', .multicombine=TRUE,
                            .init=list(list(), list()), .options.snow=opts,
                            .packages=c("boot")) %dopar% {
                            results <- maic2s.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                      resamples)
                            return(results)
                          }
    close(pb)
    means <- unlist(maic2s.results[[1]])
    variances <- unlist(maic2s.results[[2]])
    save(means, file=paste0("Results/2SMAIC/means_", file.id, ".RData"))
    save(variances, file=paste0("Results/2SMAIC/variances_", file.id, ".RData"))  
    ### Truncated one-stage matching-adjusted indirect comparison
    t.maic.results <- foreach(j=1:N_sim, .combine='comb', .multicombine=TRUE,
                              .init=list(list(), list()), .options.snow=opts,
                              .packages=c("boot")) %dopar% {
                              results <- maic.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                      resamples, truncation=TRUE,
                                                      trunc.cutoff=0.95)
                              return(results)
                            }
    close(pb)
    means <- unlist(t.maic.results[[1]])
    variances <- unlist(t.maic.results[[2]])
    save(means, file=paste0("Results/T-MAIC/means_", file.id, ".RData"))
    save(variances, file=paste0("Results/T-MAIC/variances_", file.id, ".RData"))
    ### Truncated two-stage matching-adjusted indirect comparison
    t.maic2s.results <- foreach(j=1:N_sim, .combine='comb', .multicombine=TRUE,
                                .init=list(list(), list()), .options.snow=opts,
                                .packages=c("boot")) %dopar% {
                                results <- maic2s.wrapper(IPD.AC[[j]], ALD.BC[[j]],
                                                          resamples, truncation=TRUE,
                                                          trunc.cutoff=0.95)
                                return(results)
                              }
    close(pb)
    means <- unlist(t.maic2s.results[[1]])
    variances <- unlist(t.maic2s.results[[2]])
    save(means, file=paste0("Results/T-2SMAIC/means_", file.id, ".RData"))
    save(variances, file=paste0("Results/T-2SMAIC/variances_", file.id, ".RData"))
}

stopCluster(cluster)