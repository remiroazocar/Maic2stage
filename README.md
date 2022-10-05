# Two-stage matching-adjusted indirect comparison: Code

### Antonio Remiro-Azócar

### *remiroantonio@gmail.com*

### *2022*

This repository contains the R code used for my paper [Two-stage matching-adjusted indirect comparison][1]. 

## Utilizing the Scripts

In order to use this repository, the user must first download a copy to their local machine. The user must set the working directory to the location where the download was made. To run the pipeline, the user should then open and run the scripts in the following order:

|          Script          | Explanation                                                  |
| :----------------------: | ------------------------------------------------------------ |
|       `gen_data.R`       | Specifies the settings of the simulation study and saves them in `"./sim_settings.RData"`. Generates the data for the simulation study according to the settings (saving the data to the `"./Data/"` subdirectory). |
| `covariate_adjustment.R` | Performs the matching-adjusted indirect comparison methods on the simulated data (saving the point estimates and variances to the `"./Results/"` subdirectory) |
|       `analysis.R`       | Processes the results of the simulation study and computes and graphs the relevant performance metrics (saved to the `"./Analysis/"` subdirectory) |

The `functions.R` script contains a user-defined function for weight estimation in MAIC and functions to evaluate the performance measures of interest. The file `./Analysis/scenarios.csv` records the parameter values or settings for each scenario and the key performance measures associated with each, as presented in Figure 1 of the paper. 

The `./Example` subdirectory features example `R` code implementing matching-adjusted indirect comparison (MAIC) and two-stage matching-adjusted indirect comparison (2SMAIC), as per the supplementary material of the article. The code can be easily adapted to perform the truncated versions of the methods (T-MAIC and T-2SMAIC). 

In the simulation study, the `doSNOW` package is used to parallelize the performance of the methods, distributing the tasks to different cores of the computer. 

The code was prepared in `RStudio` using `R` version `3.6.3` in a Windows architecture, with a 64-bit operating system. The following packages and versions were used:

* `boot 1.3-28` required for use of the non-parametric bootstrap in all methods 
* `doSNOW 1.0.19` used in combination with `foreach()` to start up local clusters that distribute parallel tasks to different cores
* `dplyr 1.0.7` for data manipulation in the simulated data generation
* `ggplot2 3.3.5` to plot the simulation study results (Figure 1 in the article)
* `ggridges 0.5.3` to plot the simulation study results (Figure 1 in the article)
* `gridExtra 2.3` to plot the simulation study results (Figure 1 in the article)
* `MASS 7.3-55` to simulate correlated covariates in the data-generating process for the simulation study, drawing from a multivariate normal distribution 
* `parallel 3.6.3` to detect the number of CPU cores

[1]: https://doi.org/10.1186/s12874-022-01692-9
