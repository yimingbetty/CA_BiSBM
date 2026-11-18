# CA_BiSBM

R code for simulation studies of CA-BiSBM.

`method.R` implements the clustering methods. `generatedata.R`, `generatedata_nonormal.R`, and `generatedata_poisson.R` generate simulation data under Gaussian, non-Gaussian, and Poisson settings.

Install the required R packages: MASS, clue, mclust, ggplot2, and NAC. Keep all five scripts in the same directory, set it as the R working directory, and run `source("simulation.R")`.

Each method is run for 100 repetitions per setting. Results and figures are saved in `sim_outputs/`.