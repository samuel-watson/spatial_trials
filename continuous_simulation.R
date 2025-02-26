# PACKAGES AND PRE-REQUISITES

require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
require(Matrix)
source("src/solarized.R")
# load required functions
Rcpp::sourceCpp("src/perm_test.cpp")
source("src/continuous_fn.R")

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- TRUE # save newly created data
example_plot <- TRUE # whether to generate an example plot of a generated dataset

#SIMULATION PARAMETERS
# set simulation parameters
# if data already exists with these parameters and use_data == TRUE then it will instead be loaded
n_seed <- 10
n_child <- 100
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
misspec <- TRUE # if true then a different function is used to simulate data and fit the model
del_e <- 0.25 # Delta_E
n_locs <- 16 # Number of sampled intervention locations
max_del <- max_de(n_locs) # Upper bound on Delta_E
beta <- -0.3 # max effect size
n_iter <- 1000 # number of simulation iterations
n_perm <- 1000 # number of permutation test/ bootstrap iterations
print_fit_progress <- TRUE # if true it will print the intermediate model fitting steps

source("src/run_continuous_simulation.R")