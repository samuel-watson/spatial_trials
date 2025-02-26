# PACKAGES AND PRE-REQUISITES
require(glmmrBase)
require(ggplot2)
require(sf)
require(patchwork)
require(ggforce)
source("src/solarized.R")

Rcpp::sourceCpp("src/perm_test.cpp")
source("src/discrete_fn.R")

# flags for data
use_data <- TRUE # use the simulated saved data
save_data <- TRUE # save newly created data
example_plot <- TRUE # whether to generate a plot of an example data set (See Figure S2 Supplementary Information)

#SIMULATION PARAMETERS
beta <- 0 # max absolute intervention effect
del_e <- max_dist <- 0.2 # delta -E
del_i <- 0.1 # ignore if using one-way model
n_locs <- 8 # number of intervention sites
radius <- 0.2 # radius of intervention areas
n_seed <- 10 # number of seed locations for observations
n_child <- 150 # number of children per seed
cov_pars <- c(0.25,0.5) # G.P. variance, length scale
max_del <- max_de(n_locs) # upper bound on del_e
adjust <- TRUE # whether to adjust the analysis
oneway <- TRUE # whether a two-way or one-way spillover model, one-way is external only
misspec <- FALSE # if true then the simulation model is a different function to the analysis model
sim_two <- FALSE # if true and oneway is true, then the simulation model will be two way, while model fitting is oneway
n_iter <- 1000 # number of simulation iterations
n_perm <- 1000 # number of permutation test/ bootstrap iterations
print_fit_progress <- TRUE # if true it will print the intermediate model fitting steps

source("src/run_discrete_simulation.R")



