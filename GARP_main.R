# Codes accompanying "Graph-Aligned Random Partition Model (GARP)"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2) # version 3.4.2
theme_set(theme_bw(base_size = 14))
set.seed(123)

# Load functions
source("GARP_fcts.R")

# Load data
y = get(load("Data.RData"))
P = ncol(y)
N = nrow(y)

# Data scatter plot (Figure 1 in the manuscript)
Plot1 = pre_plot(y)
Plot1

set.seed(123)
# Run the MCMC ------------------------------------------------------------
run_MCMC = TRUE

# GARP hyperparameters
# Random partition parameters
p_s       = 0.5 # change to p_v
gamma_GN  = 0.5
alpha_Dir = 0.5
## NIG hyperparameters
mu0       = colMeans(y) 
kappa0    = 0.001
nu0       = 100 
Lambda0   = diag(rep(20,P))
Sigma0    = solve(Lambda0)

# MCMC quantities
Niter        = 1e5
burnin       = 4e3
thin         = 10
verbose_step = max(round(R/20),1)
