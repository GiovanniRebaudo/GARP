# Codes accompanying "Graph-Aligned Random Partition Model (GARP)"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MASS)    # version 7.3-58.2
library(ggplot2) # version 3.4.2
theme_set(theme_bw(base_size = 14))
library(viridis) # version 0.6.2
set.seed(123)

# Load functions
source("GARP_fcts.R")

# Load data
y = get(load("Data.RData"))
P = ncol(y)
N = nrow(y)

# Data scatter plot (Figure 1 in the main manuscript)
Plot_1 = pre_plot(y)
Plot_1

# Gaussian edge contour plot (Figure S.1 in the supplemantary materials)
Plot_S1 =edge_countorplot(verices = rbind(c(-2,-2), c(3,3)))
Plot_S1

set.seed(123)
# Run the MCMC ------------------------------------------------------------
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

# MCMC quantities
Niter        = 4000
run_MCMC     = FALSE
if(run_MCMC){
  output = GARP_MCMC(data    = data,
                     mu0     = mu0, 
                     kappa0  = kappa0,
                     nu0     = nu0,
                     Lambda0 = Lambda0,
                     p_s     = p_s,
                     Niter   = Niter,
                     all     = TRUE,
                     Plot    = TRUE)
  # save(output,file="output.RData")
} else {
  load("output.RData")
}
attach(output)


# Plot


