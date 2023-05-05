# Codes accompanying "Graph-Aligned Random Partition Model (GARP)"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.14
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MASS)          # version 7.3-58.2
library(ggplot2)       # version 3.4.2
theme_set(theme_bw(base_size = 14))
library(viridis)       # version 0.6.2
library(salso)         # version 0.3.29
library(reshape)       # version 0.8.9
library(Cairo)         # version 1.6-0
library(scales)        # version 1.2.1
library(mvtnorm)       # version 1.1-3
library(LaplacesDemon) # version 16.1.6

# Load functions
source("GARP_fcts.R")

# Load data
y = get(load("../Data-and-Results/Data.RData"))
P = ncol(y)
N = nrow(y)

# Data scatter plot (Figure 1 in the main manuscript)
Plot_1 = pre_plot(y)
# CairoPNG(filename = '../Image/Mice_Data.png', width = 500, height = 400)
Plot_1
# invisible(dev.off())

# Gaussian edge contour plot (Figure S.1 in the supplementary materials)
Plot_S1 = edge_countorplot(vertices  = rbind(c(-2,-2), c(3,3)),
                           data.grid = expand.grid(X=seq(-3,4,length.out=800), 
                                                   Y=seq(-3,4,length.out=800)))
Plot_S1

# Run the MCMC -----------------------------------------------------------------
# GARP hyperparameters
# Random partition parameters
p_s       = 0.5 # change to p_v
gamma_GN  = 0.5
alpha_Dir = 0.5
## NIG hyperparameters
mu0       = colMeans(y) 
kappa0    = 0.001
nu0       = 100 
Lambda0   = diag(rep(15,P))

# MCMC quantities
Niter     = 10000
run_MCMC  = FALSE
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
#  pt1 = proc.time() compute time
  output = GARP_MCMC(data    = data,
                     mu0     = mu0, 
                     kappa0  = kappa0,
                     nu0     = nu0,
                     Lambda0 = Lambda0,
                     p_s     = p_s,
                     Niter   = Niter,
                     all     = TRUE,
                     Plot    = TRUE,
                     acc_p   = TRUE)
  # pt2 = proc.time()
  # pt2-pt1
  # save(output, file="../Data-and-Results/output.RData")
} else {
  load("../Data-and-Results/output.RData")
}
attach(output)

# Thinning
thin         = 2
burn.in      = Niter/2
seq_thin     = seq(from=burn.in, to=Niter, by=thin)
Niter_ps     = length(seq_thin)

# Assign cells to vertex/edge phases
p_s_i        = colMeans(stable_out[seq_thin,])
is_i_stable  = (p_s_i>0.5)
(N_S_map     = sum(is_i_stable))
(N_T_map     = N-N_S_map)

# Point estimate vertex-clustering
clust_VI_stable       = salso(cl_memb_all_out[seq_thin,is_i_stable], loss=VI())
uni_clust_data_stable = unique(clust_VI_stable)

# Number of vertex clusters 
c_clust_data_stable       = length(uni_clust_data_stable)

# Frequencies of vertex clusters
freq_clust_VI_stable      = double(c_clust_data_stable)
for (i in 1:c_clust_data_stable){
  uni_clust_stable        = uni_clust_data_stable[i]
  freq_clust_VI_stable[i] = sum(clust_VI_stable==uni_clust_stable)
}
freq_clust_VI_stable

# Compute probability of co-clustering of cells assigned to vertices
dissimlar_stable = psm(cl_memb_all_out[seq_thin,is_i_stable])

# Posterior probabilities of co-clustering of obs assigned to vertices.
# (Figure 2b in the main manuscript)
Plot_2b = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                           N_S_map          = N_S_map)


# CairoPNG(filename = '../Image/Prob_Coclus_obs_Mice_Data_Orange.png', 
#          width = 500, height = 400)
Plot_2b
# invisible(dev.off())

# Edge assignments
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_edge = GARP_Edge(y                   = y,
                          is_i_stable         = is_i_stable, 
                          c_clust_data_stable = c_clust_data_stable,
                          kappa0              = kappa0,
                          nu0                 = nu0,
                          Lambda0             = Lambda0,
                          Niter               = Niter,
                          Plot                = TRUE)
  
  # save(output_edge, file="../Data-and-Results/output_edge.RData")
} else {
  load("output_edge.RData")
}
attach(output_edge)


Plot_2a = Plot_result_GARP(y                   = y,
                           is_i_stable         = is_i_stable, 
                           clust_VI_stable     = clust_VI_stable,
                           mu_stable_map       = mu_stable_map,
                           Map_k_edge          = Map_k_edge,
                           cl_memb_edge_out    = cl_memb_edge_out
)

# CairoPNG(filename = '../Image/Inference_Scatter_Mice.png', 
# width = 500, height = 400)
Plot_2a
# invisible(dev.off())

# Find markers  ----------------------------------------------------------------
# load data before dimensionality reduction
load("../Data-and-Results/Mice_original_data.rda")
