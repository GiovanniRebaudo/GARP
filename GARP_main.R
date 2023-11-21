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
library(reshape2)      # version 1.4.4
library(Cairo)         # version 1.6-0
library(scales)        # version 1.2.1
library(plyr)          # version 1.8.8
library(mvtnorm)       # version 1.1-3
library(LaplacesDemon) # version 16.1.6
library(bayesm)        # version 3.1-5
library(mniw)          # version 1.0.1
library(xtable)        # version 1.8-4
library(ggpubr)        # version 0.6.0


# Load functions
source("GARP_fcts.R")

# Load data
y = get(load("./Data-and-Results/Data.RData"))
P = ncol(y)
N = nrow(y)

# If you want to save the plot
Save_Plot = TRUE

# Data scatter plot (Figure 1 in the main manuscript)
Plot_1 = pre_plot(y)
if(Save_Plot){
  CairoPNG(filename = './Image/Mice_Data.png', width = 500, height = 400)}
Plot_1
if(Save_Plot){invisible(dev.off())}

# Gaussian edge contour plot 
# (Figure S.1 in the supplementary materials)
Plot_S1 = edge_countorplot(vertices  = rbind(c(-2,-2), c(3,3)),
                           data.grid = expand.grid(X=seq(-3,4,length.out=800), 
                                                   Y=seq(-3,4,length.out=800)))

if(Save_Plot){
  CairoPNG(filename = './Image/Gauss_Edge_Countor.png', width=500, height=400)}
Plot_S1
if(Save_Plot){invisible(dev.off())}

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
 # pt1 = proc.time() # compute time
  output = GARP_MCMC(data    = y,
                     mu0     = mu0, 
                     kappa0  = kappa0,
                     nu0     = nu0,
                     Lambda0 = Lambda0,
                     p_s     = p_s,
                     Niter   = Niter,
                     Plot    = FALSE,
                     acc_p   = FALSE)
  # pt2 = proc.time()
  # pt2-pt1
  # save(output, file="./Data-and-Results/output.RData")
} else {
  load("./Data-and-Results/output.RData")
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

if(Save_Plot){
  CairoPNG(filename = './Image/Prob_Coclus_obs_Mice_Data.png', 
           width = 500, height = 400)}
Plot_2b
if(Save_Plot){invisible(dev.off())}

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
  
  # save(output_edge, file="./Data-and-Results/output_edge.RData")
} else {
  load("./Data-and-Results/output_edge.RData")
}
attach(output_edge)

# Scatter-plot of the scRNA data with GARP point estimate.
# (Figure 2a in the main manuscript)
Plot_2a = Plot_result_GARP(y                   = y,
                           is_i_stable         = is_i_stable, 
                           clust_VI_stable     = clust_VI_stable,
                           mu_stable_map       = mu_stable_map,
                           Map_k_edge          = Map_k_edge,
                           cl_memb_edge_out    = cl_memb_edge_out
)
if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_Mice.png', 
                       width = 500, height = 400)}
# Change the color palette
Plot_2a
set_palette(Plot_2a, "jco")
if(Save_Plot){invisible(dev.off())}

# Posterior distribution of the number of vertex: K_v
# (Table 4 - GARP - in the main manuscript)
Table_4_GARP = Freq_Kv(cl_samp=cl_memb_stable_out[seq_thin,])
# xtable(t(Table_4_GARP),digits = 4)

# Find markers  ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("clusterExperiment")
# BiocManager::install("scran")
# BiocManager::install("scuttle")


library(clusterExperiment) # version 2.18.2
library(dplyr)             # version 1.1.1
# load data before dimensionality reduction
load("./Data-and-Results/Mice_original_data.rda")
set.seed(20)

# Original genes and cells assigned to vertices
core_stable     = core[,is_i_stable]
core_stable_mat = assay(core_stable)

levels(clust_VI_stable) = c("V3","V2","V1","V4")
mapvalues(clust_VI_stable, from = 1:4, to = c(3, 2, 1, 4))

N_DE= 6
cluster.markers = scran::findMarkers(x=core_stable_mat, pval.type = "any", 
                                      groups=clust_VI_stable, "binom")

markers      = rownames(cluster.markers[[1]][1:N_DE,])
markers_cell = cbind(t(core_stable_mat[markers,]), clust_VI_stable)

Plot_Mark = markers_cell %>% data.frame %>%
  group_by(clust_VI_stable) %>%
  summarise_at(dplyr::all_of(markers), list(mean))

BH_pvalue = cluster.markers[[1]][1:N_DE,]$FDR
Table_3   = t(rbind(Plot_Mark, c(0, BH_pvalue)))

colnames(Table_3)    = c(paste("Vertex", 1:4), "FDR")

# Average within vertex-cluster gene expressions and FDRs 
# in the selected top 6 biomarkers.
# (Table 3 in the main manuscript)
Table_3[-1,]
# xtable(Table_3[-1,], display=c(rep("f",5),"e"))

# Heatmap of the log mean expressions in top 6 DE genes in the main phases.
# (Figure 3 in the main manuscript)
Plot_3  = Plot_heat_DE_main(Plot_Mark = Plot_Mark)
Plot_3  
if(Save_Plot){
ggsave(filename='./Image/heatmap_subset.png', plot=Plot_3, device="png",
       width = 15, height = 10, units="cm")}

# Heat-map log genetic expressions in top 6 DE genes in all cells 
# ordered by main phases
# (Left panel of figure 4 in the main manuscript)
Plot_4L = Plot_heat_DE(markers_cell    = markers_cell,
                       clust_VI_stable = clust_VI_stable)
Plot_4L
if(Save_Plot){ggsave(filename='./Image/Heatmap_all_cells_DE.png', plot=Plot_4L, 
       device="png", width = 15, height = 17, units="cm")}

# Boxplot genetic expressions (after log(· + 1) transformation) 
# in the top 6 DE genes in all cells in the different main phases
# (Right panel of figure 4 in the main manuscript)
Plot_4R = Boxplot_DE(markers_cell = markers_cell) + ggpubr::color_palette("jco")
if(Save_Plot){ggsave(filename='./Image/Boxplot_DE.png', plot=Plot_4R, 
                     device="png", width = 15, height = 17, units="cm")}

# Run the MCMC of RPM with independent atoms (non graph-aligned) ---------------
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_ind = GARP_MCMC(data    = data,
                         mu0     = mu0, 
                         kappa0  = kappa0,
                         nu0     = nu0,
                         Lambda0 = Lambda0,
                         p_s     = 1,
                         Niter   = Niter,
                         Plot    = TRUE,
                         acc_p   = FALSE)
  # save(output_ind, file="./Data-and-Results/output_ind.RData")
} else {
  load("./Data-and-Results/output_ind.RData")
}
attach(output_ind)

# Point estimate vertex-clustering
clust_VI_stable       = salso(cl_memb_all_out[seq_thin,], loss=VI())
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
dissimlar_stable = psm(cl_memb_all_out[seq_thin,])

# RPM with independent atoms. 
# Posterior probabilities of co-clustering of obs assigned to vertices.
# (Figure 5b in the main manuscript)
Plot_5b = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                           N_S_map          = N)

if(Save_Plot){CairoPNG(filename = './Image/Prob_Coclus_obs_Mice_Data_Orange.png', 
                       width = 500, height = 400)}
Plot_5b
if(Save_Plot){invisible(dev.off())}

# RPM with independent atoms. 
# Scatter-plot of the data and clustering point estimate.
# (Figure 5a in the main manuscript)

Plot_5a = Plot_result_GARP(y                   = y,
                           is_i_stable         = rep(1, nrow(y)), 
                           clust_VI_stable     = clust_VI_stable
)

if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_Mice_Orange.png',
                       width = 500, height = 400)}
Plot_5a+ggpubr::color_palette("jco")
if(Save_Plot){invisible(dev.off())}

# Posterior distribution of the number of vertex: K_v
# (Table 4 - independent atoms RPM - in the main manuscript)
Table_4_b = Freq_Kv(cl_samp=cl_memb_stable_out[seq_thin,])
# xtable(t(Table_4_b),digits = 4)

