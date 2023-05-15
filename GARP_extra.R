# Extra Codes accompanying "Graph-Aligned Random Partition Model (GARP)"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi)    # version 0.14
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
library(dplyr)         # version 1.1.1

# Load functions
source("GARP_fcts.R")
set.seed(123)

# If you want to save the plot
Save_Plot = TRUE

## Sim 1
P = 2
ver_1_mean = c(7,0)
ver_2_mean = c(-3,6)
ver_3_mean = c(3,5)
ver_4_mean = c(-4,-5)
ver_5_mean = c(2,-4)

ver_1_var = diag(rep(1,2))/4
ver_2_var = diag(rep(1,2))/4
ver_3_var = diag(rep(1,2))/4
ver_4_var = diag(rep(1,2))/4
ver_5_var = diag(rep(1,2))/4

n_s = 200
n_e = 100

data_sim_ver_1 = mvrnorm(n_s, ver_1_mean, ver_1_var)
data_sim_ver_2 = mvrnorm(n_s, ver_2_mean, ver_2_var)
data_sim_ver_3 = mvrnorm(n_s, ver_3_mean, ver_3_var)
data_sim_ver_4 = mvrnorm(n_s, ver_4_mean, ver_4_var)
data_sim_ver_5 = mvrnorm(n_s, ver_5_mean, ver_5_var)

True_par          = Edge_Parameters(unrot_means = rbind(ver_5_mean,ver_4_mean))
data_sim_edg_1 = mvrnorm(n_e, True_par$mean_edge, True_par$var_edge)
True_par          = Edge_Parameters(unrot_means = rbind(ver_3_mean,ver_2_mean))
data_sim_edg_2 = mvrnorm(n_e, True_par$mean_edge, True_par$var_edge)
True_par          = Edge_Parameters(unrot_means = rbind(ver_1_mean,ver_5_mean))
data_sim_edg_3 = mvrnorm(n_e, True_par$mean_edge, True_par$var_edge)
True_par          = Edge_Parameters(unrot_means = rbind(ver_1_mean,ver_3_mean))
data_sim_edg_4 = mvrnorm(n_e, True_par$mean_edge, True_par$var_edge)
True_par          = Edge_Parameters(unrot_means = rbind(ver_5_mean,ver_3_mean))
data_sim_edg_5 = mvrnorm(n_e, True_par$mean_edge, True_par$var_edge)

data_sim = rbind(data_sim_ver_1, data_sim_ver_2, data_sim_ver_3, data_sim_ver_4,
                 data_sim_ver_5, data_sim_edg_1, data_sim_edg_2, data_sim_edg_3,
                 data_sim_edg_4, data_sim_edg_5)

data_sim = data_sim[,c(2,1)]

xlim_sim_data = c(min(data_sim[,1]), max(data_sim[,1]))
ylim_sim_data = c(min(data_sim[,2]), max(data_sim[,2]))
xylim_sim     = c(min(xlim_sim_data, ylim_sim_data), 
                  max(xlim_sim_data, ylim_sim_data))

Cluster = factor(c(rep("1",n_s), rep("5",n_s), rep("3",n_s), rep("4",n_s), 
                   rep("2",n_s), rep("2,4",n_e), rep("3,5",n_e), rep("1,2",n_e),
                   rep("1,3",n_e), rep("2,3",n_e)))

data_plot = cbind.data.frame(data.frame(data_sim),Cluster)
colnames(data_plot)=c("X","Y","Cluster")

data_plot_vertex = data_plot %>% filter(Cluster %in% 1:5)%>% 
  droplevels()
data_plot_edge   = data_plot %>% filter(!(Cluster %in% 1:5))%>% 
  droplevels()
colnames(data_plot_vertex) = c("X","Y","Cluster")

data_plot_vertex$Cluster         = factor(data_plot_vertex$Cluster)
levels(data_plot_vertex$Cluster) = c(paste0("V",1:5))

Segment_data           = matrix(nrow=5,ncol=4)
colnames(Segment_data) = c("x","y","xend","yend")

Segment_data[1,] = c(ver_1_mean,ver_3_mean)
Seg1 = data.frame(t(Segment_data[1,]))
Segment_data[2,] = c(ver_2_mean,ver_3_mean)
Seg2 = data.frame(t(Segment_data[2,]))
Segment_data[3,] = c(ver_5_mean,ver_4_mean)
Seg3 = data.frame(t(Segment_data[3,]))
Segment_data[4,] = c(ver_5_mean,ver_1_mean)
Seg4 = data.frame(t(Segment_data[4,]))
Segment_data[5,] = c(ver_3_mean,ver_5_mean)
Seg5 = data.frame(t(Segment_data[5,]))

Plot_S2 = ggplot() +geom_point(data=data_plot_vertex, aes(x=X,y=Y))+
  geom_point(data=data_plot_edge,aes(x=X,y=Y))+
  xlab("Dim 1")+ylab("Dim 2")+labs(color="Vertex")+
  geom_segment(data=Seg1, mapping =aes(x = y, y = x, xend = yend, yend = xend),
               col="red",alpha=3, size=1)+
  geom_segment(data=Seg2, mapping =aes(x = y, y = x, xend = yend, yend = xend),
               col="red",alpha=3, size=1)+
  geom_segment(data=Seg3, mapping =aes(x = y, y = x, xend = yend, yend = xend),
               col="red",alpha=3, size=1)+
  geom_segment(data=Seg4, mapping =aes(x = y, y = x, xend = yend, yend = xend),
               col="red",alpha=3, size=1)+
  geom_segment(data=Seg5, mapping =aes(x = y, y = x, xend = yend, yend = xend),
               col="red",alpha=3, size=1)+
  theme_bw()+theme(legend.position = "right", text = element_text(size=20))

# If you want to save the plot
Save_Plot = TRUE

if(Save_Plot){CairoPNG(filename = './Image/Data_sim.png', width = 500, 
                       height = 400)}
Plot_S2
if(Save_Plot){invisible(dev.off())}

y = data_sim
N = nrow(y)


# Sim 1 scatter plot 
# (Figure S.2 in the supplementary materials)
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
Lambda0   = diag(rep(150,P))

# MCMC quantities
Niter     = 1e4
run_MCMC  = TRUE #FALSE
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_1 = GARP_MCMC(data    = y,
                           mu0     = mu0, 
                           kappa0  = kappa0,
                           nu0     = nu0,
                           Lambda0 = Lambda0,
                           p_s     = p_s,
                           Niter   = Niter,
                           Plot    = TRUE,
                           acc_p   = FALSE)
  save(output_sim_1, file="./Data-and-Results/output_sim_1.RData")
} else {
  load("./Data-and-Results/output_sim_1.RData")
}
attach(output_sim_1)

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
# (Figure S3b in the supplementary)
Plot_S3b = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                            N_S_map          = N_S_map)

if(Save_Plot){
  CairoPNG(filename = './Image/Prob_Coclus_obs_Sim_Data.png', 
           width = 500, height = 400)}
Plot_S3b
if(Save_Plot){invisible(dev.off())}

# Edge assignments
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_1_edge = GARP_Edge(y                   = y,
                                is_i_stable         = is_i_stable, 
                                c_clust_data_stable = c_clust_data_stable,
                                kappa0              = kappa0,
                                nu0                 = nu0,
                                Lambda0             = Lambda0,
                                Niter               = Niter,
                                Plot                = TRUE)
  
  save(output_sim_1_edge, file="./Data-and-Results/output_sim_1_edge.RData")
} else {
  load("./Data-and-Results/output_sim_1_edge.RData")
}
attach(output_sim_1_edge)

Plot_S3a = Plot_result_GARP_sim(y                   = y,
                                is_i_stable         = is_i_stable, 
                                clust_VI_stable     = clust_VI_stable,
                                mu_stable_map       = mu_stable_map,
                                Map_k_edge          = Map_k_edge,
                                cl_memb_edge_out    = cl_memb_edge_out)


if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_Sim.png', 
                       width = 500, height = 400)}
Plot_S3a
if(Save_Plot){invisible(dev.off())}

# Run the MCMC of RPM with independent atoms (non graph-aligned) ---------------
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_1_ind = GARP_MCMC(data    = y,
                               mu0     = mu0, 
                               kappa0  = kappa0,
                               nu0     = nu0,
                               Lambda0 = Lambda0,
                               p_s     = 1,
                               Niter   = Niter,
                               Plot    = TRUE,
                               acc_p   = FALSE)
  save(output_sim_1_ind, file="./Data-and-Results/output_sim_1_ind.RData")
} else {
  load("./Data-and-Results/output_sim_1_ind.RData")
}
attach(output_sim_1_ind)

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
# (Figure extra in the supplementary)
Plot_extra = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                              N_S_map          = N)

if(Save_Plot){CairoPNG(filename = './Image/Prob_Coclus_obs_sim_1_ver.png', 
                       width = 500, height = 400)}
Plot_extra
if(Save_Plot){invisible(dev.off())}

# RPM with independent atoms. 
# Scatter-plot of the data and clustering point estimate.
# (Figure extra in the supplementary)

Plot_extra2 = Plot_result_GARP(y                   = y,
                               is_i_stable         = rep(1, nrow(y)), 
                               clust_VI_stable     = clust_VI_stable)

if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_sim_1_ver.png',
                       width = 500, height = 400)}
Plot_extra2
if(Save_Plot){invisible(dev.off())}

## Sim 2
rm(output_sim_1,output_sim_1_edge,output_sim_1_ind)
# Run the MCMC sim 2------------------------------------------------------------
set.seed(123)
ver_1_mean = c(0,7)
ver_2_mean = c(6,-3)
ver_3_mean = c(5,3)
ver_4_mean = c(-5,-4)
ver_5_mean = c(-4,2)

ver_1_var = diag(rep(1,2))/2
ver_2_var = diag(rep(1,2))/2
ver_3_var = diag(rep(1,2))/2
ver_4_var = diag(rep(1,2))/2
ver_5_var = diag(rep(1,2))/2

n_s = 200
n_e = 100

data_sim_ver_1 = mvrnorm(n=n_s, ver_1_mean, ver_1_var)
data_sim_ver_2 = mvrnorm(n=n_s, ver_2_mean, ver_2_var)
data_sim_ver_3 = mvrnorm(n=n_s, ver_3_mean, ver_3_var)
data_sim_ver_4 = mvrnorm(n=n_s, ver_4_mean, ver_4_var)
data_sim_ver_5 = mvrnorm(n=n_s, ver_5_mean, ver_5_var)
unif_bound = 1
unif_bias  = 0.25

data_sim_edg_1 = Edge_miss(unrot_means = rbind(ver_5_mean,ver_4_mean), 
                             bias=unif_bias, n=n_e,bound=unif_bound)
data_sim_edg_2 = Edge_miss(unrot_means =rbind(ver_3_mean,ver_2_mean), 
                                 bias=unif_bias, n=n_e,bound=unif_bound)
data_sim_edg_3 = Edge_miss(unrot_means =rbind(ver_1_mean,ver_5_mean), 
                                 bias=unif_bias, n=n_e,bound=unif_bound)
data_sim_edg_4 = Edge_miss(unrot_means=rbind(ver_1_mean,ver_3_mean), 
                                 bias=unif_bias, n=n_e,bound=unif_bound)
data_sim_edg_5 = Edge_miss(unrot_means = rbind(ver_5_mean,ver_3_mean), 
                                 bias=unif_bias, n=n_e,bound=unif_bound)

data_sim = rbind(data_sim_ver_1, data_sim_ver_2, data_sim_ver_3, data_sim_ver_4,
                 data_sim_ver_5, data_sim_edg_1, data_sim_edg_2, data_sim_edg_3,
                 data_sim_edg_4, data_sim_edg_5)

N        = nrow(data_sim)

xlim_sim_data = c(min(data_sim[,1]),max(data_sim[,1]))
ylim_sim_data = c(min(data_sim[,2]),max(data_sim[,2]))
xylim_sim     = c(min(xlim_sim_data,ylim_sim_data), 
                  max(xlim_sim_data,ylim_sim_data))

Plot_extra3 = ggplot() +geom_point(data=data.frame(data_sim), aes(x=X1, y=X2))+
  xlab("Dim 1")+ylab("Dim 2")
Plot_extra3

if(run_MCMC){
  mu0       = colMeans(data_sim) 
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_2 = GARP_MCMC(data    = data_sim,
                           mu0     = mu0, 
                           kappa0  = kappa0,
                           nu0     = nu0,
                           Lambda0 = Lambda0,
                           p_s     = p_s,
                           Niter   = Niter,
                           Plot    = TRUE,
                           acc_p   = FALSE)
  save(output_sim_2, file="./Data-and-Results/output_sim_2.RData")
} else {
  load("./Data-and-Results/output_sim_2.RData")
}
attach(output_sim_2)
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
# (Figure S.4 b in the supplementary)
Plot_S4b = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                            N_S_map          = N_S_map)

if(Save_Plot){
  CairoPNG(filename = './Image/Prob_Coclus_obs_Sim_Data_Rev.png', 
           width = 500, height = 400)}
Plot_S4b
if(Save_Plot){invisible(dev.off())}

# Edge assignments
if(run_MCMC){
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_2_edge = GARP_Edge(y                   = y,
                                is_i_stable         = is_i_stable, 
                                c_clust_data_stable = c_clust_data_stable,
                                kappa0              = kappa0,
                                nu0                 = nu0,
                                Lambda0             = Lambda0,
                                Niter               = Niter,
                                Plot                = TRUE)
  
  save(output_sim_2_edge, file="./Data-and-Results/output_sim_2_edge.RData")
} else {
  load("./Data-and-Results/output_sim_2_edge.RData")
}
attach(output_sim_2_edge)

Plot_S4a = Plot_result_GARP_sim(y                   = y,
                                is_i_stable         = is_i_stable, 
                                clust_VI_stable     = clust_VI_stable,
                                mu_stable_map       = mu_stable_map,
                                Map_k_edge          = Map_k_edge,
                                cl_memb_edge_out    = cl_memb_edge_out)


if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_Sim_Rev.png', 
                       width = 500, height = 400)}
Plot_S4a
if(Save_Plot){invisible(dev.off())}


## Sim 3
rm(output_sim_2, output_sim_2_edge)
# Run the MCMC sim 3------------------------------------------------------------
set.seed(123)
ver_1_mean = c(0,7)
ver_2_mean = c(6,-3)
ver_3_mean = c(5,3)
ver_4_mean = c(-5,-4)
ver_5_mean = c(-4,2)

ver_1_var = diag(rep(1,2))/4
ver_2_var = diag(rep(1,2))/4
ver_3_var = diag(rep(1,2))/4
ver_4_var = diag(rep(1,2))/4
ver_5_var = diag(rep(1,2))/4

n_s = 200

data_sim_ver_1 = mvrnorm(n=n_s, ver_1_mean, ver_1_var)
data_sim_ver_2 = mvrnorm(n=n_s, ver_2_mean, ver_2_var)
data_sim_ver_3 = mvrnorm(n=n_s, ver_3_mean, ver_3_var)
data_sim_ver_4 = mvrnorm(n_s, ver_4_mean, ver_4_var)
data_sim_ver_5 = mvrnorm(n=n_s, ver_5_mean, ver_5_var)

data_sim = rbind(data_sim_ver_1, data_sim_ver_2, data_sim_ver_3, 
                 data_sim_ver_4, data_sim_ver_5)

N        = nrow(data_sim)

xlim_sim_data = c(min(data_sim[,1]), max(data_sim[,1]))
ylim_sim_data = c(min(data_sim[,2]), max(data_sim[,2]))
xylim_sim     = c(min(xlim_sim_data, ylim_sim_data),
                  max(xlim_sim_data, ylim_sim_data))

data_plot           = data.frame(data_sim)

Plot_extra4 = ggplot() + geom_point(data=data_plot, aes(x=X1, y=X2)) + 
  labs(color="Vertex") + theme_bw() + xlab("Dim 1") + ylab("Dim 2") +
  theme(legend.position = "right", text = element_text(size=20))
Plot_extra4

if(run_MCMC){
  mu0       = colMeans(data_sim) 
  # Set the seed for reproducibility
  set.seed(123)
  output_sim_3 = GARP_MCMC(data    = data_sim,
                           mu0     = mu0, 
                           kappa0  = kappa0,
                           nu0     = nu0,
                           Lambda0 = Lambda0,
                           p_s     = p_s,
                           Niter   = Niter,
                           Plot    = TRUE,
                           acc_p   = FALSE)
  save(output_sim_3, file="./Data-and-Results/output_sim_3.RData")
} else {
  load("./Data-and-Results/output_sim_3.RData")
}
attach(output_sim_3)

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
# (Figure S.5 b in the supplementary)
Plot_S5b = Plot_heat_vertex(dissimlar_stable = dissimlar_stable,
                            N_S_map          = N_S_map)

if(Save_Plot){
  CairoPNG(filename = './Image/Prob_Coclus_obs_Sim_Data_rev2.png', 
           width = 500, height = 400)}
Plot_S5b
if(Save_Plot){invisible(dev.off())}

if(N_T_map>0){
  # Edge assignments
  if(run_MCMC){
    # Set the seed for reproducibility
    set.seed(123)
    output_sim_3_edge = GARP_Edge(y                   = data_sim,
                                  is_i_stable         = is_i_stable, 
                                  c_clust_data_stable = c_clust_data_stable,
                                  kappa0              = kappa0,
                                  nu0                 = nu0,
                                  Lambda0             = Lambda0,
                                  Niter               = Niter,
                                  Plot                = TRUE)
    
    save(output_sim_3_edge, file="./Data-and-Results/output_sim_3_edge.RData")
  } else {
    load("./Data-and-Results/output_sim_3_edge.RData")
  }
  attach(output_sim_3_edge)
}

Plot_S5a = Plot_result_GARP_sim(y                   = data_sim,
                                is_i_stable         = is_i_stable, 
                                clust_VI_stable     = clust_VI_stable,
                                mu_stable_map       = mu_stable_map,
                                Map_k_edge          = Map_k_edge,
                                cl_memb_edge_out    = cl_memb_edge_out)


if(Save_Plot){CairoPNG(filename = './Image/Inference_Scatter_Sim_rev2.png', 
                       width = 500, height = 400)}
Plot_S5a
if(Save_Plot){invisible(dev.off())}

