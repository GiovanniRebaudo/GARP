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
y = get(load("Data.RData"))
P = ncol(y)
N = nrow(y)

# Data scatter plot (Figure 1 in the main manuscript)
Plot_1 = pre_plot(y)
# CairoPNG(filename = '../Image/Mice_Data.png', width = 500, height = 400)
Plot_1
# invisible(dev.off())

# Gaussian edge contour plot (Figure S.1 in the supplementary materials)
Plot_S1 = edge_countorplot(vertices = rbind(c(-2,-2), c(3,3)),
                           data.grid =expand.grid(X = seq(-3, 4, length.out=800), 
                                                  Y = seq(-3, 4, length.out=800)))
Plot_S1

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
  # save(output, file="output.RData")
} else {
  load("output.RData")
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
  
  # save(output_edge, file="output_edge.RData")
} else {
  load("output_edge.RData")
}
attach(output_edge)

data_plot           = data.frame(y[is_i_stable,])
Main_phases         = factor(clust_VI_stable)
levels(Main_phases) = paste0("V", 1:c_clust_data_stable)
colnames(data_plot) = c("MDS1","MDS2")

data_plot_edge           = data.frame(y[!is_i_stable,])
colnames(data_plot_edge) = c("MDS1","MDS2")
Plot = ggplot()+ theme_bw() +
  labs(color="Vertex")+ 
  geom_point(data=data_plot,aes(x=MDS1, y=MDS2,col=Main_phases))+
  geom_point(data=data_plot_edge,aes(x=MDS1, y=MDS2), col="black")+
  theme(legend.position = "right", text = element_text(size=20))

K_T_max = c_clust_data_stable*(c_clust_data_stable-1)/2


Segment_data           = matrix(nrow=K_T_max, ncol=c_clust_data_stable)
colnames(Segment_data) = c("x","y","xend","yend")

for(ind in 1:K_T_max){
  Segment_data[ind,] = c(mu_stable_map[Map_k_edge[ind,1],1], 
                         mu_stable_map[Map_k_edge[ind,1],2],
                         mu_stable_map[Map_k_edge[ind,2],1], 
                         mu_stable_map[Map_k_edge[ind,2],2])
}

Seg1 = data.frame(t(Segment_data[1,]))
Seg2 = data.frame(t(Segment_data[2,]))
Seg3 = data.frame(t(Segment_data[3,]))
Seg4 = data.frame(t(Segment_data[4,]))
Seg5 = data.frame(t(Segment_data[5,]))
Seg6 = data.frame(t(Segment_data[6,]))

Alpha = table(cl_memb_edge_out[seq_thin,])
Alpha = Alpha/max(Alpha)*20

Plot = ggplot() +
  geom_point(data=data_plot,aes(x=MDS1, y=MDS2,col=Main_phases))+
  geom_segment(data=Seg1, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[1],size=1)+
  geom_segment(data=Seg2, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[2],size=1)+
  geom_segment(data=Seg3, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[3],size=1)+
  geom_segment(data=Seg4, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[4],size=1)+
  geom_segment(data=Seg5, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[5],size=1)+
  geom_segment(data=Seg6, mapping =aes(x = x, y = y, xend = xend, yend = yend),col="black",alpha=Alpha[6],size=1)+
  xlab("Dim 1")+ylab("Dim 2")+labs(color="Vertex")+
  theme_bw()+theme(legend.position = "right", text = element_text(size=20))

CairoPNG(filename = '../Image/Inference_Mice_Orange.png', width = 500, height = 400)
Plot
invisible(dev.off())
