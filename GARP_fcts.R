# GARP functions

## Function for preliminary scatter plot
pre_plot = function(data){
  data_plot = data.frame(Data)
  ggplot(data_plot, aes(x=MDS1,y=MDS2))+geom_point()+
  xlab("Dim 1")+ylab("Dim 2")+theme_bw()+
  theme(legend.position = "right", text = element_text(size=20))
}

## Function to compute the urn scheme of the Gnedin process
urn_GN_norm <- function(freq_minus, gamma_GN){
  #  K_S<-length(unique(freq_minus))
  tot_freq_minus = sum(freq_minus)
  unorm_out      = c((freq_minus+1)*(tot_freq_minus-K_S+gamma_GN),
                     K_S^2-K_S*gamma_GN)
  return(unorm_out/(tot_freq_minus^2+gamma_GN*tot_freq_minus))
}

## Function to compute the DM marginal likelihood 
DM_all <- function(freq_minus,beta_DM,KM){
  unorm_out = freq_minus+beta_DM/KM
  return(unorm_out/sum(unorm_out))
}

## Function to compute the edge parameters given the location parameter 
## of the two vertices
Edge_Parameters = function(unrot_means=unr, var=1, qChi_2=qChi_2_99){
  diff_mu     = unrot_means[1,]-unrot_means[2,]
  if (all(diff_mu==0)){
    rot_var     = var/qChi_2*diag(c(1, 1))
    center_mean = unr[1,]
    print("Error: distance is 0")
    stop()
  } else {
    center_mean  = colMeans(unrot_means)
    if (any(diff_mu==0)){
      dist_xy   = sqrt(sum(diff_mu^2))
      s2_11     = dist_xy^2/4
      s2_22     = 2
      if (diff_mu[1]==0){
        rot_var = var/qChi_2*diag(c(s2_22, s2_11))
      } else {
        rot_var = var/qChi_2*diag(c(s2_11, s2_22))
      }
    } else {
      dist_xy   = sqrt(sum(diff_mu^2))
      
      diff_mu1 = diff_mu[1]
      diff_mu2 = diff_mu[2]
      
      
      cos_theta = ifelse(diff_mu1>0, sqrt(diff_mu1^2/(diff_mu1^2+diff_mu2^2)), -sqrt(diff_mu1^2/(diff_mu1^2+diff_mu2^2)))
      sin_theta = ifelse(diff_mu2>0, sqrt(diff_mu2^2/(diff_mu1^2+diff_mu2^2)), -sqrt(diff_mu2^2/(diff_mu1^2+diff_mu2^2)))
      
      rot_mat       = matrix(c(cos_theta,sin_theta,-sin_theta,cos_theta),
                             nrow=2,ncol=2)
      # 
      # rot_var   = var/qChi_2*rot_mat %*% diag(c(dist_xy/2,1)^2) %*% t(rot_mat)
      # rot_var   = var*rot_mat %*% diag(c(dist_xy^2/(4*qChi_2)+1,1/qChi_2)) %*% t(rot_mat)
      s2_11     = dist_xy^2/4
      s2_22     = 2
      rot_var   = var/qChi_2*rot_mat %*% diag(c(s2_11,s2_22)) %*% t(rot_mat)
    }
  }
  return(list(mean_edge=center_mean,var_edge=rot_var))
}

## Function to implement MCMC
GARP_MCMC = function(data  = data,
                     a1    = a1,
                     a2    = a2,
                     Niter = Niter, 
                     burnin = burnin, 
                     thin   = thin){
  
}
