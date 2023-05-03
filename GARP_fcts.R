# GARP functions

## Function for preliminary scatter plot
pre_plot = function(data){
  data_plot = data.frame(Data)
  ggplot(data_plot, aes(x=MDS1,y=MDS2))+geom_point()+
  xlab("Dim 1")+ylab("Dim 2")+theme_bw()+
  theme(legend.position = "right", text = element_text(size=20))
}

## Function to compute the urn scheme of the Gnedin process
urn_GN_norm <- function(freq_minus, gamma_GN, K_S = K_S){
  #  K_S<-length(unique(freq_minus))
  tot_freq_minus = sum(freq_minus)
  unorm_out      = c((freq_minus+1)*(tot_freq_minus-K_S+gamma_GN),
                     K_S^2-K_S*gamma_GN)
  return(unorm_out/(tot_freq_minus^2+gamma_GN*tot_freq_minus))
}

## Function to compute the DM marginal likelihood 
urn_Dir_all_norm_div <- function(freq_minus,beta_DM,KM){
  unorm_out = freq_minus+beta_DM/KM
  return(unorm_out/sum(unorm_out))
}

## Function to compute the edge parameters given the location parameter 
## of the two vertices
Edge_Parameters = function(unrot_means=unr, qChi_2=qchisq(0.99,2)){
  diff_mu     = unrot_means[1,]-unrot_means[2,]
  if (all(diff_mu==0)){
    rot_var     = qChi_2*diag(c(1, 1))
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
        rot_var = qChi_2*diag(c(s2_22, s2_11))
      } else {
        rot_var = qChi_2*diag(c(s2_11, s2_22))
      }
    } else {
      dist_xy   = sqrt(sum(diff_mu^2))
      
      diff_mu1 = diff_mu[1]
      diff_mu2 = diff_mu[2]
      
      
      cos_theta = ifelse(diff_mu1>0, sqrt(diff_mu1^2/(diff_mu1^2+diff_mu2^2)), 
                         -sqrt(diff_mu1^2/(diff_mu1^2+diff_mu2^2)))
      sin_theta = ifelse(diff_mu2>0, sqrt(diff_mu2^2/(diff_mu1^2+diff_mu2^2)), 
                         -sqrt(diff_mu2^2/(diff_mu1^2+diff_mu2^2)))
      
      rot_mat       = matrix(c(cos_theta,sin_theta,-sin_theta,cos_theta),
                             nrow=2,ncol=2)
      s2_11     = dist_xy^2/4
      s2_22     = 2
      rot_var   = qChi_2*rot_mat %*% diag(c(s2_11,s2_22)) %*% t(rot_mat)
    }
  }
  return(list(mean_edge=center_mean,var_edge=rot_var))
}

# Function to plot Gaussian edge contour plot
edge_countorplot = function(vertices = rbind(c(-2,-2), c(3,3)),
                    data.grid =expand.grid(X = seq(-3, 4, length.out=800), 
                                           Y = seq(-3, 4, length.out=800))){
  par    <- Edge_Parameters(unrot_means=vertices, qChi_2=qchisq(0.99,2))
  m      <- par$mean_edge
  sigma  <- par$var_edge
  q.samp <- cbind(data.grid, prob=mvtnorm::dmvnorm(data.grid, mean=m, 
                                                   sigma=sigma))
  Plot=ggplot(q.samp, aes(x=X, y=Y, z=prob)) + 
    geom_contour(colour = "red", show.legend=F) +
    coord_fixed(xlim = c(-2, 3), ylim = vertices[,2], ratio = 1) +
    geom_point(data=data.frame(cbind(vertices,c(1,1))),
               aes(x=X1,y=X2,z=X3), size=3)+
    xlab("Dim 1")+ylab("Dim 2")+labs(color="Vertex")+
    theme_bw()+theme(legend.position = "right", text = element_text(size=20))
}





## Function to implement MCMC
GARP_MCMC = function(data      = data,
                     mu0       = mu0, 
                     kappa0    = kappa0,
                     nu0       = nu0,
                     Lambda0   = Lambda0,
                     p_s       = p_s,
                     Niter     = Niter,
                     Plot      = TRUE,
                     all       = TRUE){
  if(all==TRUE){
    
    verbose_step = max(round(Niter/20),1)
    
    # Quantities where we save MCMC output
    stable_out           = matrix(nrow=Niter, ncol=N)
    Mu_stable_out        = array(dim=c(N,P,Niter))
    Sigma_stable_out     = array(dim=c(P,P,N,Niter))
    Mu_edge_out          = array(dim=c(N,P,Niter))
    Sigma_edge_out       = array(dim=c(P,P,N,Niter))
    cl_memb_stable_out   = matrix(nrow=Niter, ncol=N)
    cl_memb_all_out      = matrix(nrow=Niter, ncol=N)
    n_clus_out           = matrix(nrow=Niter, ncol=2)
    
    # Compute fixed quantities for MCMC
    Sigma0      = solve(Lambda0)
    nu_new      = nu0+1
    kappa_new   = kappa0+1
    Scale_t_new = Sigma0*kappa_new/((nu0-P+1)*kappa0)
    df_t_new    = nu0-P+1
    
    # Compute outside MCMC p_knew: marginal likelihood of the single observations
    dmvt_new       = dmvt_new2 = double(N)
    for (i in 1:N){
      y_i          = y[i,]
      dmvt_new[i]  = LaplacesDemon::dmvt(y_i, mu=mu0, S=Scale_t_new, df=df_t_new)
    }
    
    # Cluster memberships initialization
    K_T            = 0
    K_S            = 10
    K_T_max        = K_S*(K_S-1)/2
    cl_memb_stable = kmeans(y, K_S,iter.max=25, nstart = 150)$cluster
    cl_memb_stable = plyr::mapvalues(cl_memb_stable, from=unique(cl_memb_stable), 
                                     to=1:K_S)
    stable         = rep(1,N)
    if(Plot){
      plot(y, col=cl_memb_stable, pch=16, xlab="dim 1", ylab="dim 2",
           main="Equal Scale", xlim=c(min(y),max(y)), ylim=c(min(y),max(y)))
    }
    
    # Unique kernel parameters Bivariate Gaussian
    Sigma_stable   = array(0,dim=c(P,P,K_S))
    Mu_stable      = matrix(0,nrow=K_S, ncol=P)
    
    # Initialize stable parameters
    for (ks in 1:K_S){
      i_ks             = which(cl_memb_stable==ks)
      n_ks             = length(i_ks)
      y_iks            = y[i_ks,]
      if(n_ks==1) {
        mean.y_ks      = y_iks
        S_ks           = 0
        kappa_kks      = kappa_new
        nu_ks          = nu_new
      } else if(n_ks==0){
        print("error K_S")
        stop()
      } else {
        mean.y_ks      = colMeans(y_iks)
        S_ks           = var(y_iks)*(n_ks-1)
        kappa_ks       = kappa0+n_ks
        nu_ks          = n_ks+nu0
      }
      
      Lambda_Sig_ks       = solve(Lambda0 + S_ks + 
                                    kappa0*n_ks/kappa_ks*(mean.y_ks-mu0)%*%
                                    t((mean.y_ks-mu0)))
      mu_muks             = (n_ks*mean.y_ks + kappa0*mu0)/kappa_ks
      
      Sigma_ks            = bayesm::rwishart(nu=nu_ks, V=Lambda_Sig_ks)$IW
      Sigma_stable[,,ks]  = Sigma_ks
      Mu_stable[ks,]      = mvrnorm(n=1, mu=mu_muks, Sigma=Sigma_ks/kappa_ks)
    }
    
    # Compute edge parameters
    K_T_max     = K_S*(K_S-1)/2
    Map_k_edge  = matrix(nrow=K_T_max,ncol=2)
    Sigma_edge  = array(0,dim=c(P,P,K_T_max))
    Mu_edge     = matrix(0,nrow=K_T_max, ncol=P)
    
    kt = 1
    for (k_s1 in 2:K_S){
      for (k_s2 in 1:(k_s1-1)){
        Map_k_edge[kt,]  = c(k_s1,k_s2)
        unr              = Mu_stable[c(k_s1,k_s2),]
        par_out          = Edge_Parameters(unrot_means=unr, qChi_2 = qChi_2_99)
        Mu_edge[kt,]     = par_out$mean_edge
        Sigma_edge[,,kt] = par_out$var_edge
        kt               = kt+1
      }
    }
    
    cl_memb_edge_mat = matrix(ncol=2,nrow=0)
    cl_memb_all_mat  = cbind(cl_memb_stable,cl_memb_stable)
    cl_memb_edge     = integer(0)
    
    cl_memb_all      = cl_memb_stable
    sorted_obs = 1:N
    
  # Gibbs
  for (r in 1:Niter){
    # Sample cluster membership and sample
    N_S          = sum(stable)
    N_T          = N - N_S
    K_S_plus     = K_S
    K_T_max      = K_S*(K_S-1)/2
    K_T_max_plus = K_T_max
    
    ind_empty          = F
    
    # sample cell belonging to transition phases before for efficiency
    obs_edge_old     = which(stable==0)
    obs_stable_old   = which(stable==1)
    
    sorted_obs = c(obs_edge_old,obs_stable_old)
    
    n_edge_all = integer(K_T_max_plus)
    
    # If there exist occupied transition phase
    if(N_T>0){
      stable_i_old          = 0
      
      # i belongs to a transition phase in the previous iteration 
      # (the first ones to be sampled)
      for (i in 1:N_T){
        i_edge              = obs_edge_old[i]
        y_i                 = y[i_edge,]
        
        obs_stable = which(stable==1)
        obs_edge   = which(stable==0)
        
        cl_memb_stable      = cl_memb_all[obs_stable]
        cl_memb_edge        = -cl_memb_all[obs_edge]
        
        n_stable_minus_i    = tabulate(cl_memb_stable)
        
        n_edge_all          = integer(K_T_max_plus)
        n_edge_minus_i      = tabulate(cl_memb_edge)
        K_pos               = length(n_edge_minus_i)
        n_edge_all[1:K_pos] = n_edge_minus_i
        
        n_edge_minus_i      = n_edge_all
        
        cl_memb_edge_i_old                 = -cl_memb_all[i_edge]
        n_edge_minus_i[cl_memb_edge_i_old] = n_edge_all[cl_memb_edge_i_old] -1
        
        # we use dmNorm. It can be faster to write our function to avoid 
        # repeating matrix operation on Sigma_mat every time
        b_minus_i_stable = p_s*urn_GN_norm(freq_minus=n_stable_minus_i, 
                                           gamma_GN, K_S) *
          c(mniw::dmNorm(y_i, mu=Mu_stable, Sigma=Sigma_stable), dmvt_new[i_stable])
        
        b_minus_i_edge = (1-p_s)*urn_Dir_all_norm_div(n_edge_minus_i,alpha_Dir,KM = K_S)*
          mniw::dmNorm(y_i, mu=Mu_edge, Sigma=Sigma_edge)
        
        b_minus_i_all     = c(b_minus_i_stable,b_minus_i_edge)
        
        cl_memb_i         = sample.int(n=K_T_max_plus+K_S_plus+1, size=1, 
                                       prob=b_minus_i_all, useHash = F)
        stable_i          = ifelse(cl_memb_i<K_S_plus+2,1,0)
        
        # If the cell i was previously in a transition phase moves to a stable one
        if (stable_i == 1){
          stable[i_edge] = stable_i
          
          cl_memb_all[i_edge]      = cl_memb_i
          cl_memb_all_mat[i_edge,] = c(cl_memb_i,cl_memb_i)
          
          # If we open a new stable cluster
          if (cl_memb_i == (K_S_plus+1)){
            # Sample new stable parameter
            Lambda_Sig_new         = solve(Lambda0 + kappa0/kappa_new * (y_i-mu0)%*%t(y_i-mu0))
            Mean_mu_new            = (kappa0*mu0+y_i)/kappa_new
            
            Sigma_stable_new       = bayesm::rwishart(nu=nu_new , V=Lambda_Sig_new)$IW
            Mu_stable_new          = mvtnorm::rmvnorm(n=1, mean=Mean_mu_new, sigma=Sigma_stable_new/kappa_new)
            
            Mu_stable              = rbind(Mu_stable, Mu_stable_new)
            Sigma_stable           = array(c(Sigma_stable, Sigma_stable_new), dim = c(2, 2, K_S_plus+1))
            
            # Compute new edge parameters associated 
            # (we are computing also useless one: it can be made slightly faster)
            # Save new parameter values
            kt                       = K_T_max_plus
            Mu_edge_temp             = Mu_edge
            Map_k_edge_temp          = Map_k_edge
            Sigma_edge_temp          = Sigma_edge
            
            Map_k_edge               = matrix(nrow=K_T_max_plus+K_S_plus, ncol=P)
            Mu_edge                  = matrix(nrow=K_T_max_plus+K_S_plus, ncol=P)
            Sigma_edge               = array(0,dim=c(P,P,K_T_max_plus+K_S_plus))
            
            Mu_edge[1:K_T_max_plus,]     = Mu_edge_temp
            Map_k_edge[1:K_T_max_plus,]  = Map_k_edge_temp
            Sigma_edge[,,1:K_T_max_plus] = Sigma_edge_temp
            
            for (k_s in 1:K_S_plus){
              kt                 = kt+1
              Map_k_edge[kt,]    = c(K_S_plus+1,k_s)
              unr                = rbind(Mu_stable[k_s,],Mu_stable_new)
              par_out            = Edge_Parameters(unrot_means=unr, 
                                                   qChi_2=qChi_2_99)
              Mu_edge[kt,]       = par_out$mean_edge
              Sigma_edge[,,kt]   = par_out$var_edge
            }
            
            K_S                    = K_S+1
            K_S_plus               = K_S_plus+1
            K_T_max                = K_S*(K_S-1)/2
            K_T_max_plus           = K_S_plus*(K_S_plus-1)/2
            
            n_edge_all            = integer(K_T_max_plus)
            n_edge_occ            = tabulate(cl_memb_edge)
            K_T_occ               = length(n_edge_occ)
            n_edge_all[1:K_T_occ] = n_edge_occ
            n_edge_minus_i        = n_edge_all
          } # If cell i remains in a transition phase
        } else if (stable_i == 0){
          cl_memb_i_edge      = cl_memb_i - (K_S_plus+1)
          
          # If cell i remains in the transition phase but changes edge                 
          if (cl_memb_i_edge!=cl_memb_edge_i_old){
            
            cl_memb_all[i_edge]      = -cl_memb_i_edge
            cl_memb_all_mat[i_edge,] = Map_k_edge[cl_memb_i_edge,]
          }
        }
      }
      K_T = length(unique(cl_memb_edge))
    } else if (N_T==0){
      cl_memb_edge = integer(0)
      K_T          = 0
    }
    
    # if i belongs to a stable phase in the previous iteration
    stable_i_old            = 1
    for (i in 1:N_S){
      break_ind              = F
      i_stable               = obs_stable_old[i]
      y_i                    = y[i_stable,]
      
      obs_stable = which(stable==1)
      obs_edge   = which(stable==0)
      
      cl_memb_stable         = cl_memb_all[obs_stable]
      cl_memb_edge           = -cl_memb_all[obs_edge]
      
      n_stable_minus_i       = tabulate(cl_memb_stable)
      n_edge_all             = integer(K_T_max_plus)
      if(length(cl_memb_edge)>0){
        n_edge_minus_i        = tabulate(cl_memb_edge)
        if(ind_empty){
          n_i_pos_stable       = which(n_stable_minus_i>0)
          n_i_pos_edge         = which((Map_k_edge[,1] %in% n_i_pos_stable)&
                                         (Map_k_edge[,2]%in%n_i_pos_stable))
          
          K_pos                = length(n_edge_minus_i)
          n_edge_all[1:K_pos]  = n_edge_minus_i
        } else {
          K_pos               = length(n_edge_minus_i)
          n_edge_all[1:K_pos] = n_edge_minus_i
        }
      }
      n_edge_minus_i         = n_edge_all
      
      cl_memb_stable_i_old                   = cl_memb_all[i_stable]
      n_stable_i                             = n_stable_minus_i[cl_memb_stable_i_old] 
      n_stable_minus_i[cl_memb_stable_i_old] = n_stable_i - 1
      # if we are emptying a stable phase
      if (n_stable_i==1){
        cl_memb_edge_mat = cl_memb_all_mat[obs_edge,]
        # if the stable phase has an edge associated cannot be emptied
        if(any(cl_memb_edge_mat==cl_memb_stable_i_old)){
          break_ind = T
        } else {
          # print(K_S-1)
          # We have emptied a stable phase
          K_S     = K_S-1
          K_T_max = K_S*(K_S-1)/2
          
          ind_empty = T
        }
      }
      
      if(!break_ind){
        # we use dmNorm. It can be faster to write our function to avoid 
        # repeating matrix operation on Sigma_mat every time 
        if(ind_empty){
          # It can be made more efficient
          b_minus_i_stable   = double(K_S_plus+1)
          n_i_pos_stable     = which(n_stable_minus_i>0)
          
          b_minus_i_edge     = double(K_T_max_plus)
          n_i_pos_edge       = which((Map_k_edge[,1] %in% n_i_pos_stable)&
                                       (Map_k_edge[,2]%in%n_i_pos_stable))
          
          b_minus_i_stable[c(n_i_pos_stable,K_S_plus+1)] = p_s*urn_GN_norm(
            freq_minus=n_stable_minus_i[n_i_pos_stable], gamma_GN, K_S) *
            c(mniw::dmNorm(y_i, mu=Mu_stable[n_i_pos_stable,], 
                           Sigma=Sigma_stable[,,n_i_pos_stable]),
              dmvt_new[i_stable])
          # It can be made more efficient
          b_minus_i_edge[n_i_pos_edge] = (1-p_s)*urn_Dir_all_norm_div(n_edge_minus_i[n_i_pos_edge],alpha_Dir,KM = K_S)*
            mniw::dmNorm(y_i, mu=Mu_edge[n_i_pos_edge,], Sigma=Sigma_edge[,,n_i_pos_edge])
          
        } else {
          b_minus_i_stable = p_s*urn_GN_norm(freq_minus=n_stable_minus_i, gamma_GN, K_S) *
            c(mniw::dmNorm(y_i, mu=Mu_stable, Sigma=Sigma_stable), dmvt_new[i_stable])
          # It can be made more efficient
          b_minus_i_edge = (1-p_s)*urn_Dir_all_norm_div(n_edge_minus_i,alpha_Dir,KM = K_S)*
            mniw::dmNorm(y_i, mu=Mu_edge, Sigma=Sigma_edge)
        }
        
        b_minus_i_all     = c(b_minus_i_stable,b_minus_i_edge)
        
        cl_memb_i         = sample.int(n=K_T_max_plus+K_S_plus+1,size=1, prob=b_minus_i_all, useHash = F)
        stable_i          = ifelse(cl_memb_i<(K_S_plus+2),1,0)
        
        # If the cell i was previously in a stable phase (and remain stable)
        if (stable_i == 1){
          # If cell i change (stable) cluster
          if(cl_memb_stable_i_old!=cl_memb_i){
            
            cl_memb_all[i_stable]  = cl_memb_i
            cl_memb_all_mat[i_stable,] = c(cl_memb_i,cl_memb_i)
            # If we open a new stable cluster
            if (cl_memb_i == (K_S_plus+1)){
              # print(K_S+1)
              # Sample new stable parameter
              Lambda_Sig_new         = solve(Lambda0 + kappa0/kappa_new * (y_i-mu0)%*%t(y_i-mu0))
              Mean_mu_new            = (kappa0*mu0+y_i)/kappa_new
              
              Sigma_stable_new       = bayesm::rwishart(nu=nu_new , V=Lambda_Sig_new)$IW
              Mu_stable_new          = mvtnorm::rmvnorm(n=1, mean=Mean_mu_new, sigma=Sigma_stable_new/kappa_new)
              
              Mu_stable              = rbind(Mu_stable, Mu_stable_new)
              Sigma_stable           = array(c(Sigma_stable, Sigma_stable_new), dim = c(2, 2, K_S_plus+1))
              
              # Compute new edge parameters associated 
              # (we are computing also useless one: it can be made slightly faster)
              kt                       = K_T_max_plus
              Mu_edge_temp             = Mu_edge
              Map_k_edge_temp          = Map_k_edge
              Sigma_edge_temp          = Sigma_edge
              
              Map_k_edge               = matrix(nrow=K_T_max_plus+K_S_plus, ncol=P)
              Mu_edge                  = matrix(nrow=K_T_max_plus+K_S_plus, ncol=P)
              Sigma_edge               = array(0,dim=c(P,P,K_T_max_plus+K_S_plus))
              
              Mu_edge[1:K_T_max_plus,]     = Mu_edge_temp
              Map_k_edge[1:K_T_max_plus,]  = Map_k_edge_temp
              Sigma_edge[,,1:K_T_max_plus] = Sigma_edge_temp
              
              for (k_s in 1:K_S_plus){
                kt                 = kt+1
                Map_k_edge[kt,]    = c(cl_memb_i,k_s)
                unr                = rbind(Mu_stable[k_s,],Mu_stable_new)
                par_out            = Edge_Parameters(unrot_means=unr, qChi_2 = qChi_2_99)
                
                Mu_edge[kt,]       = par_out$mean_edge
                Sigma_edge[,,kt]   = par_out$var_edge
              }
              # Save new parameter values
              K_S                    = K_S+1
              K_S_plus               = K_S_plus+1
              K_T_max                = K_S*(K_S-1)/2
              K_T_max_plus           = K_S_plus*(K_S_plus-1)/2
            }
          }
          # if cell i moves from a stable to a transition phase
        } else if (stable_i == 0){
          
          stable[i_stable]           = stable_i
          cl_memb_i_edge             = cl_memb_i - (K_S_plus+1)
          
          cl_memb_all[i_stable]      = -cl_memb_i_edge
          cl_memb_all_mat[i_stable,] = Map_k_edge[cl_memb_i_edge,]
          
        }
      }
    }
    
    # Reordering in each iteration
    obs_stable                 = which(stable==1)
    cl_memb_stable             = cl_memb_all[obs_stable]
    cl_memb_stable_uni         = unique(cl_memb_stable)
    cl_memb_stable             = plyr::mapvalues(cl_memb_stable,from=cl_memb_stable_uni,to=1:K_S)
    
    cl_memb_all_mat            = plyr::mapvalues(cl_memb_all_mat,from=cl_memb_stable_uni,to=1:K_S)
    
    # Check how to remap function that connect edge and stable phases
    obs_edge                      = which(stable==0)
    cl_memb_edge                  = cl_memb_all[obs_edge]
    cl_memb_edge_uni              = unique(cl_memb_edge)
    cl_memb_edge_mat              = cl_memb_all_mat[obs_edge,]
    cl_memb_edge_mat_uni          = unique(cl_memb_edge_mat)
    
    cl_memb_all[obs_stable]       = cl_memb_stable
    
    # Sample unique kernel parameters Bivariate Gaussian
    # Reorder parameters
    Sigma_stable   = Sigma_stable[,,cl_memb_stable_uni]
    Mu_stable      = Mu_stable[cl_memb_stable_uni,]
    
    
    # Compute associated edge parameters
    K_T_max     = K_S*(K_S-1)/2
    Map_k_edge  = matrix(nrow=K_T_max,ncol=2)
    Sigma_edge  = array(0,dim=c(P,P,K_T_max))
    Mu_edge     = matrix(0,nrow=K_T_max, ncol=P)
    kt = 1
    for (k_s1 in 2:K_S){
      for (k_s2 in 1:(k_s1-1)){
        Map_k_edge[kt,]  = c(k_s1,k_s2)
        unr              = Mu_stable[c(k_s1,k_s2),]
        par_out          = Edge_Parameters(unrot_means=unr, qChi_2 = qChi_2_99)
        Mu_edge[kt,]     = par_out$mean_edge
        Sigma_edge[,,kt] = par_out$var_edge
        kt               = kt+1
      }
    }
    
    # Finish to remap
    K_T                           = nrow(cl_memb_edge_mat_uni)
    if (K_T>0){
      cl_memb_edge_mat_uni        = unique(t(apply(cl_memb_edge_mat_uni, 1,function(x){sort(x,decreasing = T)})))
      K_T                         = nrow(cl_memb_edge_mat_uni)
      cl_memb_edge_uni_new        = double(K_T)
      for(k_t in 1:K_T){
        cl_memb_edge_uni_new[k_t] = which(apply(Map_k_edge, 1, function(x) {all(x==cl_memb_edge_mat_uni[k_t,])}))
      }
      cl_memb_edge                = plyr::mapvalues(cl_memb_edge,from=cl_memb_edge_uni,to=cl_memb_edge_uni_new)
    }
    
    
    cl_memb_all[obs_edge]   = -cl_memb_edge
    
    
    # MH step
    for (ks in 1:K_S){
      i_ks             = which(cl_memb_all==ks)
      n_ks             = length(i_ks)
      y_iks            = y[i_ks,]
      if(n_ks==1) {
        mean.y_ks      = y_iks
        S_ks           = 0
        kappa_ks       = kappa_new
        nu_ks          = nu_new
      } else {
        mean.y_ks         = colMeans(y_iks)
        S_ks              = var(y_iks)*(n_ks-1) 
        kappa_ks          = kappa0+n_ks
        nu_ks             = n_ks+nu0
      }
      Lambda_Sig_ks       = solve(Lambda0 + S_ks + kappa0*n_ks/kappa_ks*
                                    (mean.y_ks-mu0)%*%t((mean.y_ks-mu0)))
      mu_muks             = (n_ks*mean.y_ks + kappa0*mu0)/kappa_ks
      
      Sigma_ks_prop       = bayesm::rwishart(nu=nu_ks, V=Lambda_Sig_ks)$IW
      Mu_ks_prop          = mvrnorm(n=1, mu=mu_muks, Sigma=Sigma_ks/kappa_ks)
      
      Mu_edge_prop        = Mu_edge
      Sigma_edge_prop     = Sigma_edge
      # Compute acceptance probability
      edge_ks             = which(((Map_k_edge[,1]==ks)+(Map_k_edge[,2]==ks))==1)
      edge_ks             = edge_ks[edge_ks %in% cl_memb_edge]
      
      if (length(edge_ks)>0){
        log_rho_acc         = 0
        for (kedge in edge_ks){
          
          i_kedge        = which(cl_memb_all==-kedge)
          n_kedge        = length(i_kedge)
          y_i_kedge      = y[i_kedge,]
          
          if(n_kedge==1) {
            mean_y_kedge      = y_i_kedge
          } else {
            mean_y_kedge      = colMeans(y_i_kedge)
          }
          
          # Compute proposed new edge parameters
          k_s_k_s2                 = Map_k_edge[kedge,]
          k_s2                     = k_s_k_s2[k_s_k_s2!=ks]
          
          unr                      = Mu_stable[k_s_k_s2,]
          par_out                  = Edge_Parameters(unrot_means=unr, qChi_2=qChi_2_99)
          
          Mu_kedge_prop            = par_out$mean_edge
          Sigma_kedge_prop         = par_out$var_edge
          
          Mu_edge_prop[kedge,]     = Mu_kedge_prop
          Sigma_edge_prop[,,kedge] = Sigma_kedge_prop
          
          Mu_edge_old              = Mu_edge[kedge,]
          Sigma_edge_old           = Sigma_edge[,,kedge]
          
          log_rho_acc       = log_rho_acc + n_kedge/2*
            ( t((mean_y_kedge-Mu_edge_old)%*%solve(Sigma_edge_old)%*%
                  (mean_y_kedge-Mu_edge_old)) -
                t((mean_y_kedge-Mu_kedge_prop)%*%solve(Sigma_kedge_prop)%*%
                    (mean_y_kedge-Mu_kedge_prop))+
                log(det(Sigma_edge_old))-log(det(Sigma_kedge_prop)) )
        }
        
        move              = (log(runif(1)) < log_rho_acc)
        
        if(move){
          Sigma_stable[,,ks]  = Sigma_ks_prop
          Mu_stable[ks,]      = Mu_ks_prop
          
          Mu_edge             = Mu_edge_prop
          Sigma_edge          = Sigma_edge_prop
        }
      } else {
        Sigma_stable[,,ks]  = Sigma_ks_prop
        Mu_stable[ks,]      = Mu_ks_prop
      }
    }
    
    # Sample stable random proportion p_s
    N_S = length(cl_memb_stable)
    N_T = N-N_S
    
    # Save output
    # If you want n_clust
    n_clus_out[r,]                  = c(K_S,K_T)
    Mu_stable_out[1:K_S,,r]         = Mu_stable
    if(K_T>0){
      index_edge_out                = unique(cl_memb_edge)
      if(max(index_edge_out)<N){
        Mu_edge_out[index_edge_out,,r]     = Mu_edge[index_edge_out,]
        Sigma_edge_out[,,index_edge_out,r] = Sigma_edge[,,index_edge_out]
      }
    }
    Sigma_stable_out[,,1:K_S,r]     = Sigma_stable
    cl_memb_stable_out[r,stable==1] = cl_memb_stable
    cl_memb_all_out[r,]             = cl_memb_all
    stable_out[r,]                  = stable
    
    if(r%%verbose_step == 0){
      
      print(r) # print iteration
      
      if(Plot==T){
        col_clus = cl_memb_all
        K_ma = max(K_S,K_T)
        if(K_T==0){
          col_clus = viridis(K_ma)[col_clus]
        } else if(K_T==1){
          col_clus = ifelse(col_clus<0,1,viridis(K_ma)[abs(col_clus)])
        } else if (K_T>1){
          col_clus = ifelse(col_clus<0,1,
                            viridis(K_ma)[abs(col_clus)])
        }
        plot(y, col=col_clus, pch=16)
      }
    }
  }
  # Output
  return(list(n_clus_out         = n_clus_out,
              Mu_stable_out     = Mu_stable_out,
              Sigma_stable_out   = Sigma_stable_out,
              Mu_edge_out        = Mu_edge_out,
              Sigma_edge_out     = Sigma_edge_out,
              cl_memb_stable_out = cl_memb_stable_out,
              cl_memb_all_out    = cl_memb_all_out,
              stable_out         = stable_out
  ))
  }
}
