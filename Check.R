## Clear the workspace
rm(list = ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(Cairo)



Edge_Parameters_old = function(unrot_means=unr, qChi_2=qChi_2_99){
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
      rot_var   = qChi_2*rot_mat %*% diag(c(s2_11,s2_22)) %*% t(rot_mat)
    }
  }
  return(list(mean_edge=center_mean,var_edge=rot_var))
}


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

unr=rbind(c(-2,-2), c(3,3))
par=Edge_Parameters(unrot_means=unr, qChi_2=qchisq(0.99,2))
m <- par$mean_edge
sigma <- par$var_edge
data.grid <- expand.grid(X = seq(-3, 4, length.out=800), Y = seq(-3, 4, length.out=800))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))
Plot=ggplot(q.samp, aes(x=X, y=Y, z=prob)) + 
  geom_contour(colour = "red", show.legend=TRUE) +
  coord_fixed(xlim = c(-2, 3), ylim = unr[,2], ratio = 1) +
  xlab("Dim 1")+ylab("Dim 2")+
  geom_point(data=data.frame(cbind(unr,c(1,1))),aes(x=X1,y=X2,z=X3), size=3)+
  xlab("Dim 1")+ylab("Dim 2")+labs(color="Vertex")+
  theme_bw()+theme(legend.position = "right", text = element_text(size=20))
Plot

vertices = cbind(unr,c("mu_1","mu_2"),c(1,1))
colnames(vertices) = c("X1","X2","X3","X4")
vertices           = data.frame(vertices)
vertices$X1        = as.double(vertices$X1)+c(0.1,-0.42)
vertices$X2        = as.double(vertices$X2)+c(0.05,-0.1)
vertices$X4        = as.double(vertices$X4)

c(mean(vertices$X1),mean(vertices$X2))

Plot + #geom_segment(aes(x = -2, y = -2, xend = 2, yend = 1), colour = "blue",size=2) + geom_segment(aes(x = -0.5, y = 0, xend = 2.8, yend = -0.4), colour = "green",size=2)+
  annotate("text", x=-1.9, y=-1.9, label=latex2exp::TeX("$\\mu_1$",output="character"),
           hjust=0, size = 6, parse = TRUE)+    
  annotate("text", x=1.7, y=0.9, label=latex2exp::TeX("$\\mu_2$",output="character"),
           hjust=0, size = 6, parse = TRUE)+    
  annotate("text", x=-0.15, y=-0.4, label=latex2exp::TeX("$\\mu_{1,2}$",output="character"),
           hjust=0, size = 6, parse = TRUE)+
  geom_point(data=data.frame(cbind(0,-0.5,1)),aes(x=X1,y=X2,z=X3), size=3) 
