#' Give kappa and e_ts and the ith row of phi_tab (+ corresponding alpha and g)

f_pseul_noconst <- function(i,e_ts,kap,phi_tab){sc_sfs <- as.vector(unname(phi_tab[i,-1]))
l_sfs <- length(sc_sfs)
e_tv <- kap*e_ts/(1+(kap-1)*e_ts)
#fact2 <- (1+kap)*e_ts
#fact2 <- fact2/(1+2*kap*e_ts)
out1 <- sum(as.integer(SFS$ts)*log((1-e_ts)*sc_sfs+e_ts*sc_sfs[l_sfs:1]))
if (e_ts>0){out1 <- out1 + Sneq$ts*log(2*kap*e_ts/(1+2*kap*e_ts))
out1 <- out1 - sum(SFS$ts)*log(1+2*kap*e_ts)}
out1 <- out1 + sum(as.integer(SFS$tv)*log((1-e_tv)*sc_sfs+e_tv*sc_sfs[l_sfs:1]))
if (e_ts>0){out1 <- out1 + Sneq$tv*log((1+kap)*e_ts/(1+2*kap*e_ts))
out1 <- out1 + sum(SFS$tv)*log((1+(kap-1)*e_ts)/(1+2*kap*e_ts))}
#Add the two additional combinatorial factors for the misid factor
if (e_ts>0){out1 <- out1 + lchoose(sum(SFS$ts)+Sneq$ts,Sneq$ts)
out1 <-out1 + lchoose(sum(SFS$tv)+Sneq$tv,Sneq$tv)}
return(out1)}



#' Read in phi table and define corresponding parameters 
phi_tab_alpha <- read.delim("powmil_E_EUR2.phi", 
                            header=FALSE)
phi_tab_alpha <- as.matrix(phi_tab_alpha)

n_alphas <- 100
alpha_steps <- 1 + (0:n_alphas)/n_alphas   
n_rhos <- 40
rho_steps <- (0:n_rhos)/n_rhos*10


misi <- seq(0,.15,length.out=16)
kappa_est <- (sum(SFS$tv)+Sneq$tv)/(2*(sum(SFS$ts)+Sneq$ts))

kappa <- c(1,2/3*kappa_est,kappa_est,3/2*kappa_est) # use kappa=1 for non-ts/tv version and 
# a naive kappa estimate
# given as  "# tv/# ts"   

beta_para <- expand.grid(g=rho_steps,alpha=alpha_steps,
                         e_ts=misi,kappa=kappa)
resb <- cbind(phirow=rep(1:nrow(phi_tab_alpha),length(misi)*length(kappa)),beta_para)  



res_v <- future_sapply(1:nrow(resb),function(j){f_pseul_noconst(resb$phirow[j],
                                                         e_ts = resb$e_ts[j],
                                                         kap = resb$kappa[j],
                                                         phi_tab = phi_tab_alpha)})

  res_list <- list(sfs_name,beta_all_k=NA,beta_k1=NA,beta_khat=NA,
                   BF_b_km=NA,BF_b_km_k1=NA,BF_b_km_khat=NA,
                   sumSFS_ts=sum(SFS$ts),sumSFS_tv=sum(SFS$tv),
                   sneqsum_ts=Sneq$ts,sneqsum_tv=Sneq$tv)
  
  whichb_kappaest <- sapply(resb$kappa,function(x){isTRUE(all.equal(x,kappa_est))})
  
  (res_list$beta_all_k <- resb[res_v==max(res_v),])
  
  (res_list$beta_k1 <- resb[resb$kappa==1,][which.max(res_v[resb$kappa==1]),])
  
  res_list$beta_khat <- resb[whichb_kappaest,][which.max(res_v[whichb_kappaest]),]
  
  #top5_b <- order(res_v,decreasing = TRUE)[1:5]
  #resb[top5_b,]
  
  #Optimize for all kappa values
  logml_beta <- max(res_v[resb$alpha!=2])
  logml_km <- max(res_v[resb$alpha==2])
  
  res_list$BF_b_km <- logml_beta-logml_km 
  
  #Optimize for kappa=1
  logml_beta <- max(res_v[resb$alpha!=2 & resb$kappa==1])
  logml_km <- max(res_v[resb$alpha==2 & resb$kappa==1])
  
  res_list$BF_b_km_k1 <- logml_beta-logml_km
  
  #Optimize for kappa=\hat{kappa}
  logml_beta <- max(res_v[resb$alpha!=2 & whichb_kappaest])
  logml_km <- max(res_v[resb$alpha==2 & whichb_kappaest])
  
  res_list$BF_b_km_khat <- logml_beta-logml_km
