
# method function ---------------------------------------------------------
library(boot)
npml_estimate_mean_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=5000
  temp_UU=1
  epi=0.001
  temp_pp=1
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      ww[k,]=mapply(function(x,y) recursion(x,y,NN[k],Y[k]),pp,UU)
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    UU[which(UU>(1-epi))]=1-epi
    UU[which(UU<epi)]=epi
    if(any(is.nan(UU))==T){
      break
    }
    if((sum((temp_UU-UU)^2)<1e-5)&sum((temp_pp-pp)^2)<1e-5){
      break
    }
    temp_pp=pp
    temp_UU=UU
  }
  return(sum(UU*pp))
}

npml_estimate_logit_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=5000
  temp_UU=1
  epi=0.001
  temp_pp=1
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      ww[k,]=mapply(function(x,y) recursion(x,y,NN[k],Y[k]),pp,UU)
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    UU[which(UU>(1-epi))]=1-epi
    UU[which(UU<epi)]=epi
    if(any(is.nan(UU))==T){
      break
    }
    if((sum((temp_UU-UU)^2)<1e-5)&sum((temp_pp-pp)^2)<1e-5){
      break
    }
    temp_pp=pp
    temp_UU=UU
  }
  return(sum(logit(UU)*pp))
}


strat_estimate_function=function(N,d){
  NN=N[d]
  result=mean(1/NN)
  return(result)
}
pooled_estimate_function=function(N,d){
  NN=N[d]
  result=length(NN)/sum(NN)
  return(result)
}

recursion <- function (pp, UU, NN, Y) {
  ww = pp*(1-UU)^(NN-Y) * UU^Y
  return(ww)
}
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=5000
  temp_UU=1
  epi=0.001
  temp_pp=1
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      ww[k,]=mapply(function(x,y) recursion(x,y,NN[k],Y[k]),pp,UU)
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    UU[which(UU>(1-epi))]=1-epi
    UU[which(UU<epi)]=epi
    if(any(is.nan(UU))==T){
      break
    }
    if((sum((temp_UU-UU)^2)<1e-5)&sum((temp_pp-pp)^2)<1e-5){
      break
    }
    temp_pp=pp
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}
logit=function(x){return (log(x/(1-x)))}

compare_mean=function(x){
  if(x[4,2]<=x[4,3]&x[4,2]<=x[4,4]){
    if(x[4,2]==x[4,3]){
      return("SP")
    }
    if(x[4,2]==x[4,4]){
      return("SN")
    }
    return("Strat")
  }
  
  if(x[4,3]<=x[4,2]&x[4,3]<=x[4,4]){
    if(x[4,3]==x[4,2]){
      return("SP")
    }
    if(x[4,3]==x[4,4]){
      return("PN")
    }
    return("Pooled")
  }
  
  if(x[4,4]<=x[4,2]&x[4,4]<=x[4,3]){
    if(x[4,4]==x[4,2]){
      return("SN")
    }
    if(x[4,4]==x[4,3]){
      return("PN")
    }
    return("Npml")
  }
  
}


compare_logit=function(x){
  if(x[12,2]<=x[12,3]&x[12,2]<=x[12,4]){
    if(x[12,2]==x[12,3]){
      return("SP")
    }
    if(x[12,2]==x[12,4]){
      return("SN")
    }
    return("Strat")
  }
  
  if(x[12,3]<=x[12,2]&x[12,3]<=x[12,4]){
    if(x[12,3]==x[12,2]){
      return("SP")
    }
    if(x[12,3]==x[12,4]){
      return("PN")
    }
    return("Pooled")
  }
  
  if(x[12,4]<=x[12,2]&x[12,4]<=x[12,3]){
    if(x[12,4]==x[12,2]){
      return("SN")
    }
    if(x[12,4]==x[12,3]){
      return("PN")
    }
    return("Npml")
  }
  
  
}

relative_mean_function=function(x){
  places=2
  temp=100*x[4,4]/x[4,3]
  return(round(temp,places))
}
relative_logit_function=function(x){
  places=2
  temp=100*x[12,4]/x[12,3]
  return(round(temp,places))
}
simulation_sd_bias_mean=function(R,x){
  places=4
  return(round(sqrt(x[3,4]/R),places))
}
simulation_sd_bias_logit=function(R,x){
  places=4
  return(round(sqrt(x[11,4]/R),places))
}
# three methods comparasion -----------------------------------------------

#theta_bar_vector(mean),M_vector(precision),K_vector(sample size) R=200(replicates)
theta_bar_vector=c(0.1)
M_vector=c(2)
K_vector=c(20)
delta_low=0.005
delta_high=0.995
for(i1 in 1:length(theta_bar_vector)){
  for(i2 in 1:length(M_vector)){
    for(i3 in 1:length(K_vector)){
      
      
      # here is the setting of the true value
      places=4
      R=100
      theta_bar=theta_bar_vector[i1]
      M=M_vector[i2]
      alpha=theta_bar*M
      beta=M*(1-theta_bar)
      K=K_vector[i3]
      cuts=1:K/(K+1)
      theta=qbeta(cuts,alpha,beta)
      theta[which(theta<delta_low)]=delta_low
      theta[which(theta>delta_high)]=delta_high
      true_theta_mean=mean(theta)
      true_theta_recip=mean(1/theta)
      true_theta_logit=mean(logit(theta))
      true_theta_var=var(theta)*(K-1)/K
      true_theta_recip_var=var(1/theta)*(K-1)/K
      true_theta_logit_var=var(logit(theta))*(K-1)/K
      
      
      
      #here we get the esimate of the three different methods
      stratified_estimate_results=NULL
      pooled_estimate_results=NULL
      npml_estimate_UU=NULL
      npml_estimate_pp=NULL
      boot_mean_var_result=NULL
      boot_logit_var_result=NULL
      NTL=sapply(theta,function(x) (rgeom(R,x)+1)) #generate the total N for all R replicates
      for(ind in 1:R){
        
        N=NTL[ind,]
        
        
        
        #apply the function to the data
        stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
        pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
        temp=npml_estimate_function(N,c(1:K))
        npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
        npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
        boot_mean_result=boot(N,npml_estimate_mean_function,1000)
        boot_logit_result=boot(N,npml_estimate_logit_function,1000)
        boot_mean_var_result=c(boot_mean_var_result,var(boot_mean_result$t))
        boot_logit_var_result=c(boot_logit_var_result,var(boot_logit_result$t))
      }
      
      
      
      
      #here we caculate the mean estimate of the three different methods
      stratified_estimate_bias=mean(stratified_estimate_results)-true_theta_mean
      stratified_estimate_var=var(stratified_estimate_results)
      stratified_estimate_mse=stratified_estimate_bias^2+stratified_estimate_var
      pooled_estimate_bias=mean(pooled_estimate_results)-true_theta_mean
      pooled_estimate_var=var(pooled_estimate_results)
      pooled_estimate_mse=pooled_estimate_bias^2+pooled_estimate_var
      npml_estimate_results=apply(npml_estimate_pp*npml_estimate_UU,2,sum)
      npml_estimate_bias=mean(npml_estimate_results)-true_theta_mean
      npml_estimate_var=var(npml_estimate_results)
      npml_estimate_mse=npml_estimate_bias^2+npml_estimate_var
      result_mean=list(true_theta_mean=round(true_theta_mean,places),true_theta_var=round(true_theta_var,places),strat_estimate_bias=round(stratified_estimate_bias,places),stratified_estimate_var=round(stratified_estimate_var,places),stratified_estimate_mse=round(stratified_estimate_mse,places),pooled_estimate_bias=round(pooled_estimate_bias,places),pooled_estimate_var=round(pooled_estimate_var,places),pooled_estimate_mse=round(pooled_estimate_mse,places),npml_estimate_bias=round(npml_estimate_bias,places),npml_estimate_var=round(npml_estimate_var,places),npml_estimate_mse=round(npml_estimate_mse,places))
      
      
      
      #here we caculate the recip estimate of the three defferent methods
      stratified_estimate_recip_bias=mean(1/stratified_estimate_results)-true_theta_recip
      stratified_estimate_recip_var=var(1/stratified_estimate_results)
      stratified_estimate_recip_mse=stratified_estimate_recip_bias^2+stratified_estimate_recip_var
      pooled_estimate_recip_bias=mean(1/pooled_estimate_results)-true_theta_recip
      pooled_estimate_recip_var=var(1/pooled_estimate_results)
      pooled_estimate_recip_mse=pooled_estimate_recip_bias^2+pooled_estimate_recip_var
      npml_estimate_recip_results=apply(npml_estimate_pp*(1/npml_estimate_UU),2,sum)
      npml_estimate_recip_bias=mean(npml_estimate_recip_results)-true_theta_recip
      npml_estimate_recip_var=var(npml_estimate_recip_results)
      npml_estimate_recip_mse=npml_estimate_recip_bias^2+npml_estimate_recip_var
      result_recip=list(true_theta_recip=round(true_theta_recip,places),true_theta_recip_var=round(true_theta_recip_var,places),strat_estimate_recip_bias=round(stratified_estimate_recip_bias,places),stratified_estimate_recip_var=round(stratified_estimate_recip_var,places),stratified_estimate_recip_mse=round(stratified_estimate_recip_mse,places),pooled_estimate_recip_bias=round(pooled_estimate_recip_bias,places),pooled_estimate_recip_var=round(pooled_estimate_recip_var,places),pooled_estimate_recip_mse=round(pooled_estimate_recip_mse,places),npml_estimate_recip_bias=round(npml_estimate_recip_bias,places),npml_estimate_recip_var=round(npml_estimate_recip_var,places),npml_estimate_recip_mse=round(npml_estimate_recip_mse,places))
      
      
      
      #here we caculate the logit estimate of the three different methods
      temp=(npml_estimate_UU==1)
      new_npml_estimate_UU=NULL
      new_npml_estimate_pp=NULL
      for(i in 1:ncol(temp)){
        if(sum(temp[,i]==F)==K){
          new_npml_estimate_UU=cbind(new_npml_estimate_UU,npml_estimate_UU[,i])
          new_npml_estimate_pp=cbind(new_npml_estimate_pp,npml_estimate_pp[,i])
        }
      }
      stratified_estimate_logit_bias=mean(logit(stratified_estimate_results))-true_theta_logit
      stratified_estimate_logit_var=var(logit(stratified_estimate_results))
      stratified_estimate_logit_mse=stratified_estimate_logit_bias^2+stratified_estimate_logit_var
      pooled_estimate_logit_bias=mean(logit(pooled_estimate_results))-true_theta_logit
      pooled_estimate_logit_var=var(logit(pooled_estimate_results))
      pooled_estimate_logit_mse=pooled_estimate_logit_bias^2+pooled_estimate_logit_var
      npml_estimate_logit_results=apply(new_npml_estimate_pp*(logit(new_npml_estimate_UU)),2,sum)
      npml_estimate_logit_bias=mean(npml_estimate_logit_results)-true_theta_logit
      npml_estimate_logit_var=var(npml_estimate_logit_results)
      npml_estimate_logit_mse=npml_estimate_logit_bias^2+npml_estimate_logit_var
      result_logit=list(true_theta_logit=round(true_theta_logit,places),true_theta_logit_var=round(true_theta_logit_var,places),strat_estimate_logit_bias=round(stratified_estimate_logit_bias,places),stratified_estimate_logit_var=round(stratified_estimate_logit_var,places),stratified_estimate_logit_mse=round(stratified_estimate_logit_mse,places),pooled_estimate_logit_bias=round(pooled_estimate_logit_bias,places),pooled_estimate_logit_var=round(pooled_estimate_logit_var,places),pooled_estimate_logit_mse=round(pooled_estimate_logit_mse,places),npml_estimate_logit_bias=round(npml_estimate_logit_bias,places),npml_estimate_logit_var=round(npml_estimate_logit_var,places),npml_estimate_logit_mse=round(npml_estimate_logit_mse,places))
      
      assign(paste("result_",theta_bar,"_",M,"_",K,"_list",sep=""),list(result_mean,result_recip,result_logit))
      assign(paste("result"),list(result_mean,result_recip,result_logit))
      method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
      strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
      pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
      npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
      true_value=c("E","Var")
      P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
      recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
      logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
      assign(paste("result_",theta_bar,"_",M,"_",K,sep=""),data.frame(method,strat,pooled,npml))
      assign(paste("true_",theta_bar,"_",M,"_",K,sep=""),data.frame(true_value,P,recip_P,logit_P))
      objectname=c(paste("result_",theta_bar,"_",M,"_",K,sep=""),paste("true_",theta_bar,"_",M,"_",K,sep=""))
      
      print(c(i1,i2,i3))
      
      
      
    }
    
  }
  
}




save.image("/home/student/hzhang1/R/project_Tom/result_Aug_16_2015_boot.Rdata")