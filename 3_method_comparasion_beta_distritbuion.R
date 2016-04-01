#Case 0.1-5-10 the beta distribution with mean 0.1   M=5(alpha+beta)  K=10(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.1
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=10
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}



stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
N=rep(0,K)
for(i in 1:K){
  N[i]=rgeom(1,theta[i])+1
}
stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
temp=npml_estimate_function(N,c(1:K))
npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.1_5_10_list=list(result_mean,result_recip,result_logit)
result=result_0.1_5_10_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.1_5_10=data.frame(method,strat,pooled,npml)
true_0.1_5_10=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.1_5_10,digits=4)
xtable(result_0.1_5_10,digits=4)
save(result_0.1_5_10,file="~/Documents/study/project/cristian/project1/result/result_0.1_5_10.Rdata")
save(true_0.1_5_10,file="~/Documents/study/project/cristian/project1/result/true_0.1_5_10.Rdata")







#Case 0.1-5-20 the beta distribution with mean 0.1   M=5(alpha+beta)  K=20(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.1
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=20
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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



stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.1_5_20_list=list(result_mean,result_recip,result_logit)
result=result_0.1_5_20_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.1_5_20=data.frame(method,strat,pooled,npml)
true_0.1_5_20=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.1_5_20,digits=4)
xtable(result_0.1_5_20,digits=4)
save(result_0.1_5_20,file="~/Documents/study/project/cristian/project1/result/result_0.1_5_20.Rdata")
save(true_0.1_5_20,file="~/Documents/study/project/cristian/project1/result/true_0.1_5_20.Rdata")






#Case 0.1-5-50 the beta distribution with mean 0.1   M=5(alpha+beta)  K=50(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.1
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=50
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
  return(list(UU,pp))
}
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.1_5_50_list=list(result_mean,result_recip,result_logit)
result=result_0.1_5_50_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.1_5_50=data.frame(method,strat,pooled,npml)
true_0.1_5_50=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.1_5_50,digits=4)
xtable(result_0.1_5_50,digits=4)
save(result_0.1_5_50,file="~/Documents/study/project/cristian/project1/result/result_0.1_5_50.Rdata")
save(true_0.1_5_50,file="~/Documents/study/project/cristian/project1/result/true_0.1_5_50.Rdata")








#Case 0.1-5-100 the beta distribution with mean 0.1   M=5(alpha+beta)  K=100(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.1
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=100
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.1_5_100_list=list(result_mean,result_recip,result_logit)
result=result_0.1_5_100_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.1_5_100=data.frame(method,strat,pooled,npml)
true_0.1_5_100=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.1_5_100,digits=4)
xtable(result_0.1_5_100,digits=4)
save(result_0.1_5_100,file="~/Documents/study/project/cristian/project1/result/result_0.1_5_100.Rdata")
save(true_0.1_5_100,file="~/Documents/study/project/cristian/project1/result/true_0.1_5_100.Rdata")









#Case 0.2-5-10 the beta distribution with mean 0.2   M=5(alpha+beta)  K=10(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.2
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=10
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.2_5_10_list=list(result_mean,result_recip,result_logit)
result=result_0.2_5_10_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.2_5_10=data.frame(method,strat,pooled,npml)
true_0.2_5_10=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.2_5_10,digits=4)
xtable(result_0.2_5_10,digits=4)
save(result_0.2_5_10,file="~/Documents/study/project/cristian/project1/result/result_0.2_5_10.Rdata")
save(true_0.2_5_10,file="~/Documents/study/project/cristian/project1/result/true_0.2_5_10.Rdata")










#Case 0.2-5-20 the beta distribution with mean 0.2   M=5(alpha+beta)  K=20(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.2
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=20
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.2_5_20_list=list(result_mean,result_recip,result_logit)
result=result_0.2_5_20_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.2_5_20=data.frame(method,strat,pooled,npml)
true_0.2_5_20=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.2_5_20,digits=4)
xtable(result_0.2_5_20,digits=4)
save(result_0.2_5_20,file="~/Documents/study/project/cristian/project1/result/result_0.2_5_20.Rdata")
save(true_0.2_5_20,file="~/Documents/study/project/cristian/project1/result/true_0.2_5_20.Rdata")



#Case 0.2-5-50 the beta distribution with mean 0.2   M=5(alpha+beta)  K=50(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.2
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=50
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.2_5_50_list=list(result_mean,result_recip,result_logit)
result=result_0.2_5_50_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.2_5_50=data.frame(method,strat,pooled,npml)
true_0.2_5_50=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.2_5_50,digits=4)
xtable(result_0.2_5_50,digits=4)
save(result_0.2_5_50,file="~/Documents/study/project/cristian/project1/result/result_0.2_5_50.Rdata")
save(true_0.2_5_50,file="~/Documents/study/project/cristian/project1/result/true_0.2_5_50.Rdata")











#Case 0.2-5-100 the beta distribution with mean 0.2   M=5(alpha+beta)  K=100(sample size) R=200(replicates)
places=4
R=200
theta_bar=0.2
M=5
alpha=theta_bar*M
beta=M*(1-theta_bar)
K=100
cuts=1:K/(K+1)
theta=qbeta(cuts,alpha,beta)
true_theta_mean=mean(theta)
true_theta_recip=mean(1/theta)
true_theta_logit=mean(logit(theta))
true_theta_var=var(theta)*(K-1)/K
true_theta_recip_var=var(1/theta)*(K-1)/K
true_theta_logit_var=var(logit(theta))*(K-1)/K



#method I the stratified estimate#
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
npml_estimate_function=function(N,d){
  NN=N[d]
  KK=length(NN)
  Y=rep(1,KK)
  UU = 1:KK/(KK+1)
  ww=matrix(1, nrow = KK, ncol=KK)
  pp=rep(1/KK, length = KK)
  Niter=10000
  temp_UU=1
  
  # The recursion
  
  for(i in 1:Niter)
  {   #Niter
    
    for(k in 1:KK){
      for(j in 1:KK) {
        ww[k,j] = pp[j]*(1-UU[j])^(NN[k]-Y[k]) * UU[j]^Y[k]
      }
      
      ww[k,] = ww[k,]/sum(ww[k,])
    }  #k to K
    
    pp=apply(ww,2,sum)/sum(ww)
    
    for(j in 1:KK) {
      UU[j] = sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    if(any(is.nan(UU))==T){
      break
    }
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  if(any(is.nan(UU))==T){
    return(list(NULL,NULL))
  }
  else{
    return(list(UU,pp))
  }
}


stratified_estimate_results=NULL
pooled_estimate_results=NULL
npml_estimate_UU=NULL
npml_estimate_pp=NULL
for(ind in 1:R){
  N=rep(0,K)
  for(i in 1:K){
    N[i]=rgeom(1,theta[i])+1
  }
  stratified_estimate_results=c(stratified_estimate_results,strat_estimate_function(N,c(1:K)))
  pooled_estimate_results=c(pooled_estimate_results,pooled_estimate_function(N,c(1:K)))
  temp=npml_estimate_function(N,c(1:K))
  npml_estimate_UU=cbind(npml_estimate_UU,temp[[1]])
  npml_estimate_pp=cbind(npml_estimate_pp,temp[[2]])
}
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

result_0.2_5_100_list=list(result_mean,result_recip,result_logit)
result=result_0.2_5_100_list
method=c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
strat=c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
pooled=c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
npml=c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
true_value=c("E","Var")
P=c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
recip_P=c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
logit_P=c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
result_0.2_5_100=data.frame(method,strat,pooled,npml)
true_0.2_5_100=data.frame(true_value,P,recip_P,logit_P)
xtable(true_0.2_5_100,digits=4)
xtable(result_0.2_5_100,digits=4)
save(result_0.2_5_100,file="~/Documents/study/project/cristian/project1/result/result_0.2_5_100.Rdata")
save(true_0.2_5_100,file="~/Documents/study/project/cristian/project1/result/true_0.2_5_100.Rdata")


