###Case I###
##we set the sample size as 50. there are 25 of them from the theta_1=0.1 and 25 of them from the theta_2=0.9
##goal caculate the mean of theta and use bootstrap to estimate the variance and bca

theta_1=0.1
theta_2=0.9
K=50
Y=rep(K,1)
n=c(25,25)
N1=rgeom(n[1],theta_1)+1
N2=rgeom(n[2],theta_2)+1
N=c(N1,N2)
true_mean_theta=0.5*(theta_1+theta_2)


###method I the stratified estimate##
strat_estimate_function=function(N,d){
  NN=N[d]
  result=mean(1/NN)
  return(result)
}
strat_estimate=strat_estimate_function(N,c(1:K))
strat_estimate_bootstrap=boot(data=N,statistic=strat_estimate_function,R=10000)
strat_estimate_variance=var(strat_estimate_bootstrap$t)
strat_estimate_bias=strat_estimate-true_mean_theta
strat_estimate_mse=strat_estimate_bias^2+strat_estimate_variance

###method II the pooled estimate##
pooled_estimate_function=function(N,d){
  NN=N[d]
  result=length(NN)/sum(NN)
  return(result)
}
pooled_estimate=pooled_estimate_function(N,c(1:K))
pooled_estimate_bootstrap=boot(N,statistic=pooled_estimate_function,R=10000)
pooled_estimate_variance=var(pooled_estimate_bootstrap$t)
pooled_estimate_bias=pooled_estimate-true_mean_theta
pooled_estimate_mse=pooled_estimate_bias^2+pooled_estimate_variance

##method III the NPML estimate##
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
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  return(sum(UU*pp))
}
npml_estimate=npml_estimate_function(N,c(1:K))
npml_estimate_bootstrap=boot(data=N,statistic=npml_estimate_function,R=10000)
npml_estimate_variance=var(npml_estimate_bootstrap$t)
npml_estimate_bias=npml_estimate-true_mean_theta
npml_estimate_mse=npml_estimate_bias^2+npml_estimate_variance

##method IV the test stratified estimate##
strat_estimate_function_test=function(N,d,n){
  NN=N[d]
  L=length(n)
  result=NULL
  for(i in 1:L){
    if(i==1){
      temp=n[i]/sum(NN[1:n[i]])
      result=c(result,temp)
    }
    else{
      temp=n[i]/sum(NN[n[i-1]+1:n[i]])
      result=c(result,temp)
    }
  }
  
  
  return(mean(result))
}
strat_estimate_function_test_estimate=strat_estimate_function_test(N,c(1:K),n)
strat_estimate_function_test_bias=strat_estimate_function_test_estimate-true_mean_theta
strat_estimate_function_test_bootstrap=boot(data=N,statistic=strat_estimate_function_test,R=10000,n=n)
strat_estimate_function_test_variance=var(strat_estimate_function_test_bootstrap$t)
strat_estimate_function_test_mse=strat_estimate_function_test_variance+strat_estimate_function_test_bias^2


###Case II###
##we set the sample size as 50. there are 10 of them from the theta_1=0.15 and 10 of them from the theta_2=0.35
## 10 of them from theta_3=0.55; 10 of them from the theta_4=0.75 and 10 of them from the theta_5=0.95
##goal caculate the mean of theta and use bootstrap to estimate the variance and bca

theta_1=0.15
theta_2=0.35
theta_3=0.55
theta_4=0.75
theta_5=0.95
K=50
Y=rep(K,1)
N1=rgeom(10,theta_1)+1
N2=rgeom(10,theta_2)+1
N3=rgeom(10,theta_3)+1
N4=rgeom(10,theta_4)+1
N5=rgeom(10,theta_5)+1
N=c(N1,N2,N3,N4,N5)
true_mean_theta=0.2*(theta_1+theta_2+theta_3+theta_4+theta_5)


###method I the stratified estimate##
strat_estimate_function=function(N,d){
  NN=N[d]
  result=mean(1/NN)
  return(result)
}
strat_estimate=strat_estimate_function(N,c(1:K))
strat_estimate_bootstrap=boot(data=N,statistic=strat_estimate_function,R=10000)
strat_estimate_variance=var(strat_estimate_bootstrap$t)
strat_estimate_bias=strat_estimate-true_mean_theta
strat_estimate_mse=strat_estimate_bias^2+strat_estimate_variance

###method II the pooled estimate##
pooled_estimate_function=function(N,d){
  NN=N[d]
  result=length(NN)/sum(NN)
  return(result)
}
pooled_estimate=pooled_estimate_function(N,K)
pooled_estimate_bootstrap=boot(N,statistic=pooled_estimate_function,R=10000)
pooled_estimate_variance=var(pooled_estimate_bootstrap$t)
pooled_estimate_bias=pooled_estimate-true_mean_theta
pooled_estimate_mse=pooled_estimate_bias^2+pooled_estimate_variance

##method III the NPML estimate##
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
    if(sum((temp_UU-UU)^2)<1e-5){
      break
    }
    
    temp_UU=UU
  }
  return(sum(UU*pp))
}
npml_estimate=npml_estimate_function(N,c(1:K))
npml_estimate_bootstrap=boot(data=N,statistic=npml_estimate_function,R=10000)
npml_estimate_variance=var(npml_estimate_bootstrap$t)
npml_estimate_bias=npml_estimate-true_mean_theta
npml_estimate_mse=npml_estimate_bias^2+npml_estimate_variance