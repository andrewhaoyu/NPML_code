#recursion function is part of the npml algorithm
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
  Niter=10000
  temp_UU=1
  epi=0.001
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
    UU[which(UU==1)]=1-epi
    UU[which(UU==0)]=epi
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


for(i1 in 3:3){
  for(i2 in 1:length(M_vector)){
    for(i3 in 1:length(K_vector)){
      print(xtable(get(paste("true_",theta_bar_vector[i1],"_",M_vector[i2],"_",K_vector[i3],sep="")),caption=paste("theta=",theta_bar_vector[i1],",","M=",M_vector[i2],",","K=",K_vector[i3],sep=""),digits=4))     
      print(xtable(get(paste("result_",theta_bar_vector[i1],"_",M_vector[i2],"_",K_vector[i3],sep="")),caption=paste("theta=",theta_bar_vector[i1],",","M=",M_vector[i2],",","K=",K_vector[i3],sep=""),digits=4))     
      
    }
  }
}
print( xtable(get(paste("true_",theta_bar_vector[i1],"_",M_vector[i2],"_",K_vector[i3],sep="")),digits=4) )



estimate=c("E(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,"logit(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"E(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,"logit(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"E(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,"logit(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"E(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,"logit(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"E(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA,"logit(P)",NA,NA,NA,NA,NA,NA,NA,NA,NA)
theta_0.1=c("M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","theta=0.2","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","theta=0.3","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","theta=0.4","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","theta=0.5","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000","M=0.05","M=0.25","M=0.5","M=1","M=1.5","M=2","M=5","M=10","M=25","M=10000")                         

K_150=K_100=K_50=K_20=K_10=rep(NA,104)
for(i3 in 1:length(K_vector)){
  K_result=rep(NA,104)
  for(i1 in 1:length(theta_bar_vector)){
    for(i2 in 1:length(M_vector)){
      temp=get(paste("result_",theta_bar_vector[i1],"_",M_vector[i2],"_",K_vector[i3],sep=""))
      
        K_result[i2+(i1-1)*20+(i1-1)]=compare_mean(temp)
        K_result[i2+(i1-1)*20+(i1-1)+10]=compare_logit(temp)
      
      
      
    }
  }
  assign(paste("K_",K_vector[i3],sep=""),K_result)
}
Temp_result=data.frame(estimate,theta_0.1,K_10,K_20,K_50,K_100)




xtable(Temp_result)








save(list=objectname[1],file=paste("~/Documents/study/project/cristian/project1/result/",objectname[1],".Rdata",sep=""))

save(list=objectname[2],file=paste("~/Documents/study/project/cristian/project1/result/",objectname[2],".Rdata",sep=""))


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






plot the mse trend ------------------------------------------------------
  # plot the mse trend with three different styles:
  #   K and M fixed, theta grows
  #   K and theta fixed, M grows
  #   M and theta fixed, K grows
  pdf('K_M_fixed_theta_grow_mean', width = 8.3, height = 11.7)  ## Device with dimensions of A4 paper
par(omi = rep(.5, 4))                      ## 1/2 inch outer margins
par(mfrow = c(5,3)) 
##K and M fixed, theta grows
for(i1 in 1:length(K.Vector)){
  for(i2 in 1:(length(M.Vector)-1)){
    jdx=which(colnames(Temp.Result)==paste("K_",K.Vector[i1],sep=""))
    idx=which(Temp.Result[,2]==paste("M=",M.Vector[i2],sep=""))
    plot(Temp.Result[idx[(2*(1:length(K.Vector))-1)],jdx],ylab="relative_mse_mean",main=paste("K_",K.Vector[i1],",M_",M.Vector[i2],",theta_grow,"))
    abline(h=100,col="red")
  }
}
for(i1 in 1:length(Theta.Bar.Vector)){
  for(i2 in 1:(length(M.Vector)-1)){
    jdx <- which(colnames(Temp.Result)==paste("K_",K.Vector[i1],sep=""))
    idx <- which(Temp.Result[,2]==paste("M=",M.Vector[i2],sep=""))
    plot(Temp.Result[idx[(2*(1:length(K.Vector)))],jdx],ylab="relative_mse_logitmean",main=paste("K_",K.Vector[i1],",M_",M.Vector[i2],",theta_grow,"))
    abline(h=100,col="red")
  }
}
##K and theta fixed,M grows
for(i1 in 1:length(Theta.Bar.Vector)){
  for(i2 in 1:length(K.Vector)){
    plot(Temp.Result[(i1+2*length(M.Vector)*(i1-1)):(i1-1+length(M.Vector)-1+2*length(M.Vector)*(i1-1)),i2+2],ylab="realative_mse_mean",main=paste("theta_",Theta.Bar.Vector[i1],",K_",K.Vector[i2],",M_grows"))
    abline(h=100,col="red")
  }
}
for(i1 in 1:length(Theta.Bar.Vector)){
  for(i2 in 1:length(K.Vector)){
    plot(Temp.Result[(i1+2*length(M.Vector)*(i1-1)+length(M.Vector)):(i1-1+length(M.Vector)-1+2*length(M.Vector)*(i1-1)+length(M.Vector)),i2+2],ylab="realative_mse_logitmean",main=paste("theta_",Theta.Bar.Vector[i1],",K_",K.Vector[i2],",M_grows"))
    abline(h=100,col="red")
  }
}
##theta and M fixed, K grows
for(i1 in 1:length(Theta.Bar.Vector)){
  for(i2 in 1:(length(M.Vector)-1)){
    plot(as.numeric(Temp.Result[((2*length(M.Vector)+1)*(i1-1)+i2),3:(2+length(K.Vector))]),ylab="relative_mse_mean",main=paste("theta_",Theta.Bar.Vector[i1],",M_",M.Vector[i2],",K_grows"))
    abline(h=100,col="red")
  }
}
for(i1 in 1:length(Theta.Bar.Vector)){
  for(i2 in 1:(length(M.Vector)-1)){
    plot(as.numeric(Temp.Result[((2*length(M.Vector)+1)*(i1-1)+i2+length(M.Vector)),3:(2+length(K.Vector))]),ylab="relative_mse_logitmean",main=paste("theta_",Theta.Bar.Vector[i1],",M_",M.Vector[i2],",K_grows"))
    abline(h=100,col="red")
  }
}
dev.off() #save the pdf into the local directory

