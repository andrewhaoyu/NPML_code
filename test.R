rm(list=ls()) #clear the environment

commanarg <- commandArgs(trailingOnly = TRUE)
i1 <- as.numeric(commanarg[1])
i2 <- as.numeric(commanarg[2])

# method function ---------------------------------------------------------

library(bbmle)
BetaGeometricLikehood <- function(Logit.ThetaBar.Par,Log.M.Par){
  a <- Logit.ThetaBar.Par
  b <- Log.M.Par
  ThetaBar.Par <- exp(a)/(exp(a)+1)
  M.Par <- exp(b)
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- log(alpha)+lgamma(alpha+beta)+lgamma(beta+N-1)-lgamma(beta)-lgamma(alpha+beta+N)
  return(-sum(result))
}
ParaBetaEstimateFunction <- function(N){
  fit <-  mle2(BetaGeometricLikehood,method="SANN",start = list(Logit.ThetaBar.Par=logit(True.Theta.Mean),Log.M.Par=M))
  a<-coef(fit)[[1]]
  return(exp(a)/(exp(a)+1))
}

#this is the function for the stratified estimate
StratEstimateFunction <- function(N,d) {
  NN <- N[d] #set d as parameters for convenience to use {boot} package later
  result <- mean(1/NN)
  return(result)
}



#this is the function for the pooled estimate
PooledEstimateFunction <- function(N,d) {
  NN <- N[d]
  result <- length(NN)/sum(NN)
  return(result)
}



#this recursion function is used to replace the loop inside NPML function
recursion <- function (pp, UU, NN, Y) {
  ww <- pp*(1-UU)^(NN-Y) * UU^Y
  return(ww)
}
#this is the function for the NPML estimate
NPMLEstimateFunction <- function(N,d) {
  #beginning setting of NPML:
  # NN is the geometric sampling result
  # UU is the probability mass point of NPML
  # UU begins with 1:K/(K+1)
  # ww is the posterior probability distribution of ith sample coming from UU_jth
  # YY is the a vector with all elements equals 1,
  # YY setting is easy for negative bionomial case
  # pp is the probability distrubution of UU
  #recursion steps of NPML:
  # 500 Niter steps for the EM algorithm
  # a standard to break the loop,
  # if sum(UU_nextstep-UU_now)^2<1e-5&sum(ww_nextstep-ww_now)^2<1e-5&steps>25
  # set a epi standard for UU, if UU>(1-epi),UU=1-epi; if UU<epi, UU=epi;
  #results of the NPML:
  # we return a lis of UU and pp for the NPML
  NN <- N[d] 
  KK <- length(NN)  
  Y <- rep(1,KK)
  UU <- 1:KK/(KK+1) 
  
  ww <- matrix(1, nrow = KK, ncol=KK)
  pp <- rep(1/KK, length = KK)
  Niter <- 500
  temp_UU=1
  epi=0.001
  temp_pp=1
  # The recursion
  
  for(i in 1:Niter) {   
    #Niter
    
    for(k in 1:KK) {
      ww[k,] <- mapply(function(x,y) recursion(x,y,NN[k],Y[k]),pp,UU)
      ww[k,] <- ww[k,]/sum(ww[k,])
    }  #k to K
    pp <- apply(ww,2,sum)/sum(ww)
    for(j in 1:KK) {
      UU[j] <- sum(ww[,j]*Y)/sum(ww[,j]*NN)
    }
    UU[which(UU>(1-epi))] <- 1-epi
    UU[which(UU<epi)] <- epi
    if(any(is.nan(UU))==T) {
      break
    }
    if((sum((temp_UU-UU)^2)<1e-5)&sum((temp_pp-pp)^2)<1e-5&(i>50)) {
      break
    }
    temp_pp <- pp
    temp_UU <- UU
  }
  
  
  return(list(UU,pp))
  
}



#caculate the logit transform
logit <- function(x) {
  return (log(x/(1-x)))
}



#compare the mse of the mean estimate results of the three methods
#SP means Strat and Pooled is same, SN means Strat and NPML same, PN means Pooled and NPML same
CompareMean <- function(x) {
  if(x[4,2]<=x[4,3]&x[4,2]<=x[4,4]) {
    if(x[4,2]==x[4,3]){
      return("SP")
    }
    if(x[4,2]==x[4,4]){
      return("SN")
    }
    return("Strat")
  }
  
  if(x[4,3]<=x[4,2]&x[4,3]<=x[4,4]) {
    if(x[4,3]==x[4,2]) {
      return("SP")
    }
    if(x[4,3]==x[4,4]) {
      return("PN")
    }
    return("Pooled")
  }
  
  if(x[4,4]<=x[4,2]&x[4,4]<=x[4,3]) {
    if(x[4,4]==x[4,2]) {
      return("SN")
    }
    if(x[4,4]==x[4,3]) {
      return("PN")
    }
    return("Npml")
  }
  
}



#compare the mse of estimate of logit transform of three methods
CompareLogit <- function(x) {
  if(x[12,2]<=x[12,3]&x[12,2]<=x[12,4]) {
    if(x[12,2]==x[12,3]) {
      return("SP")
    }
    if(x[12,2]==x[12,4]) {
      return("SN")
    }
    return("Strat")
  }
  
  if(x[12,3]<=x[12,2]&x[12,3]<=x[12,4]) {
    if(x[12,3]==x[12,2]) {
      return("SP")
    }
    if(x[12,3]==x[12,4]) {
      return("PN")
    }
    return("Pooled")
  }
  
  if(x[12,4]<=x[12,2]&x[12,4]<=x[12,3]) {
    if(x[12,4]==x[12,2]) {
      return("SN")
    }
    if(x[12,4]==x[12,3]) {
      return("PN")
    }
    return("Npml")
  }
  
  
}



#compare the relative mse of the three different methods
RelativeMeanFunction <- function(x) {
  places <- 2
  temp <- 100*x[4,4]/x[4,3]
  return(round(temp,places))
}
RelativeLogitFunction <- function(x) {
  places <- 2
  temp <- 100*x[12,4]/x[12,3]
  return(round(temp,places))
}



#Caculate the standard devivation of the simulation bias
SimulationSdBiasMean <- function(R,x) {
  places <- 4
  return(round(sqrt(x[3,4]/R),places))
}
SimulationSdBiasLogit <- function(R,x) {
  places <- 4
  return(round(sqrt(x[11,4]/R),places))
}





# calculate three methods results -----------------------------------------------
#setting of the true parameters:
# Theta.Bar.Vector(mean),M.Vector(precision),K.Vector(sample size) R(replicates)
# Theta.Low and Theta.High are the boundaries of theta,
# if theta is less than Theta.Low, then theta equals Theta.Low; same case for Theta.High
Theta.Bar.Vector <- c(0.1,0.2,0.3,0.4,0.5)
M.Vector <- c(0.05,0.25,0.5,1,1.5,2,5,10,25,10000)
K.Vector <- c(10,20,50,75,100,150)
Theta.Low <- 0.005
Theta.High <- 0.995

    for(i3 in 1:length(K.Vector)) {
      
      
      
      places <- 5
      R <- 500
      #get alpha and beta from tehta_bar and M
      Theta.Bar <- Theta.Bar.Vector[i1]
      M <- M.Vector[i2]
      alpha <- Theta.Bar*M
      beta <- M*(1-Theta.Bar)
      K <- K.Vector[i3]
      #choose the 1:K/(K+1) quantile of beta distribution
      cuts <- 1:K/(K+1)
      theta <- qbeta(cuts,alpha,beta)
      theta[which(theta<Theta.Low)] <- Theta.Low
      theta[which(theta>Theta.High)] <- Theta.High
      #calculate the true results of theta mean, theta recip and theta logit
      True.Theta.Mean <- mean(theta)
      True.Theta.Recip <- mean(1/theta)
      True.Theta.Logit <- mean(logit(theta))
      True.Theta.Var <- var(theta)*(K-1)/K
      True.Theta.Recip.Var <- var(1/theta)*(K-1)/K
      True.Theta.Logit.Var <- var(logit(theta))*(K-1)/K
      
      
      
      #get the esimate of the three different methods
      Stratified.Estimate.Results <- NULL
      Pooled.Estimate.Results <- NULL
      Npml.Estimate.UU <- NULL
      Npml.Estimate.pp <- NULL
      Para.Beta.Estimate.Result <- NULL
      set.seed(123456)
      NTL <- sapply(theta,function(x) (rgeom(R,x)+1)) #generate the total N for all R replicates
      
      
      for(ind in 1:R){
       
        N <- NTL[ind,]
        #apply the three method function to the data
        Para.Beta.Estimate.Result <- c(Para.Beta.Estimate.Result,ParaBetaEstimateFunction(N))
        
        Stratified.Estimate.Results <- c(Stratified.Estimate.Results,StratEstimateFunction(N,c(1:K)))
        Pooled.Estimate.Results <- c(Pooled.Estimate.Results,PooledEstimateFunction(N,c(1:K)))
        temp <- NPMLEstimateFunction(N,c(1:K))
        Npml.Estimate.UU <- cbind(Npml.Estimate.UU,temp[[1]])
        Npml.Estimate.pp <- cbind(Npml.Estimate.pp,temp[[2]])
        
        
        
        
        
        
        
        
      }
      
      
      
      
      #caculate the mean estimate of the three different methods
      #caculate the bias and variance of the mean estimate and get the mse
      Stratified.Estimate.Bias <- mean(Stratified.Estimate.Results)-True.Theta.Mean
      Stratified.Estimate.Var <- var(Stratified.Estimate.Results)
      Stratified.Estimate.Mse <- Stratified.Estimate.Bias^2+Stratified.Estimate.Var
      Pooled.Estimate.Bias <- mean(Pooled.Estimate.Results)-True.Theta.Mean
      Pooled.Estimate.Var <- var(Pooled.Estimate.Results)
      Poole.Estimate.Mse <- Pooled.Estimate.Bias^2+Pooled.Estimate.Var
      Npml.Estimate.Results <- apply(Npml.Estimate.pp*Npml.Estimate.UU,2,sum)
      Npml.Estimate.Bias <- mean(Npml.Estimate.Results)-True.Theta.Mean
      Npml.Estimate.Var <- var(Npml.Estimate.Results)
      Npml.Estimate.Mse <- Npml.Estimate.Bias^2+Npml.Estimate.Var
      Result.Mean <- list(True.Theta.Mean=round(True.Theta.Mean,places),True.Theta.Var=round(True.Theta.Var,places),Stratified.Estimate.Bias=round(Stratified.Estimate.Bias,places),Stratified.Estimate.Var=round(Stratified.Estimate.Var,places),Stratified.Estimate.Mse=round(Stratified.Estimate.Mse,places),Pooled.Estimate.Bias=round(Pooled.Estimate.Bias,places),Pooled.Estimate.Var=round(Pooled.Estimate.Var,places),Poole.Estimate.Mse=round(Poole.Estimate.Mse,places),Npml.Estimate.Bias=round(Npml.Estimate.Bias,places),Npml.Estimate.Var=round(Npml.Estimate.Var,places),Npml.Estimate.Mse=round(Npml.Estimate.Mse,places))
      
      
      
      #caculate the recip estimate of the three defferent methods
      #caculate the bias and variance of the recip estimate and get the mse
      Stratified.Estimate.Recip.Bias <- mean(1/Stratified.Estimate.Results)-True.Theta.Recip
      Stratified.Estimate.Recip.Var <- var(1/Stratified.Estimate.Results)
      Stratified.Estimate.Recip.Mse <- Stratified.Estimate.Recip.Bias^2+ Stratified.Estimate.Recip.Var
      Pooled.Estimate.Recip.Bias <- mean(1/Pooled.Estimate.Results)-True.Theta.Recip
      Pooled.Estimate.Recip.Var <- var(1/Pooled.Estimate.Results)
      Pooled.Estimate.Recip.Mse <- Pooled.Estimate.Recip.Bias^2+Pooled.Estimate.Recip.Var
      Npml.Estimate.Recip.Results <- apply(Npml.Estimate.pp*(1/Npml.Estimate.UU),2,sum)
      Npml.Estimate.Recip.Bias <- mean(Npml.Estimate.Recip.Results)-True.Theta.Recip
      Npml.Estimate.Recip.Var <- var(Npml.Estimate.Recip.Results)
      Npml.Estimate.Recip.Mse <- Npml.Estimate.Recip.Bias^2+Npml.Estimate.Recip.Var
      Result.Recip <- list(True.Theta.Recip=round(True.Theta.Recip,places),True.Theta.Recip.Var=round(True.Theta.Recip.Var,places),strat_estimate_recip_bias=round(Stratified.Estimate.Recip.Bias,places), Stratified.Estimate.Recip.Var=round( Stratified.Estimate.Recip.Var,places),Stratified.Estimate.Recip.Mse=round(Stratified.Estimate.Recip.Mse,places),Pooled.Estimate.Recip.Bias=round(Pooled.Estimate.Recip.Bias,places),Pooled.Estimate.Recip.Var=round(Pooled.Estimate.Recip.Var,places),Pooled.Estimate.Recip.Mse=round(Pooled.Estimate.Recip.Mse,places),Npml.Estimate.Recip.Bias=round(Npml.Estimate.Recip.Bias,places),Npml.Estimate.Recip.Var=round(Npml.Estimate.Recip.Var,places),Npml.Estimate.Recip.Mse=round(Npml.Estimate.Recip.Mse,places))
      
      
      
      #caculate the logit estimate of the three different methods
      #caculate the bias and variance of the logit estimate and get the mse
      temp <- (Npml.Estimate.UU==1)
      New.Npml.Estimate.UU <- NULL
      New.Npml.Estimate.pp <- NULL
      for(i in 1:ncol(temp)){
        if(sum(temp[,i]==F)==K){
          New.Npml.Estimate.UU <- cbind(New.Npml.Estimate.UU,Npml.Estimate.UU[,i])
          New.Npml.Estimate.pp <- cbind(New.Npml.Estimate.pp,Npml.Estimate.pp[,i])
        }
      }
      Stratified.Estimate.Logit.Bias <- mean(logit(Stratified.Estimate.Results))-True.Theta.Logit
      Stratified.Estimate.Logit.Var <- var(logit(Stratified.Estimate.Results))
      Stratified.Estimate.Logit.Mse <- Stratified.Estimate.Logit.Bias^2+Stratified.Estimate.Logit.Var
      Pooled.Estimate.Logit.Bias <- mean(logit(Pooled.Estimate.Results))-True.Theta.Logit
      Pooled.Estimate.Logit.Var <- var(logit(Pooled.Estimate.Results))
      Pooled.Estimate.Logit.Mse <- Pooled.Estimate.Logit.Bias^2+Pooled.Estimate.Logit.Var
      Npml.Estimate.Logit.Results <- apply(New.Npml.Estimate.pp*(logit(New.Npml.Estimate.UU)),2,sum)
      Npml.Estimate.Logit.Bias <- mean(Npml.Estimate.Logit.Results)-True.Theta.Logit
      Npml.Estimate.Logit.Var <- var(Npml.Estimate.Logit.Results)
      Npml.Estimate.Logit.Mse <- Npml.Estimate.Logit.Bias^2+Npml.Estimate.Logit.Var
      Result.Logit <- list(True.Theta.Logit=round(True.Theta.Logit,places),True.Theta.Logit.Var=round(True.Theta.Logit.Var,places),strat_estimate_logit_bias=round(Stratified.Estimate.Logit.Bias,places),Stratified.Estimate.Logit.Var=round(Stratified.Estimate.Logit.Var,places),Stratified.Estimate.Logit.Mse=round(Stratified.Estimate.Logit.Mse,places),Pooled.Estimate.Logit.Bias=round(Pooled.Estimate.Logit.Bias,places),Pooled.Estimate.Logit.Var=round(Pooled.Estimate.Logit.Var,places),Pooled.Estimate.Logit.Mse=round(Pooled.Estimate.Logit.Mse,places),Npml.Estimate.Logit.Bias=round(Npml.Estimate.Logit.Bias,places),Npml.Estimate.Logit.Var=round(Npml.Estimate.Logit.Var,places),Npml.Estimate.Logit.Mse=round(Npml.Estimate.Logit.Mse,places))
      
      
      
      #build a data frame for the result
      assign(paste("result_",Theta.Bar,"_",M,"_",K,"_list",sep=""),list(Result.Mean,Result.Recip,Result.Logit))
      assign(paste("result"),list(Result.Mean,Result.Recip,Result.Logit))
      method <- c("E(P)","Bias","Var","MSE","E(1/P)","Bias","Var","MSE","E(logit(P))","Bias","Var","MSE")
      strat <- c(NA,as.numeric(result[[1]][3]),as.numeric(result[[1]][4]),as.numeric(result[[1]][5]),NA,as.numeric(result[[2]][3]),as.numeric(result[[2]][4]),as.numeric(result[[2]][5]),NA,as.numeric(result[[3]][3]),as.numeric(result[[3]][4]),as.numeric(result[[3]][5]))
      pooled <- c(NA,as.numeric(result[[1]][6]),as.numeric(result[[1]][7]),as.numeric(result[[1]][8]),NA,as.numeric(result[[2]][6]),as.numeric(result[[2]][7]),as.numeric(result[[2]][8]),NA,as.numeric(result[[3]][6]),as.numeric(result[[3]][7]),as.numeric(result[[3]][8]))
      npml <- c(NA,as.numeric(result[[1]][9]),as.numeric(result[[1]][10]),as.numeric(result[[1]][11]),NA,as.numeric(result[[2]][9]),as.numeric(result[[2]][10]),as.numeric(result[[2]][11]),NA,as.numeric(result[[3]][9]),as.numeric(result[[3]][10]),as.numeric(result[[3]][11]))
      True.Value <- c("E","Var")
      P <- c(as.numeric(result[[1]][1]),as.numeric(result[[1]][2]))
      Recip.P <- c(as.numeric(result[[2]][1]),as.numeric(result[[2]][2]))
      Logit.P <- c(as.numeric(result[[3]][1]),as.numeric(result[[3]][2]))
      assign(paste("result_",Theta.Bar,"_",M,"_",K,sep=""),data.frame(method,strat,pooled,npml))
      assign(paste("true_",Theta.Bar,"_",M,"_",K,sep=""),data.frame(True.Value,P,Recip.P,Logit.P))
      objectname <- c(paste("result_",Theta.Bar,"_",M,"_",K,sep=""),paste("true_",Theta.Bar,"_",M,"_",K,sep=""))
      Para.Beta.Estimate.Bias <- mean(Para.Beta.Estimate.Result)-True.Theta.Mean
      Para.Beta.Estimate.Var <- var(Para.Beta.Estimate.Result)
      Para.Beta.Estimate.mse <- Para.Beta.Estimate.Bias^2+ Para.Beta.Estimate.Var
      assign(paste("Para.Beta.Result_",Theta.Bar,"_",M,"_",K,sep=""),list(bias=round(Para.Beta.Estimate.Bias,places),var=round(Para.Beta.Estimate.Var,places),mse=round(Para.Beta.Estimate.mse,places)))
      print(c(i1,i2,i3))
      
      
      
    }
 



save.image(paste0("/home/student/hzhang1/R/project_Tom/Sep_30/",Theta.Bar.Vector[i1],"_",M.Vector[i2]))