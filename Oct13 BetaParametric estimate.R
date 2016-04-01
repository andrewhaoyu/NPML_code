rm(list=ls()) #clear the environment

commanarg <- commandArgs(trailingOnly = TRUE)
i1 <- as.numeric(commanarg[1])
i2 <- as.numeric(commanarg[2])

# method function ---------------------------------------------------------

library(bbmle)
BetaGeometricLikehood <- function(N,par.vec){
  ThetaBar.Par <- par.vec[1]
  M.Par <- par.vec[2]
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- log(alpha)+lgamma(alpha+beta)+lgamma(beta+N-1)-lgamma(beta)-lgamma(alpha+beta+N)
  return(-sum(result))
}



ParametricBeta <- function(Theta.Bar,M,N){
  theta.grid <- seq(0.7*Theta.Bar,1.3*Theta.Bar,by = (0.6*Theta.Bar)/1000)
  M.grid <- seq(0.7*M,1.3*M,by=(0.6*M)/1000)
  grid <- expand.grid(theta.grid,M.grid)
  myresult <- apply(grid,1,function(x) BetaGeometricLikehood(N,x))
  id <- which(myresult ==min(myresult))
  return(grid[id,])
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
Theta.Bar.Vector <- c(0.1)
M.Vector <- c(1)
K.Vector <- c(150)
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
        Para.Beta.Estimate.Result <- c(Para.Beta.Estimate.Result,unlist(ParametricBeta(True.Theta.Mean,M,N))[1])
      }
      
      
      Para.Beta.Estimate.Bias <- mean(Para.Beta.Estimate.Result)-True.Theta.Mean
      Para.Beta.Estimate.Var <- var(Para.Beta.Estimate.Result)
      Para.Beta.Estimate.mse <- Para.Beta.Estimate.Bias^2+ Para.Beta.Estimate.Var
      assign(paste("Para.Beta.Result_",Theta.Bar,"_",M,"_",K,sep=""),list(bias=round(Para.Beta.Estimate.Bias,places),var=round(Para.Beta.Estimate.Var,places),mse=round(Para.Beta.Estimate.mse,places)))
      print(c(i1,i2,i3))
      
      
      
    }
 



