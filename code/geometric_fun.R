
StratEstimateFunction <- function(N,d,cen) {
  NN <- N[d] #set d as parameters for convenience to use {boot} package later
  cend <- cen[d]
  MLE <- rep(0,length(NN))
  MLE[cend==1] <- 1/NN[cend==1]
  
  result <- mean(MLE)
  return(result)
}



#this is the function for the pooled estimate
PooledEstimateFunction <- function(N,d,cen) {
  NN <- N[d]
  K <- length(NN)
  cend <- cen[d]
  result <- sum(cend)/(sum(NN)-K+sum(cend))
  return(result)
}



#this recursion function is used to replace the loop inside NPML function
recursion <- function (pp, UU, NN, Y,cend) {
  ww <- pp*(1-UU)^(NN-Y) * UU^(Y*cend)
  return(ww)
}
#this is the function for the NPML estimate
NPMLEstimateFunction <- function(N,d,cen) {
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
  cend = cen[d]
  Y <- rep(1,KK)
  UU <- 1:KK/(KK+1)
  
  ww <- matrix(1, nrow = KK, ncol=KK)
  pp <- rep(1/KK, length = KK)
  Niter <- 500
  temp_UU=1
  epi=0.001
  temp_pp=1
  tol <- 1e-5
  # The recursion
  
  for(i in 1:Niter) {
    #Niter
    
    for(k in 1:KK) {
      ww[k,] <- mapply(function(x,y) recursion(x,y,NN[k],Y[k],cend[k]),pp,UU)
      ww[k,] <- ww[k,]/sum(ww[k,])
    }  #k to K
    pp <- apply(ww,2,sum)/sum(ww)
    for(j in 1:KK) {
      UU[j] <- sum(ww[,j]*Y*cend)/sum(ww[,j]*(NN-Y+cend))
    }
    UU[which(UU>(1-epi))] <- 1-epi
    UU[which(UU<epi)] <- epi
    if(any(is.nan(UU))==T) {
      break
    }
    if((sum((temp_UU-UU)^2)<tol)&sum((temp_pp-pp)^2)<tol) {
      break
    }
    temp_pp <- pp
    temp_UU <- UU
  }
  
  
  return(list(UU,pp))
  
}






#this recursion function is used to replace the loop inside NPML function

#this is the function for the NPML estimate


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

BetaGeometricLikehood <- function(x){
  
  ThetaBar.Par<- x[1]
  M.Par <- x[2]
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- lgamma(alpha+beta)+lgamma(alpha+cen)+lgamma(beta+N-1)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+N+cen-1)
  return(sum(result))
}
# n <- 1000
# begin_point1 <- runif(n,0.005,0.1)
# begin_point2 <- runif(n,150,250)
# begin_point <- cbind(begin_point1,begin_point2)
# try_result <- rep(0,n)
# 
# for(i in 1:n){
#   M <- 5
#   fit <-  optim(par = begin_point[i,],fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.8,10000),method="L-BFGS-B",control=list(fnscale=-1))
#   try_result[i] <- BetaGeometricLikehood(fit$par)
# 
# }
# idx <- which.max(try_result)




ParaBetaEstimateFunction <- function(N){
  fit <-  optim(par = begin,fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.5,10000),method="L-BFGS-B",control=list(fnscale=-1))
  
  return(fit$par[1])
}
LogL.Derivatives <- function(x){
  return(c(LogL.Dmu(x[1],x[2]),LogL.DM(x[1],x[2])))
}
LogL.Dmu <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par)*M.Par-M*LogL.Dbeta(ThetaBar.Par,M.Par)
  return(result)
}
LogL.DM <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par)*ThetaBar.Par-(1-ThetaBar.Par)*LogL.Dbeta(ThetaBar.Par,M.Par)
}
LogL.Dalpha <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/(exp(lgamma(alpha+beta)))+digamma(alpha+cen)/exp(lgamma(alpha+cen))-digamma(alpha)/exp(lgamma(alpha))
  -digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  return(sum(result))
}

LogL.Dbeta <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/exp(lgamma(alpha+beta))+digamma(beta+N-1)/exp(lgamma(beta+N-1))-
    digamma(beta)/exp(lgamma(beta))-digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  sum(result)
}


NPMLEstimateFunction_mean <- function(N,d,cen) {
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
  cend = cen[d]
  Y <- rep(1,KK)
  UU <- 1:KK/(KK+1)
  
  ww <- matrix(1, nrow = KK, ncol=KK)
  pp <- rep(1/KK, length = KK)
  Niter <- 500
  temp_UU=1
  epi=0.001
  temp_pp=1
  tol <- 1e-5
  # The recursion
  
  for(i in 1:Niter) {
    #Niter
    
    for(k in 1:KK) {
      ww[k,] <- mapply(function(x,y) recursion(x,y,NN[k],Y[k],cend[k]),pp,UU)
      ww[k,] <- ww[k,]/sum(ww[k,])
    }  #k to K
    pp <- apply(ww,2,sum)/sum(ww)
    for(j in 1:KK) {
      UU[j] <- sum(ww[,j]*Y*cend)/sum(ww[,j]*(NN-Y+cend))
    }
    UU[which(UU>(1-epi))] <- 1-epi
    UU[which(UU<epi)] <- epi
    if(any(is.nan(UU))==T) {
      break
    }
    if((sum((temp_UU-UU)^2)<tol)&sum((temp_pp-pp)^2)<tol) {
      break
    }
    temp_pp <- pp
    temp_UU <- UU
  }
  
  
  
  
  return(sum(UU*pp))
  
  
}
