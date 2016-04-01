rm(list=ls()) #clear the environment
commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])
i2 <- as.numeric(commanarg[2])
i3 <- as.numeric(commanarg[3])
i4 <- as.numeric(commanarg[4])
# method function ---------------------------------------------------------
library(MASS)
library(nleqslv)
logit_inver <- function(x){
  exp(x)/(1+exp(x))
}
gr_u <- function(u,N,ww_vec,b,x){
  K <- length(N)
  sum <- 0
  temp_Qu <- 1-N*logit_inver(u+b*x)
  sum <- crossprod(temp_Qu,ww_vec)
  return(sum^2)
}
gr_u_d <- function(u,N,ww_vec,b,x){
  K <- length(N)
  sum <- 0
  temp_Qu <- 1-N*logit_inver(u+b*x)
  sum_1 <- crossprod(temp_Qu,ww_vec)
  temp_d_Qu <- N*logit_inver(u+b*x)/(1+exp(u+b*x))
  sum_2 <- -crossprod(temp_d_Qu,ww_vec)
  return(sum_1*sum_2)
}
gr_b <- function(uu,N,ww,b,x){
  u <- uu
  K <- length(N)
  t <- 0
  temp_Qb <- matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      temp_Qb[i,j] <- x[i]-N[i]*x[i]*logit_inver(u[j]+b*x[i])
    }
  }
  t <- sum(temp_Qb*ww)
  return(t^2)
}
NPML_logistic_function <- function(N,x){
  w <- rep(1/K,K)
  uu <- runif(K,u_new[1],u_new[K])
  b <- 0
  pp <- matrix(0,K,K)
  ww <- matrix(1,K,K)
  epi <- 0.001
  Niter <- 1000
  idx <- 0
  try <- NULL
  for(m in 1:Niter){
    
    b_old <- b
    uu_old <- uu
    #Estep
    for(i in 1:K){
      #Bayesian formula
      for(j in 1:K){
        pp[i,j] <- (1-logit_inver(uu[j]+b*x[i]))^(N[i]-1)*logit_inver(uu[j]+b*x[i])
        ww[i,j] <- w[j]*pp[i,j]
      }
      ww[i,] <- ww[i,]/sum(ww[i,])
     
    }
    for(i in 1:K){
      w[i] <- sum(ww[,i])/sum(ww)
    }
    #Mstep
    Niter_M <- 50
    epi_M <- 0.01
    jdx <- 0
    
    for(i in 1:Niter_M){
      uu_temp <-uu
      b_temp <- b
      for(j in 1:K){
        ww_vec <- ww[,j]

        ans <-  optim(u_new[j],fn=gr_u,gr=gr_u_d,method ="L-BFGS-B",lower=u_new_min,upper=u_new_max,N=N,ww_vec=ww_vec,b=b,x=x)
        uu[j] <- ans$par
      }
      ans <- nleqslv(b,fn=function(y){gr_b(uu,N,ww,y,x)})
      b <- ans$x
      
      thres_M <- sum((uu-uu_temp)^2)+(b-b_temp)^2
      if(thres_M < epi_M){
        break
      }
      jdx <- jdx+1
      
    }
    
    thres <- sum((uu_old-uu)^2)+(b-b_old)^2
    try <- c(try,thres)
    if(thres < epi){
      break
    }
    idx <- idx +1 
    
  }
return(c(crossprod(uu,w),b))
}




likelihood <- function(ww,N,x,par){
  sum <- 0
  K <- length(N)
  u <- par[1:K]
  b <- par[K+1]
  for(i in 1:K){
    for(j in 1:K){
      temp <- (u[j]+b*x[i]-N[i]*log(1+exp(u[j]+b*x[i])))*ww[i,j]
      sum <- temp+sum
    }
  }
  return(-sum)
}


standard_logistic <- function(N,x){
  
  data1 <- data.frame(x,N)
  model <- glm(cbind(1,N-1)~x,family = binomial(logit),data=data1)
  return(model$coefficients)
}
theta_set <- function(mu_u,sd_u,K,r,beta){
  z_vec <- qnorm((1:K)/(K+1))
  u_vec <- (z_vec-mean(z_vec))/sd(z_vec)
  u_new <- u_vec*sd_u+mu_u
  if(r==0){
    x <- rep(c(-1,1,1,-1,1,-1,-1,1),K/8)
    theta <- logit_inver(u_new+x*beta)
  }
  if(r==0.4){
    x <-  c(rep(-1,K/8),rep(c(-1,1,1,-1),3*K/16),rep(1,K/8))
    theta <- logit_inver(u_new+x*beta)
  }
  if(r==0.8){
    x <- c(rep(-1,K/2),rep(1,K/2))
    theta <- logit_inver(u_new+x*beta)
  }
  return(c(theta,u_new,x))
}
R <- 200

K_vector <- c(16,32,48,64,80,96,112)
sd_u_vector <- c(0,1,3,5)
r_vector <- c(0,0.4,0.8)
beta_vector <- c(0,1,2)
    set.seed(123)
    K <- K_vector[i1]
    sd_u <- sd_u_vector[i2]
    r <- r_vector[i3]
    beta <- beta_vector[i4]
    mu_u <- 0
    data <- theta_set(mu_u,sd_u,K,r,beta)
    theta <- data[1:K]
    u_new <- data[(K+1):(2*K)]
    u_new_min <- min(u_new)
    u_new_max <- max(u_new)
    x <- data[(2*K+1):(3*K)]
    NTL <- sapply(theta,function(x) (rgeom(R,x)+1)) #generate the total N for all R replicates

    standard_logistic_results <- apply(NTL,1,function(y){standard_logistic(y,x)})
    true_u <- mean(u_new)
    standard_logistic_u_mse <- mean((standard_logistic_results[1,]-true_u)^2)
    standard_logistic_b_mse <- mean((standard_logistic_results[2,]-beta)^2)
   NPML_logistic_results <- apply(NTL,1,function(y){NPML_logistic_function(y,x)})
    NPML_logistic_u_mse <- mean((NPML_logistic_results[1,]-true_u)^2)
    NPML_logistic_b_mse <- mean((NPML_logistic_results[2,]-beta)^2)
    mse_results <- c(standard_logistic_u_mse,standard_logistic_b_mse,NPML_logistic_u_mse,NPML_logistic_b_mse,K,sd_u,r,beta)
    save(mse_results,file=paste0("/home/student/hzhang1/R/project_Tom/logistic_regression/simulation2/mse_",K,"_",sd_u,"_",r,"_",beta))





