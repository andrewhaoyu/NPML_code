logit_inver <- function(x){
  exp(x)/(1+exp(x))
}
logit <- function(x){
  log(x/(1-x))
}
#####logistic regression
# standard_logistic <- function(N,x){
# 
#   data1 <- data.frame(x,N)
#   model <- glm(cbind(1,N-1)~x,family = binomial(logit),data=data1,epsilon=1e-14,maxit=1000)
#   return(model$coefficients)
# }
gr_u_logistic <- function(u,N,b,x){a
  K <- length(N)
  sum1 <- sum(1-N*logit_inver(u+b*x))
  return(sum1)
}

gr_b_logistic <- function(u,N,b,x){
  
  K <- length(N)
  sum1 <- sum(x-N*x*logit_inver(u+b*x))
  return(sum1)
}
logistic_lik <- function(u,N,b,x){
  sum1 <- sum(u+b*x-n*log(1+exp(u+b*x)))
  return(sum)
} 

newton_inter <- function(param,x,A,z){
  b <- crossprod(x,z)
  M <- wcrossprod(x,x,A)
  if(det(M)==0){
    param.new <- rep(NaN,length(param))
    return(param.new)
  }
  temp <- solve(M,b)
  
  param.new <- temp+param
  return(param.new)
}


standard_logistic <- function(N,x,param.start,tol,maxit){
  converge <- F
  beta_old <- param.start
  K <- length(N)
  x <- cbind(rep(1,K),x)
  predictor <- x%*%beta_old
  A <- diag(1/(1+exp(predictor)))
  z <- 1-N*logit_inver(predictor)
  
  for(i in 1:maxit){
    
    
    beta_new <- newton_inter(beta_old,x,A,z)
    
    #threshold <- sum((beta_new-beta_old)^2)
    if(sum(is.nan(beta_new))!=0|abs(sum(beta_new))>100){
      beta_new <- beta_old
      break
    }
    threshold1 <- sum(abs(beta_new-beta_old)/(abs(beta_new)+0.1))
    threshold2 <- sum((beta_new-beta_old)^2)
    threshold <- max(threshold1,threshold2)
    if(threshold < tol){
      converge <- T
      break
    }
    
    beta_old <- beta_new
    predictor <- x%*%beta_old
    A <- diag(1/(1+exp(predictor)))
    z <- 1-N*logit_inver(predictor)
  }
  
  return(c(beta_new,converge))
}

#####logistic regression warning record
# standard_logistic_warnings <- function(N,x){
# #   tryCatch(standard_logistic(N,x),warning=function(w){return(c(NA,NA))})
# # }

######generate mechanism
generating <- function(mu_u,sd_u,K,r,beta,R,observation){
  z_vec <- qnorm((1:K)/(K+1))
  u_vec <- (z_vec-mean(z_vec))/sd(z_vec)
  u_new <- u_vec*sd_u+mu_u
  u_new[u_new<(u_new_min)] <- u_new_min
  u_new[u_new>u_new_max] <- u_new_max
  if(observation==1){
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
    NTL <- sapply(theta,function(x) (rgeom(R,x)+1)) #generate the total N for all R replicates
    result <- list(u_new = u_new,x=x,theta=theta,NTL=NTL)
    return(result)
  }
  if(r==0){
    x <- rep(c(-1,1,1,-1,1,-1,-1,1),K/8)
    x1 <- x-0.5
    x2 <- x+0.5
    theta1 <- logit_inver(u_new+x1*beta)
    theta2 <- logit_inver(u_new+x2*beta)
  }
  if(r==0.4){
    x <-  c(rep(-1,K/8),rep(c(-1,1,1,-1),3*K/16),rep(1,K/8))
    x1 <- x-0.5
    x2 <- x+0.5
    theta1 <- logit_inver(u_new+x1*beta)
    theta2 <- logit_inver(u_new+x2*beta)
  }
  if(r==0.8){
    x <- c(rep(-1,K/2),rep(1,K/2))
    x1 <- x-0.5
    x2 <- x+0.5
    theta1 <- logit_inver(u_new+x1*beta)
    theta2 <- logit_inver(u_new+x2*beta)
  }
  x <- as.vector(rbind(x1,x2))
  theta <- as.vector(rbind(theta1,theta2))
  NTL <- sapply(theta,function(x) (rgeom(R,x)+1)) #generate the total N for all R replicates
  result <- list(u_new = u_new,x=x,theta=theta,NTL=NTL)
  return(result)
}

#### LogLikehood for NPML
# NPMLEloglik <- function(param,data){
#   KK <- length(data)/2
#   N <- data[1:KK]
#   x <- data[(KK+1):length(data)]
#   UU <- param[1:KK]
#   pp <- param[(KK+1):(length(param)-1)]
#   b <- param[length(param)]
#   result <- 0
#   for(i in 1:KK){
#     result <- result + log(sum(pp,dgeom((N[i]-1),logit_inver(UU+b*x[i]))))
#   }
#   return(-result)
# }

