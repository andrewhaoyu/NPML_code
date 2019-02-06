n = 500
x1 = rnorm(n,1,4)
x2 = rnorm(n,31,4.2)

u = rnorm(n,0,2)
beta1 = 0.01
beta2 = 0.02
p = logit_inver(u+beta1*x1+beta2*x2)
y = rgeom(n,p)+1
y_sm = y+1
uu0 = seq(min(log((1/y_sm)/(1-1/y_sm))),
          max(log((1/y_sm)/(1-1/y_sm))),
          (max(log((1/y_sm)/(1-1/y_sm)))-min(log((1/y_sm)/(1-1/y_sm))))/(n-1))
beta0 = c(0,0)
x = cbind(x1,x2)

#model.NPMLlog <- NPMLLogFun(y=y,x=cbind(x1,x2),uu0,beta0)
ObsLikfun <- function(y,x,uu_old,beta_old,w){
  n <- length(y)
  result <- rep(0,n)
  xb = x%*%beta_old
  for(i in 1:n){
    result[i] <- log(sum(w*logit_inver(uu_old+xb[i])*(1-logit_inver(uu_old+xb[i]))^(y[i]-1)))
  }
  return(sum(result))
  
}

ComLikfun <- function(y,x,uu_old,beta_old,ww){
  xb = x%*%beta_old
  n <- length(y)
  temp_mat = matrix(0,n,n)
  for(i in 1:n){
    temp_mat[i,] = (uu_old+xb[i]-y[i]*log(1+exp(uu_old+xb[i])))
  }
  return(sum(temp_mat*ww))
}










y=y;
x=cbind(x1,x2);
x <- as.matrix(x)
step = 1500
uu_old = uu0
beta_old = beta0
tol = 0.0001
n = length(y)
w = rep(1/n,n)
#set the step length of gradient decent to aviod unconvergence
var.x = apply(x,2,var)
alpha_x = rep(1/n,length(beta_old))
for(i in 1:length(beta_old)){
  if(var.x[i]>=1){
    alpha_x[i] = alpha_x[i]/var.x[i]
  }
}
LikeliResult <- rep(0,step)
StepResult <- rep(0,step)
#library(PAV)
for(l in 1:step){
  uu_beta_old <- c(uu_old,beta_old)
  print(uu_beta_old)
  #print(uu_beta_old)
  ww = Estep(uu_old,beta_old,x,y,w)
  LikeliResult[l] <- ObsLikfun(y,x,uu_old,beta_old,w)
  #rowSums(ww)
  Mstep_result = Mstep(uu_old,beta_old,x,y,ww,alpha_x)
  #Mstep_result = Mstep2(uu_old,beta_old,x,y,ww)
  uu_new = Mstep_result[[1]]
  beta_new = Mstep_result[[2]]
  #StepResult[l] = Mstep_result[[3]]
  uu_beta_new <- c(uu_new,beta_new)
  error <- max(abs(uu_beta_new-uu_beta_old))
  if(error<tol){
    break
  }
  uu_old <- uu_new
  beta_old <- beta_new
  w = colSums(ww)/sum(ww)
}
plot(LikeliResult[1:l])
min(diff(LikeliResult[1:l]))>=0

Mstep <- function(uu_old,beta_old,x,y,ww,alpha_x){
  step <- 200
  n <- length(y)
  alpha <- 1/(n*10)
  tol <- alpha/10
  ComLik0 = ComLikfun(y,x,uu_old,beta_old,ww)
  #Likelihood_old <- Likelihoodfun(y,x,uu_old,beta_old,ww)
  for(l in 1:step){
    
    uu_beta_old <- c(uu_old,beta_old)
    #print(uu_beta_old)
    uu_new = uu_old+alpha*gr_u_fun(uu_old,y,ww,beta_old,x,n)
    beta_new <- beta_old+(alpha_x/100)*gr_b_fun(uu_new,y,ww,beta_old,x,n)
    ComLik_new = ComLikfun(y,x,uu_new,beta_new,ww)
    if(ComLik_new>ComLik0){
      break
    }else if(ComLik_new<=ComLik0){
      uu_new = uu_old
      beta_new = beta_old
      alpha = alpha/2
      alpha_x = alpha_x/2
    }
    # uu_beta_new <- c(uu_new,beta_new)
    # error <- max(abs(uu_beta_new-uu_beta_old))
    # if(error<tol){
    #   break
    # }
    # uu_old <- uu_new
    # beta_old <- beta_new
  }
  return(list(uu_new,beta_new,l))
}


#Second derivative function
SeDFunction <- function(uu_old,beta_old,x,y,ww){
  n <- length(uu_old)
  q <- length(beta_old)
  ##the second derivatives is a n+q n+q matrix
  result <- matrix(0,n+q,n+q)
  temp_mat <- matrix(0,n,n)
  ##first calculate the second derivative for uu_old
  for(i in 1:n){
    p = logit_inver(u+xb[i])
    temp_mat[i,] <- y[i]*p*(1-p)
  }
  temp_mat_ww = temp_mat*ww
  diag(result)[1:n] <- colSums(temp_mat_ww)
  ##the second derivative for uu_old beta_old
  for(k in 1:q){
    x_temp <- x[,k]
    #use the element times in R
    result[1:n,n+k] <- x_temp%*%(temp_mat_ww)
    result[n+k,1:n] <- result[1:n,n+k]
  }
  ##the second derivative for beta_old beta_old
  for(k in 1:q){
    x_temp <- x[,k]
    #use the element times in R
    result[n+k,n+k] <- x_temp%*%(temp_mat_ww)%*%x_temp
    
  }
 return(result) 
  
}




Mstep2 <- function(uu_old,beta_old,x,y,ww){
  step <- 200
  n <- length(y)
  
  tol = 0.001
  xb = x%*%beta_old
  
  #Likelihood_old <- Likelihoodfun(y,x,uu_old,beta_old,ww)
  for(l in 1:step){
    
    uu_beta_old <- c(uu_old,beta_old)
    print(uu_beta_old)
    M <- SeDFunction(uu_old,beta_old,x,y,ww)
    score <- c(gr_u_fun(uu_old,y,ww,beta_old,x,n),
               gr_b_fun(uu_old,y,ww,beta_old,x,n))
    uu_beta_new <- uu_beta_old+solve(M)%*%score
    uu_new <- uu_beta_new[1:length(y)]
    beta_new <- uu_beta_new[(length(y)+1):length(uu_beta_new)]
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
  }
  return(list(uu_new,beta_new,l))
}












