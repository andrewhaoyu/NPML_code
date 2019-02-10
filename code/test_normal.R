#-------------------------------------------------------------------
# Update Date: 02/03/2019
# Goal: test gradient descent and EM algorithm on mixture normal example
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
n <- 1000
mu1 = 1
mu2 = 4
p = 0.3
p1 = rbinom(n,1,p)
y = p1*rnorm(n,mu1,1)+(1-p1)*rnorm(n,mu2,1)
hist(y)


######gradient decent
mu1_start = 0
mu2_start = 5
p0 = 0.5
gr_p <- function(y,mu1_old,mu2_old,p_old){
  temp = (dnorm(y,mu1_old)-dnorm(y,mu2_old))/
        (p_old*dnorm(y,mu1_old)+(1-p_old)*dnorm(y,mu2_old))
  return(sum(temp))
}
gr_mu1 = function(y,mu1_old,mu2_old,p_old){
  temp = (p_old*dnorm(y,mu1_old)*(y-mu1_old))/
    (p_old*dnorm(y,mu1_old)+(1-p_old)*dnorm(y,mu2_old))
  return(sum(temp))
} 
gr_mu2 = function(y,mu1_old,mu2_old,p_old){
  temp = ((1-p_old)*dnorm(y,mu2_old)*(y-mu2_old))/
    (p_old*dnorm(y,mu1_old)+(1-p_old)*dnorm(y,mu2_old))
  return(sum(temp))
} 

mu1_old = mu1_start
mu2_old = mu2_start
p_old = p0
step = 300
alpha = 1/(10*n)
tol = 0.001
for(l in 1:step){
mu_p_vec_old <- c(mu1_old,mu2_old,p_old)
print(mu_p_vec_old)
p_new = p_old + alpha*gr_p(y,mu1_old,mu2_old,p_old)
mu1_new = mu1_old + alpha*gr_mu1(y,mu1_old,mu2_old,p_new)
mu2_new = mu2_old + alpha*gr_mu2(y,mu1_new,mu2_old,p_new)
mu_p_vec_new <- c(mu1_new,mu2_new,p_new)
error = max(abs(mu_p_vec_new-mu_p_vec_old))
if(error<tol){
  break
}
mu1_old = mu1_new
mu2_old = mu2_new
p_old = p_new
}


#######EM algorithm
#######E step
w = rep(0,n)
PostFun <- function(y,mu1_old,mu2_old,p_old){
  result <- dnorm(y,mu1_old)*p_old/
    (dnorm(y,mu1_old)*p_old+dnorm(y,mu2_old)*(1-p_old))
  return(result)
}
#######M step
MstepFun <- function(y,w){
  mu1_new = sum(y*w)/sum(w)
  mu2_new = sum(y*(1-w))/sum(1-w)
  p_new = sum(w)/length(y)
  return(list(mu1_new,mu2_new,p_new))
}
#######M step
ObsLihfun <- function(y,mu1_old,mu2_old,p_old){
  return(sum(log(dnorm(y,mu1_old)*p_old+dnorm(y,mu2_old)*(1-p_old))))
}
#######EM step
mu1_start = 0
mu2_start = 5
p0 = 0.5
mu1_old = mu1_start
mu2_old = mu2_start
p_old = p0
step = 300
tol = 0.001
likelihood_result <- rep(0,step)
for(l in 1:step){
  mu_p_vec_old <- c(mu1_old,mu2_old,p_old)
  w <- PostFun(y,mu1_old,mu2_old,p_old)
  likelihood_result[l] <- ObsLihfun(y,mu1_old,mu2_old,p_old)
  print(mu_p_vec_old)
  Mstep_result <- MstepFun(y,w)
  
  p_new = Mstep_result[[3]]
  mu1_new = Mstep_result[[1]]
  mu2_new = Mstep_result[[2]]
  mu_p_vec_new <- c(mu1_new,mu2_new,p_new)
  error = max(abs(mu_p_vec_new-mu_p_vec_old))
  if(error<tol){
    break
  }
  mu1_old = mu1_new
  mu2_old = mu2_new
  p_old = p_new
}
plot(likelihood_result[1:l])
min(diff(likelihood_result[1:l]))>0


####test the MM algorithm
####beta_new = beta_old +4*(t(x)x)^-1%*%t(x)(Y-P)
n <- 1000
x1 = rnorm(n,0,1)
x2 = rnorm(n,0,2)
p = logit_inver(1+2*x1+4*x2)
beta_old = c(0,0,0)
y = rbinom(n,1,p)
#glm(y~cbind(x1,x2),family=binomial(logit))
x <- cbind(1,x1,x2)
time1 = proc.time()
for(h in 1:1000){
  n <- 1000
  x1 = rnorm(n,0,1)
  x2 = rnorm(n,0,2)
  p = logit_inver(1+2*x1+4*x2)
  beta_old = c(0,0,0)
  y = rbinom(n,1,p)
  #glm(y~cbind(x1,x2),family=binomial(logit))
  x <- cbind(1,x1,x2)
  xx.inver.xt <- solve(t(x)%*%x)%*%t(x)
  tol = 0.001
  for(l in 1:1000){
    #print(beta_old)
    p_old <- logit_inver(x%*%beta_old)
    beta_new = beta_old + 4*xx.inver.xt%*%(y-p_old)
    error = max(abs(beta_new-beta_old))
    if(error<tol){
      break
    }
    beta_old = beta_new
  }
  
}
time1 = proc.time()-time1 

#gradient descent
time2 = proc.time()
for(h in 1:1000){
  n <- 1000
  x1 = rnorm(n,0,1)
  x2 = rnorm(n,0,2)
  p = logit_inver(1+2*x1+4*x2)
  beta_old = c(0,0,0)
  y = rbinom(n,1,p)
  #glm(y~cbind(x1,x2),family=binomial(logit))
  x <- cbind(1,x1,x2)
  alpha = 1/n
  tol = 0.001
  for(l in 1:1000){
    #print(beta_old)
    p_old <- logit_inver(x%*%beta_old)
    beta_new = beta_old + alpha*t(x)%*%(y-p_old)
    error = max(abs(beta_new-beta_old))
    if(error<tol){
      break
    }
    beta_old = beta_new
  }
  
}
time2 = proc.time()-time2

#glm 
time3 = proc.time()
for(h in 1:1000){
  n <- 1000
  x1 = rnorm(n,0,1)
  x2 = rnorm(n,0,2)
  p = logit_inver(1+2*x1+4*x2)
  beta_old = c(0,0,0)
  y = rbinom(n,1,p)
  x <- cbind(1,x1,x2)
  w = matrix(0,n,n)
  tol = 0.001
  for(l in 1:1000){
   # print(beta_old)
    p_old <- logit_inver(x%*%beta_old)
    diag(w) = p_old*(1-p_old)
    beta_new = beta_old + solve(t(x)%*%w%*%x)%*%t(x)%*%(y-p_old)
    #beta_new = beta_old + 4*xx.inver.xt%*%(y-p_old)
    error = max(abs(beta_new-beta_old))
    if(error<tol){
      break
    }
    beta_old = beta_new
  }
  
 
 
  
}
time3 = proc.time()-time3

