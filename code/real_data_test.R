#-------------------------------------------------------------------
# Update Date: 02/06/2019
# Goal: test the best starting point for NPML
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
rm(list=ls())
args <- commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])
library(sas7bdat)
#setwd('/Users/zhangh24/GoogleDrive/project/Tom/mixture_approach_estimate_population_value/mixture_approach')
setwd('/spin1/users/zhangh24/mixture_approach')
data <- read.sas7bdat('./data/LIFE_DATA/dailycycle.sas7bdat')
data.baseline <- read.sas7bdat('./data/LIFE_DATA/baseline.sas7bdat')
library(data.table)
data.sum <- as.data.frame(fread('./data/LIFE_DATA/summary_table.csv',header=T))
data.sum <- data.sum[-1,]
data.sum[,4] <- as.numeric(data.sum[,4])
data.sum[,3] <- as.numeric(data.sum[,3])
data.sum[,2] <- as.numeric(data.sum[,2])
data.sum$still_in_study <- data.sum[,4]-data.sum[,3]-data.sum[,2]
colnames(data.sum)[6] <- "Still in study"
n.sub <- length(table(data$ID))
data.clean <- data.sum[,c(1,2,3,6)]
data.m <- melt(data.clean,id="Number of menstrual cycles")
colnames(data.m)[1] <- "nmc"
colnames(data.m)[2] <- "Current_status"

ID <- unique(data$ID)
ID <- sort(ID)

obs <- rep(0,n.sub)
N <- rep(0,n.sub)
for(i in 1:n.sub){
  print(i)
  idx <- which(data$ID==ID[i])
  obs[i] <- max(data[idx,]$preg,na.rm=T)
  N[i] <- max(data[idx,]$method5,na.rm=T)
  
}
data.temp <- cbind(ID,obs,N)

data.com <- merge(data.temp,data.baseline,by.x="ID",by.y="ID")
library(tidyverse)
library(dplyr)
data.com <- data.com %>% mutate(
  age_average = (Age_m+Age_f)/2,
  age_diff = abs(Age_m-Age_f)/2
)

###############take out the first half cycle

data.clean <- data.com
data.clean$N <- data.com$N+1
dim(data.clean)
n.couple <- nrow(data.clean)
n.cycle <- sum(data.clean$N)
Y <- rep(0,n.cycle)
age_averge.cycle <- rep(0,n.cycle)
age_diff.cycle <- rep(0,n.cycle)
ID.cycle <- rep(0,n.cycle)
temp <- 0
for(i in 1:n.couple){
  print(i)
  if(data.clean$obs[i]==1){
    Y[temp] = 1
  }
  age_averge.cycle[temp+(1:data.clean$N[i])] <- data.clean$age_average[i]
  age_diff.cycle[temp+(1:data.clean$N[i])] <- data.clean$age_diff[i]
  ID.cycle[temp+(1:data.clean$N[i])] <- data.clean$ID[i]
  temp <- temp+data.clean$N[i]
}

ObsLikfun <- function(y,x,uu_old,beta_old,w,obs){
  n <- length(y)
  result <- rep(0,n)
  xb = x%*%beta_old
  for(i in 1:n){
    result[i] <- log(sum(w*logit_inver(uu_old+xb[i])^obs[i]*(1-logit_inver(uu_old+xb[i]))^(y[i]-obs[i])))
  }
  return(sum(result))
  
}

ComLikfun <- function(y,x,uu_old,beta_old,ww,obs){
  xb = x%*%beta_old
  n <- length(y)
  temp_mat = matrix(0,n,n)
  for(i in 1:n){
    temp_mat[i,] = ((uu_old+xb[i])*obs[i]-y[i]*log(1+exp(uu_old+xb[i])))
  }
  return(sum(temp_mat*ww))
}





Mstep <- function(uu_old,beta_old,x,y,ww,alpha_x,obs){
  step <- 200
  n <- length(y)
  alpha <- 1/(n*10)
  tol <- alpha
  ComLik0 = ComLikfun(y,x,uu_old,beta_old,ww,obs)
  #Likelihood_old <- Likelihoodfun(y,x,uu_old,beta_old,ww)
  for(l in 1:step){
    
    uu_beta_old <- c(uu_old,beta_old)
    #print(uu_beta_old)
    uu_new = uu_old+alpha*gr_u_fun(uu_old,y,ww,beta_old,x,n,obs)
    beta_new <- beta_old+(alpha_x/100)*gr_b_fun(uu_new,y,ww,beta_old,x,n,obs)
    ######test the Likelihood for a few steps
    if(l%%10==0){
      ComLik_new = ComLikfun(y,x,uu_new,beta_new,ww,obs)
      if(ComLik_new>ComLik0){
        break
      }else if(ComLik_new<=ComLik0){
        uu_new = uu_old
        beta_new = beta_old
        alpha = alpha/2
        alpha_x = alpha_x/2
      }
    }
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
  }
  return(list(uu_new,beta_new,l))
}


###############logistic regression adjusting for age average and age difference
library(PAV)
model.logistic <- glm(Y~age_averge.cycle+age_diff.cycle)
summary(model.logistic)
confint(model.logistic)
tl.list <- c(0.0001,0.001,0.005,0.01)
tl <- tl.list[i1]
max_likelihood <- rep(0,length(tl))
beta_result <- matrix(0,length(tl),2)
mu_result <- matrix(0,length(tl),1)
max_step <- rep(0,length(tl))
#try different starting point
set.seed(i1)
for(s in 1:length(tl)){
  n <- nrow(data.clean)
  y=data.clean$N;
  x=cbind(data.clean$age_average,data.clean$age_diff);
  x <- as.matrix(x)
  step = 5000
  y_sm = y
  y_sm[y_sm==1] = y_sm[y_sm==1] + tl[s]
  # uu_old = seq(min(log((1/y_sm)/(1-1/y_sm))),
  #              max(log((1/y_sm)/(1-1/y_sm))),
  #              (max(log((1/y_sm)/(1-1/y_sm)))-min(log((1/y_sm)/(1-1/y_sm))))/(n-1))
  #uu_old = log((1/y_sm)/(1-1/y_sm))
  uu_old = rnorm(n,0,6)
  beta_old = rnorm(2,0,0.1)
  tol = 1e-04
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
    ww = Estep(uu_old,beta_old,x,y,w,obs)
    LikeliResult[l] <- ObsLikfun(y,x,uu_old,beta_old,w,obs)
    #rowSums(ww)
    Mstep_result = Mstep(uu_old,beta_old,x,y,ww,alpha_x,obs)
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
  max_likelihood[s] <-   LikeliResult[l]
  beta_result[s,] <- beta_new
  mu_result[s] <- crossprod(uu_new,w)
  max_step[s] = l
}

result = list(max_likelihood,beta_result,mu_result,max_step)
save(result,file = paste0("./result/npml_test",i1,".rdata"))

#plot(LikeliResult[1:l])

#test result load
setwd('/spin1/users/zhangh24/mixture_approach')
n <- 1000

likelihood_result <- rep(n,0)
beta_result <- matrix(0,n,2)
mu_result = rep(n,0)
steps = rep(n,0)
for(i1 in 1:n){
  load(paste0("./result/npml_test",i1,".rdata"))
  likelihood_result[i1] <- result[[1]]
  beta_result[i1,] <- result[[2]]
  mu_result[i1] = result[[3]]
  steps[i1] = result[[4]]

}
idx <- which.max(likelihood_result)
likelihood_result[idx]
mu_result[idx]
beta_result[idx,]
steps[idx]
