#------------------------------------------------------------------
# Goal: real data analysis results NPML adjusting any covariates boot
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
library(ggplot2)
#png("./data/LIFE_DATA/summary_plot.png",height=6,width=11,units="cm",
#   res=300,pointsize=10)
# ggplot(data=data.m, aes(x=nmc, y=value, fill=Current_status)) +
#   geom_bar(stat="identity")+
#   # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
#   #          color="white", size=3.5)+
#   #scale_fill_brewer(palette="Paired")+
#   theme_minimal()+
#   scale_x_continuous(breaks=seq(1,17,2))+
#   ggtitle("LIFE study summary")+
#   xlab("Number of menstrucal cycle")+
#   ylab("Number of people")
#dev.off()
ID <- unique(data$ID)
ID <- sort(ID)

obs <- rep(0,n.sub)
N <- rep(0,n.sub)
for(i in 1:n.sub){
  print(i)
  idx <- which(data$ID==ID[i])
  cbind(data[idx,]$preg,data[idx,]$method5)
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

###############start from the first enrollment cycle
###############so we take n+1 to indicate the enrollment cycle
###############since the enrollment cycle is coded as 0

#idx <- which(data.com$N!=0)
data.clean <- data.com
data.clean$N <- data.com$N+1
y <- data.clean$N
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



load("./result/start.Rdata")
uu_old <- start[[1]]
beta_old <- start[[2]]
set.seed(i1)
Rboot <- 5
NPMLEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(y)),length(y),replace = T)
  y_boot <- y[ind]
  obs_boot <- obs[ind]
  x_boot <- x[ind,]
  result <- NPMLLogFun(y_boot,x_boot,obs_boot,uu_old,beta_old)
  NPMLEst_boot[i] <- crossprod(result[[1]],result[[2]])
}

save(NPMLEst_boot,file=paste0("./result/real_data_logistic_NPML_boot",i1,".Rdata"))
# uu_new <- result[[1]]
# w_new <- result[[2]]
# beta_new <- result[[3]]
# crossprod(uu_new,w_new)
# start = list(uu_new,beta_new)
# save(start,file = )

# uu_old <- uu_new
# beta_old <- beta_new
# result_try <- NPMLLogFun(y,x,obs,uu_old,beta_old)
