library(sas7bdat)
setwd('/Users/zhangh24/GoogleDrive/project/Tom/mixture_approach_estimate_population_value/mixture_approach')
data <- read.sas7bdat('./data/LIFE_DATA/dailycycle.sas7bdat')
data.baseline <- read.sas7bdat('./data/LIFE_DATA/baseline.sas7bdat')
##############clean the data and get everyone's pregnancy, censored, and menstrucal cycle
n.sub <- length(table(data$ID))
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

table(obs,N)
total.cycles <- 0 
for(i in 1:N)



idx.new <- which(N!=0)
obs.new <- obs[idx.new]
N.new <- N[idx.new]


N <- N.new
cen <- obs.new
censor.rate <- sum(cen)/length(cen)
library(devtools)
install_github("andrewhaoyu/PAV")
library(PAV)

NPML.estimate <- NPMLEstimateFunction(N,cen)
UU <- NPML.estimate[[1]]
pp <- NPML.estimate[[2]]
mean.estimate <- sum(UU%*%pp)
var.estimate <- sum((UU-mean.estimate)^2%*%pp)

M <- mean.estimate*(1-mean.estimate)/var.estimate-1











mean.estimate <- 0.175
M <- 19.58
alpha <- mean.estimate*M
beta <- (1-mean.estimate)*M

x <- seq(0,1,0.0001)

y <- dbeta(x,alpha,beta)

plot(x,y,type='l',main="beta distribution density estimate")
























preg.probability.est <- NPML.estimate[[3]]

UU.cut <- seq(0,1,0.0001)
pp.cut <- rep(0,length(UU.cut))
for(i in 1:length(UU.cut)){
  idx <- which(UU<=UU.cut[i])
  pp.cut[i] <- sum(pp[idx])
}



png("./result/cdf_NPML.png",width = 8,height = 6,units = "in",res = 600)
plot(UU.cut,pp.cut,type="l",xlab="probability mass estimate for pregnancy probability",
     ylab="cumulative probability",main="cumulative curve of NPML estimate pregnancy probability")
dev.off()

ID <- ID[idx.new]
data.clean <- data.frame(ID,N,cen)
colnames(data.clean) <- c("ID","men_cycles","pregnant")
library(tidyverse)
#data.clean %>% filter(men_cycles==0&pregnant==0)
#ggplot(data=data.clean,aes(men_cycles))+geom_histogram(colour = "darkgreen", fill = "white")+scale_x_continuous(expand = c(0,0))

source('./code/plot_theme.R')


library(ggplot2); library(scales); library(grid); library(RColorBrewer)

png("./result/men_cyc_dis.png",width = 8,height = 6,units = "in",res = 600)
ggplot(data.clean, aes(men_cycles)) +
  geom_histogram( fill="#c0392b", alpha=0.75) +
  fte_theme() +
  #theme_Publication()
  labs(x="Time to Heal(days)", y="# of patients") +
  #scale_x_continuous(breaks=seq(0,1000, by=50)) +
  #scale_y_continuous(labels=comma)+
  geom_hline(yintercept=0, size=0.4, color="black")



dev.off()









source("./code/geometric_fun.R")
set.seed(123)
library(boot)
library(ggplot2)
#library(ggthemes)
library(dplyr)
Rboot <- 300
M <- 1
K <- n.sub
#Nt <- data[data[,2]=='H',1]
#data.uncen <- data.frame(days=Nt)
Nt <- N
cen <- cen
NPMLEst <- NPMLEstimateFunction_mean(Nt,c(1:K),cen)


