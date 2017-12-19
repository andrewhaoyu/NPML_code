library(sas7bdat)
setwd('/users/hzhang1/mixture_approach')
data <- read.sas7bdat('./data/LIFE_DATA/dailycycle.sas7bdat')
data.baseline <- read.sas7bdat('./data/LIFE_DATA/baseline.sas7bdat')

n.sub <- length(table(data$ID))
ID <- unique(data$ID)
ID <- sort(ID)

cen <- rep(0,n.sub)
N <- rep(0,n.sub)
for(i in 1:n.sub){
  print(i)
  idx <- which(data$ID==ID[i])
  cen[i] <- max(data[idx,]$preg,na.rm=T)
  N[i] <- max(data[idx,]$method5,na.rm=T)
  
}

censor.rate <- sum(cen)/length(cen)
#install_github("andrewhaoyu/PAV")
library(PAV)

NPML.estimate <- NPMLEstimateFunction(N,cen)
UU <- NPML.estimate[[1]]
pp <- NPML.estimate[[2]]
preg.probability.est <- NPML.estimate[[3]]

UU.cut <- seq(0,1,0.0001)
pp.cut <- rep(0,length(UU.cut))
for(i in 1:length(UU.cut)){
  idx <- which(UU<=UU.cut[i])
  pp.cut[i] <- sum(pp[idx])
}



png("./result/cdf_NPML.png")
plot(UU.cut,pp.cut,type="l",xlab="probability mass estimate for pregnancy probability",
     ylab="cumulative probability",main="cumulative curve of NPML estimate pregnancy probability")
dev.off()

data.clean <- data.frame(ID,N,cen)
colnames(data.clean) <- c("ID","men_cycles","pregnant")
data.clean %>% filter(men_cycles==0&pregnant==0)
ggplot(data=data,aes(days))+geom_histogram(colour = "darkgreen", fill = "white", binwidth = 25)+scale_x_continuous(expand = c(0,0))

source('/code/plot_theme.R')


library(ggplot2); library(scales); library(grid); library(RColorBrewer)

png("/Users/haoyuzhang/Dropbox/project/cristian/project1/NPML_paper/Oct_28/manuscriptv6/distribution_of_patients_all.png",width = 8,height = 6,units = "in",res = 600)
ggplot(data, aes(days)) +
  geom_histogram(binwidth=20, fill="#c0392b", alpha=0.75) +
  fte_theme() +
  #theme_Publication()
  labs(x="Time to Heal(days)", y="# of patients") +
  scale_x_continuous(breaks=seq(0,1000, by=50)) +
  scale_y_continuous(labels=comma,breaks = seq(0,30,by=5))+
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


