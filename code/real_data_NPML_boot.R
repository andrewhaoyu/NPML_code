args <- commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])
print(i1)
set.seed(i1)

library(sas7bdat)
setwd('/users/hzhang1/mixture_approach')
data <- read.sas7bdat('./data/LIFE_DATA/dailycycle.sas7bdat')
data.baseline <- read.sas7bdat('./data/LIFE_DATA/baseline.sas7bdat')
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
idx.new <- which(N!=0)
obs.new <- obs[idx.new]
N.new <- N[idx.new]


N <- N.new
cen <- obs.new
censor.rate <- sum(cen)/length(cen)
library(devtools)
#install_github("andrewhaoyu/PAV")
library(PAV)

#NPML.estimate <- NPMLEstimateFunction(N,cen)


#bootdata <- data.frame(N=N,cen=cen)
Nt <- N
Rboot <- 10
NPMLEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  NPMLEst_boot[i] <- NPMLEstimateFunction(N,cend)[[3]]
}



save(NPMLEst_boot,file=paste0("./result/realdataresult.",i1,"Rdata"))






