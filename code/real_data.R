library(sas7bdat)
setwd('/Users/zhangh24/GoogleDrive/project/Tom/mixture_approach_estimate_population_value/mixture_approach')
#setwd('/spin1/users/zhangh24/mixture_approach')
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
png("./data/LIFE_DATA/summary_plot.png",height=6,width=11,units="cm",
    res=300,pointsize=10)
ggplot(data=data.m, aes(x=nmc, y=value, fill=Current_status)) +
  geom_bar(stat="identity")+
 # geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
  #          color="white", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+
  scale_x_continuous(breaks=seq(1,17,2))+
  ggtitle("LIFE study summary")+
  xlab("Number of menstrucal cycle")+
  ylab("Number of people")
dev.off()
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




###############logistic regression adjusting for age average and age difference
model.logistic <- glm(Y~age_averge.cycle+age_diff.cycle)
summary(model.logistic)
confint(model.logistic)
##############mixed effect logistic regression allowing for sample difference
library(lme4)
model.mix.logistic <- glmer(Y~(1|ID.cycle)+age_averge.cycle+age_diff.cycle,family = binomial)
summary(model.mix.logistic)
confint(model.mix.logistic)
#############NPML model
# uu0 = seq(summary(model.logistic)$coefficients[1,1]-5,
#           summary(model.logistic)$coefficients[1,1]+5,
#           10/(n-1))
#   
# beta0 = summary(model.logistic)$coefficients[2:3,1]
#model.NPMLlog <- NPMLLogFun(y=data.clean$N,x=cbind(data.clean$age_average,data.clean$age_diff),uu0,beta0)
tl <- c(0.005)
max_likelihood <- rep(0,length(tl))
beta_result <- matrix(0,length(tl),2)
mu_result <- matrix(0,length(tl),1)
max_step <- rep(0,length(tl))
#try different starting point
for(s in 1:length(tl)){
  n <- nrow(data.clean)
  y=data.clean$N;
  x=cbind(data.clean$age_average,data.clean$age_diff);
  x <- as.matrix(x)
  step = 2000
  y_sm = y
  y_sm[y_sm==1] = y_sm[y_sm==1] + tl[s]
   uu_old = seq(min(log((1/y_sm)/(1-1/y_sm))),
                max(log((1/y_sm)/(1-1/y_sm))),
                (max(log((1/y_sm)/(1-1/y_sm)))-min(log((1/y_sm)/(1-1/y_sm))))/(n-1))
  #uu_old = log((1/y_sm)/(1-1/y_sm))
  beta_old = summary(model.logistic)$coefficients[2:3,1]
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
  plot(LikeliResult[1:l])
}
  
  


Mstep <- function(uu,beta,x,y,ww,alpha_x){
  uu_old <- uu
  beta_old <- beta
  step <- 100
  n <- length(y)
  alpha <- 1/n
  tol <- 0.001
  for(l in 1:step){
    
    uu_beta_old <- c(uu_old,beta_old)
    #print(uu_beta_old)
    uu_new = uu_old+alpha*gr_u_fun(uu_old,y,ww,beta_old,x,n)
    beta_new <- beta_old+(alpha_x/10)*gr_b_fun(uu_new,y,ww,beta_old,x,n)
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
  }
  return(list(uu_new,beta_new))
}



##############original estimate without any confounder adjustment


for(i in 1:n.sub){
  print(i)
  idx <- which(data$ID==ID[i])
  obs[i] <- max(data[idx,]$preg,na.rm=T)
  N[i] <- max(data[idx,]$method5,na.rm=T)
  
}
###########take out people with no complete menstrcucal cycle
table(obs,N)
idx.new <- which(N!=0)
obs.new <- obs[idx.new]
N.new <- N[idx.new]


N <- N.new
cen <- obs.new
#################new data to analyze (all people with complete menstrual cycle)
#################N is the number of menstrual cycle
#################cen 1 represents pregnant, 0 represents dropped out of study
new.data <- as.data.frame(table(N,cen))

risk <- rep(0,max(N))

for(i in 1:max(N)){
  if(i==1){
    idx <- which(new.data[,1]==i)
    remove <- sum(new.data[idx,3])
    risk[i] <- length(N)
  }else{
    idx <- which(new.data[,1]==i)
    remove <- sum(new.data[idx,3])
    risk[i] <- risk[i-1]-remove  
  }
}


sur.data <- cbind(table(N,cen),risk)

K.estimate <- sur.data[,2]/risk


censor.rate <- sum(cen)/length(cen)
library(devtools)
install_github("andrewhaoyu/PAV")
library(PAV)


M <- 1
Nt <- N

################Stratified estimate with no covariates adjustment
StratEst <- StratEstimateFunction(Nt,cen)
################Pooled estimate with no covariates adjustment
PooledEst <- PooledEstimateFunction(Nt,cen)
################NPML estimate with no covariates adjustment
NPML.estimate <- NPMLEstimateFunction(N,cen)
NPMLEst <- NPML.estimate[[3]]
N <- Nt
###############beta prior estimate
begin <- c(0.3804882,20.0296849)


fit <-  optim(par = begin,fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.9995,10000),method="L-BFGS-B",control=list(fnscale=-1))
BetaEst <- fit$par[1]
like <- BetaGeometricLikehood(fit$par)












#set.seed(123)
#Rboot <- 300

# n <- 5000
# begin_point1 <- runif(n,0.005,0.9)
# begin_point2 <- runif(n,1,500)
# begin_point <- cbind(begin_point1,begin_point2)
# try_result <- rep(0,n)
# for(i in 1:n){
#   fit <-  optim(par = begin_point[i,],fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.9,10000),method="L-BFGS-B",control=list(fnscale=-1))
#   try_result[i] <- BetaGeometricLikehood(fit$par)
#   
# }
# idx <- which.max(try_result)
bootdata <- data.frame(N=N,cen=cen)
# NPMLEst_boot <- boot(data.frame(N=N,cen=cen),NPMLEstimateFunction_mean,Rboot)
# PooledEst_boot <- boot(N,PooledEstimateFunction,Rboot)
# StratEst_boot <- boot(N,StratEstimateFunction,Rboot)
set.seed(123)
Rboot <- 10000
# NPMLEst_boot <- rep(0,Rboot)
# 
# for(i in 1:Rboot){
#   print(i)
#   ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
#   N <- Nt[ind]
#   cend <- cen[ind]
#   NPMLEst_boot[i] <- NPMLEstimateFunction_mean(N,c(1:K),cend)
# }

PooledEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  PooledEst_boot[i] <- PooledEstimateFunction(N,cend)
}

StratEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  StratEst_boot[i] <- StratEstimateFunction(N,cend)
}
#Rboot <- 300
BetaEst_boot <- rep(0,Rboot)

cen_original <- cen
for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cen <- cen_original[ind]

  
  BetaEst_boot[i] <- ParaBetaEstimateFunction(N)
  
}
cen <- cen_original
#BetaEst_boot_f <- BetaEst_boot[BetaEst_boot<0.1]

places <- 3
NPMLEstRound <- round(NPMLEst,places)
StratEstRound <- round(StratEst,places)
PooledEstRound <- round(PooledEst,places)
BetaEstRound <- round(BetaEst,places)
NPMLEstLow <- round(quantile(NPMLEst_boot,0.025),places)
NPMLEstHigh <- round(quantile(NPMLEst_boot,0.975),places)
StratEstLow <- round(quantile(StratEst_boot,0.025),places)
StratEstHigh <- round(quantile(StratEst_boot,0.975),places)
PooledEstLow <- round(quantile(PooledEst_boot,0.025),places)
PooledEstHigh <- round(quantile(PooledEst_boot,0.975),places)
BetaEstLow <- round(quantile(BetaEst_boot,0.025),places)
BetaEstHigh <- round(quantile(BetaEst_boot,0.975),places)

NPMLresult <- paste0(NPMLEstRound,"(",NPMLEstLow,"-",NPMLEstHigh,")")
Pooledresult <- paste0(PooledEstRound,"(",PooledEstLow,"-",PooledEstHigh,")")
Stratresult <-  paste0(StratEstRound,"(",StratEstLow,"-",StratEstHigh,")")
Betaresult <- paste0(BetaEstRound,"(",BetaEstLow,"-",BetaEstHigh,")")

realdataresult <- data.frame(Strat=Stratresult,Pooled=Pooledresult,
                             NPML=NPMLresult,Beta=Betaresult)
realdataresult <- NULL
save(realdataresult,file="./result/realdataresult.Rdata")






# ggplot(data=data,aes(days))+geom_histogram(colour = "darkgreen", fill = "white", binwidth = 25)+scale_x_continuous(expand = c(0,0))
# 
# source('/code/plot_theme.R')
# 
# 
# library(ggplot2); library(scales); library(grid); library(RColorBrewer)
# 
# png("/Users/haoyuzhang/Dropbox/project/cristian/project1/NPML_paper/Oct_28/manuscriptv6/distribution_of_patients_all.png",width = 8,height = 6,units = "in",res = 600)
# ggplot(data, aes(days)) +
#   geom_histogram(binwidth=20, fill="#c0392b", alpha=0.75) +
#   fte_theme() +
#   #theme_Publication()
#   labs(x="Time to Heal(days)", y="# of patients") +
#   scale_x_continuous(breaks=seq(0,1000, by=50)) +
#   scale_y_continuous(labels=comma,breaks = seq(0,30,by=5))+
#   geom_hline(yintercept=0, size=0.4, color="black")
# 
# 
# 
# dev.off()
# 
# 
# 
# mu.mean <- rep(c(rep(0.1,3),rep(0.3,3),rep(0.5,3)),12)
# M <- rep(c(rep(0.5,9),rep(5,9),rep(9999,9)),4)
# K <- rep(c(20,75,125),36)
# MSE <- c(53,37,38,92,74,70,104,81,73,148,133,140,398,356,330,368,320,300,
#          287,248,251,552,497,475,442,393,378,54,79,83,714,768,775,1990,2068,2080,42,52,54,142,190,210,122,123,127,5,1,1,40,10,5,64,15,10,14,5,3,29,10,5,47,11,7,22,7,5,62,20,10,83,24,12,33,7,6,92,22,11,95,24,14,19,7,4,89,58,25,50,33,8,14,6,3,32,11,7,74,15,10,5,1,1,40,10,5,64,15,10)
# 
# method <- c(rep("Stratified",27),rep("Pooled",27),rep("NPML",27),rep("Beta Prior",27))
# #M[M==10000] <- 8
# 
# 
# geometric_result <- data.frame(MSE=MSE,method=method,mu.mean=mu.mean,M=M,K=K)
# geometric_result <- geometric_result%>%mutate(group=paste0(mu.mean," M = ",M))
# idx.M.5 <- which(M==5)
# idx.M.9999 <- which(M==9999)
# temp.M.5 <- geometric_result[idx.M.5,]
# geometric_result[idx.M.5,] <- geometric_result[idx.M.9999,]
# geometric_result[idx.M.9999,] <- temp.M.5
# idx1 <- which(method=="Pooled"&mu.mean==0.5&M==0.5)
# geometric_result[idx,] <- NA
# idx2 <- which(method=="Pooled"&M==0.5&mu.mean==0.3)
# geometric_result[idx2,] <- NA
# idx3 <- which(method=="Pooled"&M==0.5&mu.mean==0.5)
# geometric_result[idx3,] <- NA
# idx4 <- which(method=="Stratified"&M==9999&mu.mean==0.1)
# geometric_result[idx4,] <- NA
# idx5 <- which(method=="Stratified"&M==5&mu.mean==0.1)
# geometric_result[idx5,] <- NA
# idx6 <- which(method=="Stratified"&M==5&mu.mean==0.3)
# geometric_result[idx6,] <- NA
# idx7 <- which(method=="Stratified"&M==9999&mu.mean==0.5)
# geometric_result[idx7,] <- NA
# idx8 <- which(method=="Stratified"&M==5&mu.mean==0.5)
# geometric_result[idx8,] <- NA
# na.id <- c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8)
# geometric_result$mu.mean <- paste0("expression(mu==",mu.mean,")")
# geometric_result$mu.mean[na.id] <- NA
# 
# #  geometric_result$group.f <- factor(geometric_result$group)
# #  geometric_result$group.f
# # levels(geometric_result$group.f) <- levels(geometric_result$group.f)[c(1,3,2,4,6,5,7,9,8,10)]
# # geometric_result$group.f
# # geometric_result$group <- as.character(geometric_result$group)
# 
# ##########to make M ordered
# # temp <- geometric_result[10,]
# # geometric_result[10,] <- geometric_result[37,]
# # geometric_result[37,] <- temp
# palette <- brewer.pal("Greys", n=9)
# color.background = palette[1]
# color.grid.major = palette[3]
# color.axis.text = palette[6]
# color.axis.title = palette[7]
# color.title = palette[9]
# group <- geometric_result$group
# group <- geometric_result$group
# colnames(geometric_result)[2] <- "Method"
# png("/Users/haoyuzhang/Dropbox/project/cristian/project1/NPML_paper/Oct_28/manuscriptv6/MSE_uncengeometric.png",width = 8,height = 6,units = "in",res = 600)
# ggplot(data=na.omit(geometric_result),aes(x=K,y=MSE,group=Method,color=Method))+
#   geom_line(size=0.8)+
#   geom_point(size=1.1)+
#   #geom_line()+
#   #scale_color_ptol("method") +
#   facet_wrap(~group,scales = "free",labeller = label_bquote(mu==.(group)))+
#   #facet_wrap(mu.mean~M,scales = "free",labeller=label_bquote(mu==.(mu.mean)))+
#   #facet_grid(mu.mean~M,scales = "free",labeller = label_bquote(mu==.(mu.mean),M==.(M)))+
#   #theme_Publication()+
#   theme_minimal()+
#   theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
#   theme(panel.grid.minor=element_blank()) +
#   theme(axis.ticks=element_blank()) +
#   theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
#   theme(axis.text.x=element_text(size=10,color=color.axis.text)) +
#   theme(axis.text.y=element_text(size=10,color=color.axis.text)) +
#   theme(axis.title.x=element_text(size=14,color=color.axis.title, face="bold",vjust=0)) +
#   theme(axis.title.y=element_text(size=14,color=color.axis.title, face="bold",vjust=1.25))+
#   theme(strip.text = element_text(size = 11,face="bold"))+
#   scale_x_discrete(limits=c(25,75,125))+
#   labs(x="K(Sample Size)", y="MSEx1e4")+
#   scale_color_manual(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c"))
# 
# dev.off()
# 
# 
# 
# ggplot(data=na.omit(geometric_result),aes(x=K,y=MSE,group=Method,color=Method))+
#   geom_line(size=0.8)+
#   geom_point(size=1.1)+
#   facet_wrap(~group,scales = "free",labeller = label_bquote(mu==.(group)))+
#   scale_colour_Publication()+ theme_Publication()
# 
# #theme_fivethirtyeight()
# #theme_bw()
# 
# my.label_bquote <- function (expr1 = (mu == .(x)),expr2 = (M == .(x))) 
# {
#   quoted1<- substitute(expr1)
#   quoted2 <- substitute(expr2)
#   function(variable, value) {
#     value <- as.character(value)
#     if(variable == 'mu.mean')
#       lapply(value, function(x)
#         eval(substitute(bquote(expr1, list(x = x)),list(expr1 = quoted1))))
#     else
#       lapply(value, function(x) 
#         eval(substitute(bquote(expr2, list(x = x)),list(expr2 = quoted2))))
#   }
# }
# 
# 
# 
# 
# theme_Publication <- function(base_size=14, base_family="Times") {
#   library(grid)
#   library(ggthemes)
#   (theme_foundation(base_size=base_size, base_family=base_family)
#     + theme(plot.title = element_text(face = "bold",
#                                       size = rel(1.2), hjust = 0.5),
#             text = element_text(),
#             panel.background = element_rect(colour = NA),
#             plot.background = element_rect(colour = NA),
#             panel.border = element_rect(colour = NA),
#             axis.title = element_text(face = "bold",size = rel(1)),
#             axis.title.y = element_text(angle=90,vjust =2),
#             axis.title.x = element_text(vjust = -0.2),
#             axis.text = element_text(), 
#             axis.line = element_line(colour="black"),
#             axis.ticks = element_line(),
#             panel.grid.major = element_line(colour="#f0f0f0"),
#             panel.grid.minor = element_blank(),
#             legend.key = element_rect(colour = NA),
#             legend.position = "bottom",
#             legend.direction = "horizontal",
#             legend.key.size= unit(0.2, "cm"),
#             # legend.margin = unit(0, "cm"),
#             legend.title = element_text(face="italic"),
#             plot.margin=unit(c(10,5,5,5),"mm"),
#             strip.background=element_rect(colour="white",fill="white"),
#             strip.text = element_text(face="bold")
#     ))
#   
# }
# 
# scale_fill_Publication <- function(...){
#   library(scales)
#   discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
#   
# }
# 
# scale_colour_Publication <- function(...){
#   library(scales)
#   discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
#   
# }
# 
# load("variance.result.Rdata")
# library(ggplot2)
# data <- variance.result[[3]]
# data[data[,4]==10^4,4] = 9
# data <- as.data.frame(data)
# data <- data%>%mutate(group=paste0(Theta," M = ",M))
# png("/Users/haoyuzhang/Dropbox/project/cristian/project1/NPML_paper/Oct_28/manuscriptv6/variance_uncengeometric.png",width = 9.7,height = 6,units = "in",res = 600)
# ggplot(data,aes(factor(K),Npml.Estimate.Variance))+geom_boxplot(outlier.shape = 16, outlier.size = 1)+
#   facet_wrap(~group,scales = "free",labeller = label_bquote(mu==.(group)))+
#   geom_point(aes(factor(K),True.Var),size=1.2,shape=15,colour="red")+
#   scale_colour_Publication()+theme_Publication()+
#   scale_y_continuous(limits=c(0,0.25))+
#   labs(x="K(Sample Size)", y="NPML Variance Estimate ")
# dev.off()
# 
# 
# png("/Users/haoyuzhang/Dropbox/project/cristian/project1/NPML_paper/Oct_28/manuscriptv6/p15_uncengeometric.png",width = 9.7,height = 6,units = "in",res = 600)
# ggplot(data,aes(factor(K),Npml.Estimate.below10))+geom_boxplot(outlier.shape = 16, outlier.size = 1)+
#   facet_wrap(~group,scales = "free",labeller = label_bquote(mu==.(group)))+
#   geom_point(aes(factor(K),True.below10),size=1.2,shape=15,colour="red")+
#   scale_colour_Publication()+theme_Publication()+
#   labs(x="K(Sample Size)", y=expression(paste("NPML Estimate Pr(",theta,"<0.2)"))) 
# dev.off()
# 
# 
# 


