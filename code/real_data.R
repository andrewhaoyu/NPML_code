
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
StratEst <- StratEstimateFunction(Nt,c(1:K),cen)
PooledEst <- PooledEstimateFunction(Nt,c(1:K),cen)
NPMLEst <- NPMLEstimateFunction_mean(Nt,c(1:K),cen)
N <- Nt

n <- 100
begin_point1 <- runif(n,0.005,0.2)
begin_point2 <- runif(n,100,300)
begin_point <- cbind(begin_point1,begin_point2)
try_result <- rep(0,n)
for(i in 1:n){
  fit <-  optim(par = begin_point[i,],fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.8,10000),method="L-BFGS-B",control=list(fnscale=-1))
  try_result[i] <- BetaGeometricLikehood(fit$par)
  
}
idx <- which.max(try_result)
begin <- begin_point[idx,]

BetaEst <- ParaBetaEstimateFunction(Nt)
bootdata <- data.frame(N=N,cen=cen)
# NPMLEst_boot <- boot(data.frame(N=N,cen=cen),NPMLEstimateFunction_mean,Rboot)
# PooledEst_boot <- boot(N,PooledEstimateFunction,Rboot)
# StratEst_boot <- boot(N,StratEstimateFunction,Rboot)

NPMLEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  NPMLEst_boot[i] <- NPMLEstimateFunction_mean(N,c(1:K),cend)
}

PooledEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  PooledEst_boot[i] <- PooledEstimateFunction(N,c(1:K),cend)
}

StratEst_boot <- rep(0,Rboot)

for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cend <- cen[ind]
  StratEst_boot[i] <- StratEstimateFunction(N,c(1:K),cend)
}
Rboot <- 300
BetaEst_boot <- rep(0,Rboot)

cen_original <- cen
for(i in 1:Rboot){
  print(i)
  ind <- sample(c(1:length(Nt)),length(Nt),replace = T)
  N <- Nt[ind]
  cen <- cen_original[ind]
  n <- 100
  begin_point1 <- runif(n,0.005,0.2)
  begin_point2 <- runif(n,100,300)
  begin_point <- cbind(begin_point1,begin_point2)
  try_result <- rep(0,n)
  for(j in 1:n){
    fit <-  optim(par = begin_point[j,],fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.8,10000),method="L-BFGS-B",control=list(fnscale=-1))
    try_result[j] <- BetaGeometricLikehood(fit$par)
    
  }
  idx <- which.max(try_result)
  begin <- begin_point[idx,]
  
  BetaEst_boot[i] <- ParaBetaEstimateFunction(N)
  
}
cen <- cen_original
#BetaEst_boot_f <- BetaEst_boot[BetaEst_boot<0.1]

places <- 0
NPMLEstRound <- round(NPMLEst*10^4,places)
StratEstRound <- round(StratEst*10^4,places)
PooledEstRound <- round(PooledEst*10^4,places)
BetaEstRound <- round(BetaEst*10^4,places)
NPMLEstLow <- round(quantile(NPMLEst_boot,0.025)*10^4,places)
NPMLEstHigh <- round(quantile(NPMLEst_boot,0.975)*10^4,places)
StratEstLow <- round(quantile(StratEst_boot,0.025)*10^4,places)
StratEstHigh <- round(quantile(StratEst_boot,0.975)*10^4,places)
PooledEstLow <- round(quantile(PooledEst_boot,0.025)*10^4,places)
PooledEstHigh <- round(quantile(PooledEst_boot,0.975)*10^4,places)
BetaEstLow <- round(quantile(BetaEst_boot,0.025)*10^4,places)
BetaEstHigh <- round(quantile(BetaEst_boot,0.975)*10^4,places)

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


