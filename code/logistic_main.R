rm(list=ls()) #clear the environment
commanarg <- commandArgs(trailingOnly = T)
i1 <- as.numeric(commanarg[1])
i2 <- as.numeric(commanarg[2])
i3 <- as.numeric(commanarg[3])
i4 <- as.numeric(commanarg[4])
# method function ---------------------------------------------------------
source("./code/logistic_support.R")
source("./code/logistic_NPML.R")
R <- 400
observation=1
K_vector <- c(16,32,48,64,80,96,112)
sd_u_vector <- c(0,1,3,5)
r_vector <- c(0,0.4,0.8)
beta_vector <- c(0,0.5,1)
set.seed(123)
logistic_result <- NULL


indicator <- 1
K <- K_vector[i1]
sd_u <- sd_u_vector[i2]
r <- r_vector[i3]
beta <- beta_vector[i4]
mu_u <- 0
u_new_min <- logit(0.02)+beta*1.5
u_new_max <- logit(0.98)-beta*1.5
data <- generating(mu_u,sd_u,K,r,beta,R,observation)
u_new <- data$u_new
x <- data$x
theta <- data$theta
NTL <- data$NTL


param.start <- c(mean(u_new),beta)
tol <- 1e-07
maxit <- 500
standard_logistic_results <- apply(NTL,1,function(y){standard_logistic(y,x,param.start,tol,maxit)})
idx_logistic <- standard_logistic_results[3,]

true_u <- mean(u_new)
NPML_logistic_function(NTL[1,],x,tol,maxit,K,observation)
NPML_logistic_results <- apply(NTL,1,function(y){NPML_logistic_function(y,x,tol,maxit,K,observation)})
idx_NPML <- NPML_logistic_results[3,]
result_table <- table(idx_logistic,idx_NPML)
idx_converge <- which(idx_logistic==1&idx_NPML==1)
standard_logistic_u_bias2 <- (mean(standard_logistic_results[1,idx_converge])-true_u)^2
standard_logistic_u_var <- var(standard_logistic_results[1,idx_converge])
standard_logistic_u_mse <- mean((standard_logistic_results[1,idx_converge]-true_u)^2)

standard_logistic_b_bias2 <- (mean(standard_logistic_results[2,idx_converge])-beta)^2
standard_logistic_b_var <- var(standard_logistic_results[2,idx_converge])
standard_logistic_b_mse <- mean((standard_logistic_results[2,idx_converge]-beta)^2)


NPML_logistic_u_bias2 <- (mean(NPML_logistic_results[1,idx_converge])-true_u)^2
NPML_logistic_u_var <- var(NPML_logistic_results[1,idx_converge])
NPML_logistic_u_mse <- mean((NPML_logistic_results[1,idx_converge]-true_u)^2)
NPML_logistic_b_bias2 <- (mean(NPML_logistic_results[2,idx_converge])-beta)^2
NPML_logistic_b_var <- var(NPML_logistic_results[2,idx_converge])
NPML_logistic_b_mse <- mean((NPML_logistic_results[2,idx_converge]-beta)^2)



mse_results <- c(standard_logistic_u_bias2,
                 standard_logistic_u_var,
                 standard_logistic_u_mse,
                 standard_logistic_b_bias2,
                 standard_logistic_b_var,
                 standard_logistic_b_mse,
                 NPML_logistic_u_bias2,
                 NPML_logistic_u_var,
                 NPML_logistic_u_mse,
                 NPML_logistic_b_bias2,
                 NPML_logistic_b_var,
                 NPML_logistic_b_mse,
                 K,
                 sd_u,
                 r,
                 beta)
final_result <- list(mse_results=mse_results,
                     result_table=result_table,
                     standard_logistic_results = standard_logistic_results,
                     NPML_logistic_results = NPML_logistic_results
)
save(final_result,file=paste0("/users/hzhang1/R/project_Tom/logistic_regression/simulation11/mse_",K,"_",sd_u,"_",r,"_",beta))
#save(final_result,file=paste0("mse_",K,"_",sd_u,"_",r,"_",beta))


