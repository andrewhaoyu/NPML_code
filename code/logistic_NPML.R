library(MASS)
library(nleqslv)
library(SQUAREM)
library(Rcpp)
library(car)


###E step calculation function
###let i represent person
###let j represent the population
###E_step_pp calculate the pr(N_i|Z_j) = (1-p_j)^(N_i-1)*(1-p_j)

cppFunction('NumericMatrix E_step_pp(NumericVector uu,double b,NumericVector x,NumericVector N){
            int KK = N.size();
            NumericMatrix pp(KK,KK);
            for(int i=0; i < KK; ++i){
            for(int j=0; j <KK; ++j){
            pp(i,j) = pow((1-exp(uu(j)+b*x(i))/(1+exp(uu(j)+b*x(i)))),(N(i)-1))*exp(uu(j)+b*x(i))/(1+exp(uu(j)+b*x(i)));
            }
            }
            return(pp);
            }')
######### E_step_ww perform the bayesian rule
######### Pr(Z_j|N_i) = Pr(N_i|Z_j)*Pr(Z_j)/(sum_j Pr(N_i|Z_j)*Pr(Z_j))
cppFunction('NumericMatrix E_step_ww(NumericMatrix pp, NumericVector w){
            int KK= w.size();
            NumericMatrix ww(KK,KK);
            double isum=0;
            for(int i=0; i< KK;++i ){
            isum=0;
            for(int j=0; j < KK; ++j){
            ww(i,j)=w(j)*pp(i,j);
            isum += ww(i,j); 
            }
            for(int j=0; j <KK; ++j){
            ww(i,j)=ww(i,j)/isum;
            }
            }
            return(ww);
}')


#### M step calculation function
gr_u <- function(u,N,ww_vec,b,x,K,observation){
  sum <- 0
  temp_Qu <- 1-N*logit_inver(u+b%*%x)
  sum <- crossprod(temp_Qu,ww_vec)
  return(sum^2)
  #  for(l in 1:observation){
  #   N_temp <- N[seq(l,K*observation,observation)]
  #   x_temp <- x[seq(l,K*observation,observation)]
  #   temp_Qu <- 1-N_temp*logit_inver(u+b*x_temp)
  #   sum <- crossprod(temp_Qu,ww_vec)+sum
  # }
  #return(sum^2)
}


gr_u_d <- function(u,N,ww_vec,b,x,K,observation){
  xb <- x%*%b
  temp_Qu <- 1-N*logit_inver(u+xb)
  sum1 <- crossprod(temp_Qu,ww_vec)
  temp_d_Qu <- N*logit_inver(u+xb)/(1+exp(u+xb))
  sum2 <- -crossprod(temp_d_Qu,ww_vec)
  sum <- 2*sum1*sum2
  
  return(sum)
  # for(l in 1:observation){
  #   N_temp <- N[seq(l,K*observation,observation)]
  #   x_temp <- x[seq(l,K*observation,observation)]
  #   temp_Qu <- 1-N_temp*logit_inver(u+b*x_temp)
  #   sum1 <- crossprod(temp_Qu,ww_vec)+sum1
  # }
  # for(l in 1:observation){
  #   N_temp <- N[seq(l,K*observation,observation)]
  #   x_temp <- x[seq(l,K*observation,observation)]
  #   temp_d_Qu <- N_temp*logit_inver(u+b*x_temp)/(1+exp(u+b*x_temp))
  #   sum2 <- -crossprod(temp_d_Qu,ww_vec)+sum2
  # }
  # 
  #sum <- 2*sum1*sum2
  
  #return(sum)
}
gr_b <- function(uu,N,ww,b,x,K,observation){
  u <- uu
  p <- length(b)
  t <- 0
  xb <- x%*%b
  for(k in 1:p){
    temp_Qb <- matrix(0,K,K)
    for(i in 1:K){
      for(j in 1:K){
        temp_Qb[i,j] <- x[i,k]-N[i]*x[i,k]*logit_inver(u[j]+xb)
      }
    }
    t <- sum(temp_Qb*ww)
    p[k] <- t
  }
  
  return(p)
  # for(l in 1:observation){
  #   N_temp <- N[seq(l,K*observation,observation)]
  #   x_temp <- x[seq(l,K*observation,observation)]
  #   temp_Qb <- matrix(0,K,K)
  #   for(i in 1:K){
  #     for(j in 1:K){
  #       temp_Qb[i,j] <- x_temp[i]-N_temp[i]*x_temp[i]*logit_inver(u[j]+b*x_temp[i])
  #     }
  #   }
  #   t <- sum(temp_Qb*ww)+t
  # }
  # return(t)
}

gr_b_d <- function(uu,N,ww,b,x,K,observation){
  u <- uu
  t <- 0
  for(l in 1:observation){
    N_temp <- N[seq(l,K*observation,observation)]
    x_temp <- x[seq(l,K*observation,observation)]
    temp_Qb <- matrix(0,K,K)
    for(i in 1:K){
      for(j in 1:K){
        temp_Qb[i,j] <- -N_temp[i]*x_temp[i]^2*logit_inver(u[j]+b*x_temp[i])/(1+exp(u[j]+b*x_temp[i]))
      }
    }
    t <- sum(temp_Qb*ww)+t
  }
  return(as.matrix(t))
}

#### NPML iteration

NPML_logistic_function_iter <- function(param,data){
  observation <- data[length(data)]
  K <- data[length(data)-1]
  N <- data[1:(observation*K)]
  x <- data[(observation*K+1):(length(data)-2)]
  uu <- param[1:K]
  w <- param[(K+1):(length(param)-1)]
  b <- param[length(param)]
  ww <- matrix(nrow=K,ncol=K)
  
  ### this epi was set for each M step
  epi <- tol
  pp <- matrix(1,K,K)
  for(l in 1:observation){
    N_temp <- N[seq(l,K*observation,observation)]
    x_temp <- x[seq(l,K*observation,observation)]
    pp <- E_step_pp(uu,b,x_temp,N_temp)*pp
  }
  
  ww <- E_step_ww(pp,w)
  for(i in 1:K){
    w[i] <- sum(ww[,i])/sum(ww)
  }
  ############
  #Mstep
  Niter_M <- 100
  
  
  for(i in 1:Niter_M){
    uu_temp <-uu
    b_temp <- b
    for(j in 1:K){
      ww_vec <- ww[,j]
      
      ans <-  optim(uu_temp[j],fn=gr_u,gr=gr_u_d,method ="L-BFGS-B",lower=(u_new_min-0.5),upper=(u_new_max+0.5),N=N,ww_vec=ww_vec,b=b,x=x,K=K,observation=observation)
      uu[j] <- ans$par
    }
    ans <- nleqslv(b_temp,
                   fn=function(y){gr_b(uu,N,ww,y,x,K,observation)},
                   jac=function(y){gr_b_d(uu,N,ww,y,x,K,observation)})
    b <- ans$x
    
    thres_M <- sum((uu-uu_temp)^2)+(b-b_temp)^2
    if(thres_M < epi){
      break
    }
    
  }
  
  intercept <- crossprod(uu,w)
  slope <- b
  if(abs(intercept+slope)>100){
    uu <- rep(1000,K)
    w <- rep(1/K,K)
    b <- 1000
  }
  
  param.new <- c(uu,w,b)
  return(param.new)
}



#### LogLikehood for NPML
NPMLEloglik <- function(param,data){
  observation <- data[length(data)]
  K <- data[length(data)-1]
  N <- data[1:(observation*K)]
  x <- data[(observation*K+1):(length(data)-2)]
  uu <- param[1:K]
  w <- param[(K+1):(length(param)-1)]
  b <- param[length(param)]
  
  KK <- length(data)/2
  N <- data[1:KK]
  x <- data[(KK+1):length(data)]
  UU <- param[1:KK]
  pp <- param[(KK+1):(length(param)-1)]
  b <- param[length(param)]
  
  result <- 0
  for(l in 1:observation){
    N_temp <- N[seq(l,K*observation,observation)]
    x_temp <- x[seq(l,K*observation,observation)]
    for(i in 1:K){
      result <- result + log(sum(w,dgeom((N_temp[i]-1),logit_inver(uu+b*x_temp[i]))))
    }
  }
  
  return(-result)
}



# NPML_logistic_function <- function(N,x,tol,maxit){
#   converge <- F
#   K <- length(N)
#   w <- rep(1/K,K)
#   ww <- matrix(1/K,K,K)
#   uu <- seq(u_new_min,u_new_max,(u_new_max-u_new_min)/(K-1))
#   pp <- matrix(0,K,K)
#   b <- beta
#   param.start <- c(uu,w,b)
#   data <- c(N,x)
#   param.old <- param.start
#   for(i in 1:maxit){
#     uu_old <- param.old[1:K]
#     pp_old <- param.old[(K+1):(2*K)]
#     b_old <- param.old[(2*K+1)]
#     mean_u_old <- crossprod(uu_old,pp_old)
#     sd_u_old <- crossprod((uu_old-mean_u_old)^2,pp_old)
#     param.new <- NPML_logistic_function_iter(param.old,data)
#     uu_new <- param.new[1:K]
#     pp_new <- param.new[(K+1):(2*K)]
#     b_new <- param.new[(2*K+1)]
#     mean_u_new <- crossprod(uu_new,pp_new)
#     sd_u_new <- crossprod((uu_new-mean_u_new)^2,pp_new)
#     param.old <- param.new
#     threshold <- abs(mean_u_old-mean_u_new)/(abs(mean_u_new)+0.1)+
#       abs(sd_u_old-sd_u_new)/(abs(sd_u_new)+0.1)+
#       abs(b_old-b_new)/(abs(b_new)+0.1)
#     if(threshold<tol){
#       converge <- T
#       break
#     }
#     
#   }
#   
#   # f1 <- squarem1(par = param.start, 
#   #               objfn = NPMLEloglik, 
#   #               fixptfn = NPML_logistic_function_iter, 
#   #               data = c(N,x),
#   #               control  = list(maxiter = maxit,tol=tol)
#   # )
# 
# 
#   # return(c(crossprod(UU,pp),b,converge))
#   return(c(mean_u_new,b_new,converge))
#   
# }

NPML_logistic_function <- function(N,x,tol,maxit,K,observation){
  print(indicator)
  indicator <<- indicator+1
  w <- rep(1/K,K)
  ww <- matrix(1/K,K,K)
  uu <- seq(u_new_min,u_new_max,(u_new_max-u_new_min)/(K-1))
  pp <- matrix(0,K,K)
  b <- beta
  param.start <- c(uu,w,b)
  data <- c(N,x,K,observation)
  param.old <- param.start
  
  #   uu_old <- param.old[1:K]
  #   pp_old <- param.old[(K+1):(2*K)]
  #   b_old <- param.old[(2*K+1)]
  #   mean_u_old <- crossprod(uu_old,pp_old)
  #   sd_u_old <- crossprod((uu_old-mean_u_old)^2,pp_old)
  #   param.new <- NPML_logistic_function_iter(param.old,data)
  #   uu_new <- param.new[1:K]
  #   pp_new <- param.new[(K+1):(2*K)]
  #   b_new <- param.new[(2*K+1)]
  #   mean_u_new <- crossprod(uu_new,pp_new)
  #   sd_u_new <- crossprod((uu_new-mean_u_new)^2,pp_new)
  #   threshold <- abs(mean_u_old-mean_u_new)/(abs(mean_u_new)+0.1)+
  #     abs(sd_u_old-sd_u_new)/(abs(sd_u_new)+0.1)+
  #     abs(b_old-b_new)/(abs(b_new)+0.1)
  #   if(threshold<tol){
  #     converge <- T
  #     break
  #   }
  #   
  # }
  # 
  f1 <- tryCatch(squarem1(par = param.start,
                          objfn = NPMLEloglik,
                          fixptfn = NPML_logistic_function_iter,
                          data = data,
                          control  = list(maxiter = maxit,tol=tol)
  ),
  error=function(e){
    result <- list(par=rep(1,(2*K+1)),convergence=F)
    return(result)
  })
  
  
  UU <- f1$par[1:K]
  pp <- f1$par[(K+1):(2*K)]
  b <- f1$par[2*K+1]
  converge <- f1$convergence
  intercept <- crossprod(UU,pp)
  if(abs(intercept+b)>100){
    converge <- F
  }
  return(c(crossprod(UU,pp),b,converge))
  # return(c(mean_u_new,b_new,converge))
  
}


squarem1 <- function(par, fixptfn, objfn, ... , control=list()) {
  # par = starting value of parameter vector
  # fixptfn = fixed-point iteration F(x)
  # for which the solution: F(x*) = x* is sought
  # objfn = underlying objective function which is minimized at x*
  #
  #
  K <- (length(par)-1)/2
  
  
  control.default <- list(K=1, square=TRUE, method=3, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,
                          tol=1.e-07, maxiter=1500, trace=FALSE)
  
  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)
  
  #####
  # method = 1, 2, or 3, indicating the type of steplength to be used
  # A key parameter is `step.min0'. This can be either positive or negative if the eigenvalues of the Jacobian of fixed-point mapping are all positive;
  # We set it to +1 for contraction mappings such as EM and MM algorithms
  # Must be negative, e.g. -1, if the fixed-point mapping can have negative eiegnvalues
  #####
  # parameter "objfn.inc" dictates the amount of non-monotonicity in the objective function
  # setting objfn.inc=0 would enforce monotonicity, whereas objfn.inc=Inf would be a non-monotonic scheme
  # The defalut objfn.inc=1 would enforce monotonicity far from solution, but allows for non-monotonicity closer to solution
  #
  method <- ctrl$method
  maxiter <- ctrl$maxiter
  tol <- ctrl$tol
  step.min <- ctrl$step.min0
  step.max0 <- ctrl$step.max0
  step.max <- ctrl$step.max0
  mstep <- ctrl$mstep
  objfn.inc <- ctrl$objfn.inc
  trace <- ctrl$trace
  
  if (trace) cat("Squarem-1 \n")
  
  if (missing(objfn)) stop("\n squarem2 should be used if objective function is not available \n\n")
  
  iter <- 1
  objval <- rep(NA,1)
  p <- par
  
  lold <- objfn(p, ...)
  leval <- 1
  if (trace) cat("Objective fn: ", lold, "\n")
  feval <- 0
  conv <- TRUE
  
  while (feval < maxiter) {   
    
    
    extrap <- TRUE
    p1 <- try(fixptfn(p, ...),silent=TRUE)
    
    feval <- feval + 1
    if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) stop("Error in function evaluation")
    q1 <- p1 - p
    sr2 <- crossprod(q1)
    uu <- p[1:K]
    pp <- p[(K+1):(2*K)]
    b <- p[2*K+1]
    u_mean <- crossprod(uu,pp)
    u_sd <- crossprod((uu-u_mean)^2,pp)
    uu1 <- p1[1:K]
    pp1 <- p1[(K+1):(2*K)]
    b1 <- p1[2*K+1]
    u_mean1 <- crossprod(uu1,pp1)
    u_sd1 <- crossprod((uu1-u_mean1)^2,pp1)
    set_tol1 <- (abs(u_mean1-u_mean))/(abs(u_mean1)+0.1)+
      (abs(b-b1))/(abs(b1)+0.1)
    set_tol1_temp <- sum((p-p1)^2)
    set_tol1 <- max(set_tol1,set_tol1_temp)
    # +(abs(u_sd-u_sd1))/(abs(u_sd1)+0.1)
    if (set_tol1 < tol) break
    
    p2 <- try(fixptfn(p1, ...),silent=TRUE)
    feval <- feval + 1
    if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) stop("Error in function evaluation")
    
    q2 <- p2 - p1
    sq2 <- sqrt(crossprod(q2))
    uu2 <- p2[1:K]
    pp2 <- p2[(K+1):(2*K)]
    b2 <- p2[2*K+1]
    u_mean2 <- crossprod(uu2,pp2)
    u_sd2 <- crossprod((uu2-u_mean2)^2,pp2)
    set_tol2 <- (abs(u_mean2-u_mean1))/(abs(u_mean2)+0.1)+
      (abs(b2-b1))/(abs(b2)+0.1)
    # +(abs(u_sd2-u_sd1))/(abs(u_sd2)+0.1)
    set_tol2_temp <- sum((p2-p1)^2)
    set_tol2 <- max(set_tol2,set_tol2_temp)
    if (set_tol2 < tol) break
    sv2 <- crossprod(q2-q1)
    srv <- crossprod(q1, q2-q1)
    
    alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
    
    alpha <- max(step.min, min(step.max, alpha))
    p.new <- p + 2*alpha*q1 + alpha^2*(q2-q1)
    if (abs(alpha - 1) > 0.01 ) {
      p.new <- try(fixptfn(p.new , ...),silent=TRUE)   # stabilization step
      feval <- feval + 1
    }
    
    if (class(p.new) == "try-error" | any(is.nan(p.new))) {
      p.new <- p2
      lnew <- try(objfn(p2, ...), silent=TRUE)
      leval <- leval + 1
      if (alpha == step.max) step.max <- max(step.max0, step.max/mstep)
      alpha <- 1
      extrap <- FALSE
    } else {
      if (is.finite(objfn.inc)) {
        lnew <- try(objfn(p.new, ...), silent=TRUE)
        leval <- leval + 1
      } else lnew <- lold
      if (class(lnew) == "try-error" | is.nan(lnew) |
          (lnew > lold + objfn.inc)) {
        p.new <- p2
        lnew <- try(objfn(p2, ...), silent=TRUE)
        leval <- leval + 1
        if (alpha==step.max) step.max <- max(step.max0, step.max/mstep)
        alpha <- 1
        extrap <- FALSE
      }
    }
    
    if (alpha == step.max) step.max <- mstep*step.max
    if (step.min < 0 & alpha == step.min) step.min <- mstep*step.min
    
    p <- p.new
    if (!is.nan(lnew)) lold <- lnew
    if (trace) cat("Objective fn: ", lnew, "  Extrapolation: ", extrap, "  Steplength: ", alpha, "\n")
    iter <- iter+1
    
  }
  if (feval > maxiter) conv <- FALSE
  if (is.infinite(objfn.inc)) {
    lold <- objfn(p, ...)
    leval <- leval + 1
  }
  
  
  return(list(par=p, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv))
  
}

