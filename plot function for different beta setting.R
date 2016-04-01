theta_bar_vector=c(0.1,0.2,0.3,0.4,0.5)
M_vector=c(0.05,0.25,0.5,1,1.5,2,5,10,25,10000)
K_vector=c(10,20,50,100,150)
  for(i1 in 1:length(theta_bar_vector)){
    for(i2 in 1:length(M_vector)){
      for(i3 in 1:length(K_vector)){
        places=4
        R=200
        theta_bar=theta_bar_vector[i1]
        
        
        M=M_vector[i2]
        alpha=theta_bar*M
        beta=M*(1-theta_bar)
        K=K_vector[i3]
        cuts=1:K/(K+1)
        theta=qbeta(cuts,alpha,beta)
        delta=0.005
        theta[which(theta<delta)]=delta
        plot(cuts,theta,main=paste("theta=",theta_bar,"M=",M,"K=",K),ylim=c(0,1))
      }
    }
  }



# theme -------------------------------------------------------------------

# sectionII ---------------------------------------------------------------


