aceBB<-function(x,nboot=250){
  
  diffpotmat<-matrix(rep(x,nboot),nrow=nboot,byrow=T)
  
  dirichlet_sample <- matrix( rexp(length(x) * nboot, 1) , nrow=nboot, byrow = TRUE)
  dirichlet_sample <- dirichlet_sample / rowSums(dirichlet_sample)
  
  posteriorace<-rowSums(dirichlet_sample*diffpotmat)
  
  ## sample from the posterior ##
  MCMCace<-sample(posteriorace,size=1)
  
  return(MCMCace)
}