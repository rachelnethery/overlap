bartalone<-function(xtr,ytr,xte,nburn=10000,nsamp=5000){
  
  bartps<-bart(x.train=xtr,y.train=ytr,x.test=xte,nskip=nburn,ndpost=nsamp)
  # ppd_test<-pnorm(t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=1))))
  # ppd_train<-pnorm(t(apply(bartps$yhat.train,1,function(x) rnorm(n=length(x),mean=x,sd=1))))
  ppd_test<-pnorm(bartps$yhat.test)
  ppd_train<-pnorm(bartps$yhat.train)
  ppd_test_mean<-apply(ppd_test,2,mean)
  ppd_train_mean<-apply(ppd_train,2,mean)
  
  ## individual causal effects ##
  iceavg<-rep(NA,length(ytr))
  iceavg[which(xtr[,1]==1)]<-ppd_train_mean[which(xtr[,1]==1)]-ppd_test_mean[which(xtr[,1]==1)]
  iceavg[which(xtr[,1]==0)]<-ppd_test_mean[which(xtr[,1]==0)]-ppd_train_mean[which(xtr[,1]==0)]
  
  ## average causal effects ##
  ppd_ice<-matrix(NA,nrow=nrow(ppd_test),ncol=length(ytr))
  for (j in 1:length(ytr)){
    if (xtr[j,1]==1) ppd_ice[,j]<-ppd_train[,j]-ppd_test[,j]
    else ppd_ice[,j]<-ppd_test[,j]-ppd_train[,j]
  }
  icelw<-apply(ppd_ice,2,quantile,probs=0.025)
  icehi<-apply(ppd_ice,2,quantile,probs=0.975)
  
  ## get ACE posterior using the Bayesian bootstrap ##
  ace_pd<-apply(ppd_ice,1,aceBB)
  
  return(list(iceavg,icelw,icehi,ace_pd))
}