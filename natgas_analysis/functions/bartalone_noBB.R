bartalone<-function(xtr,ytr,xte,nburn=10000,nsamp=5000){
  
  bartps<-bart(x.train=xtr,y.train=ytr,x.test=xte,nskip=nburn,ndpost=nsamp)
  ppd_test<-t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=bartps$sigma)))
  ppd_test_mean<-apply(ppd_test,2,mean)
  
  ## individual causal effects ##
  iceavg<-rep(NA,length(ytr))
  iceavg[which(xtr[,1]==1)]<-ytr[which(xtr[,1]==1)]-ppd_test_mean[which(xtr[,1]==1)]
  iceavg[which(xtr[,1]==0)]<-ppd_test_mean[which(xtr[,1]==0)]-ytr[which(xtr[,1]==0)]
  icelw<-apply(ppd_test,2,quantile,probs=0.025)
  icehi<-apply(ppd_test,2,quantile,probs=0.975)
  
  ## average causal effects ##
  ppd_ice<-matrix(NA,nrow=nrow(ppd_test),ncol=length(ytr))
  for (j in 1:length(ytr)){
    if (xtr[j,1]==1) ppd_ice[,j]<-ytr[j]-ppd_test[,j]
    else ppd_ice[,j]<-ppd_test[,j]-ytr[j]
  }
  ## get ACE posterior using the Bayesian bootstrap ##
  #ace_pd<-apply(ppd_ice,1,aceBB)
  temp<-matrix(NA,nrow=nrow(ppd_test),ncol=ncol(ppd_test))
  temp[,which(xtr[,1]==1)]<-matrix(rep(ytr[which(xtr[,1]==1)],nrow(ppd_test)),nrow=nrow(ppd_test),byrow = T)-ppd_test[,which(xtr[,1]==1)]
  temp[,which(xtr[,1]==0)]<-ppd_test[,which(xtr[,1]==0)]-matrix(rep(ytr[which(xtr[,1]==0)],nrow(ppd_test)),nrow=nrow(ppd_test),byrow = T)
  ace_pd<-rowMeans(temp)
  
  return(list(iceavg,icelw,icehi,ace_pd))
}