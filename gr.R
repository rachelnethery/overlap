## gutman and rubin's method ##
gr<-function(Y,trt,ps,X,M,qps){
  ## step 1: create subclasses based on the PS ##
  ps.<-as.numeric(cut(ps,qps,include.lowest = T,right=F))
  
  use<-1
  for (i in 1:length(unique(ps.))){
    if (length(which(ps.==i & trt==0))<3 | length(which(ps.==i & trt==1))<3) use<-0
  }
  
  if (use==0){
    return(list(use))
  }
  else{
    ## steps 2 & 3: estimate splines separately for treated and controls & sample M times from the posterior ##
    pstrans<-log(ps/(1-ps))
    psspline<-ns(pstrans,knots = qps[-c(1,length(qps))])
    xort<-NULL
    for (i in 1:ncol(X)){
      xort<-rbind(xort,lm(pstrans~X[,i])$resid)
    }
    xort<-t(xort)
    
    Ytrt<-Y[which(trt==1)]
    pssplinetrt<-psspline[which(trt==1)]
    xorttrt<-xort[which(trt==1)]
    Yctl<-Y[which(trt==0)]
    pssplinectl<-psspline[which(trt==0)]
    xortctl<-xort[which(trt==0)]
    
    keep<-sample(1:10000,M)
    
    outtrt<-as.matrix(MCMCregress(Ytrt~pssplinetrt+xorttrt,burnin = 1000,mcmc=10000,verbose=F))[keep,]
    outctl<-as.matrix(MCMCregress(Yctl~pssplinectl+xortctl,burnin = 1000,mcmc=10000,verbose=F))[keep,]
    
    ## step 4: impute the missing potential outcomes for each sample from the posterior##
    imptrt<-list()
    impctl<-list()
    for (i in 1:M){
      imptrt<-c(imptrt,list(rowSums(matrix(outctl[i,(-ncol(outctl))],nrow=length(Ytrt),ncol=ncol(outctl)-1,byrow=T)*cbind(1,pssplinetrt,xorttrt))))
      impctl<-c(impctl,list(rowSums(matrix(outtrt[i,(-ncol(outtrt))],nrow=length(Yctl),ncol=ncol(outtrt)-1,byrow=T)*cbind(1,pssplinectl,xortctl))))
    }
    imptrt<-as.matrix(as.data.frame(imptrt))
    impctl<-as.matrix(as.data.frame(impctl))
    
    ## step 6: estimate the treatment effect and its sample variance for each dataset ##
    indtrteff<-t(rbind(matrix(Ytrt,nrow=length(Ytrt),ncol=M,byrow=F)-imptrt,
                     impctl-matrix(Yctl,nrow=length(Yctl),ncol=M,byrow=F)))
    
    iceavg<-apply(indtrteff,2,mean)
    icelw<-apply(indtrteff,2,quantile,probs=0.025)
    icehi<-apply(indtrteff,2,quantile,probs=0.975)
    ace<-rowMeans(indtrteff)
    
    #gammahat<-mean(avgtrteff)
    #gammaci<-quantile(avgtrteff,probs=c(.025,.975))
    
    return(list(use,iceavg,icelw,icehi,ace))
  }
}