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
    
    keep<-sample(1:5000,M)
    
    outtrt<-as.matrix(MCMClogit(Ytrt~pssplinetrt+xorttrt,burnin = 10000,mcmc=5000,verbose=5000))[keep,]
    outctl<-as.matrix(MCMClogit(Yctl~pssplinectl+xortctl,burnin = 10000,mcmc=5000,verbose=5000))[keep,]
    
    ## step 4: impute the missing potential outcomes for each sample from the posterior##
    imptrt1<-list()
    imptrt2<-list()
    impctl1<-list()
    impctl2<-list()
    for (i in 1:M){
      imptrt1<-c(imptrt1,list(rowSums(matrix(outtrt[i,],nrow=length(Ytrt),ncol=ncol(outtrt),byrow=T)*cbind(1,pssplinetrt,xorttrt))))
      imptrt2<-c(imptrt2,list(rowSums(matrix(outctl[i,],nrow=length(Ytrt),ncol=ncol(outctl),byrow=T)*cbind(1,pssplinetrt,xorttrt))))
      impctl1<-c(impctl1,list(rowSums(matrix(outtrt[i,],nrow=length(Yctl),ncol=ncol(outtrt),byrow=T)*cbind(1,pssplinectl,xortctl))))
      impctl2<-c(impctl2,list(rowSums(matrix(outctl[i,],nrow=length(Yctl),ncol=ncol(outctl),byrow=T)*cbind(1,pssplinectl,xortctl))))
    }
    imptrt1<-as.matrix(as.data.frame(imptrt1))
    imptrt2<-as.matrix(as.data.frame(imptrt2))
    impctl1<-as.matrix(as.data.frame(impctl1))
    impctl2<-as.matrix(as.data.frame(impctl2))
    
    ## step 6: estimate the treatment effect and its sample variance for each dataset ##
    indtrteff<-t(rbind((exp(imptrt1)/(1+exp(imptrt1)))-(exp(imptrt2)/(1+exp(imptrt2))),
                     (exp(impctl1)/(1+exp(impctl1)))-(exp(impctl2)/(1+exp(impctl2)))))
    
    iceavg<-apply(indtrteff,2,mean)
    icelw<-apply(indtrteff,2,quantile,probs=0.025)
    icehi<-apply(indtrteff,2,quantile,probs=0.975)
    ace<-rowMeans(indtrteff)
    
    #gammahat<-mean(avgtrteff)
    #gammaci<-quantile(avgtrteff,probs=c(.025,.975))
    
    return(list(use,iceavg,icelw,icehi,ace))
  }
}