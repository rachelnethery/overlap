bartspl<-function(datall,RO,nburn=10000,nsamp=5000){
  
  ## this function implements BART+SPL method (for continuous outcomes)
  ## first column in datall should be the observed outcome variable
  ## second column should be the binary exposure indicator
  ## third variable should be the PS or confounder on which to base the non-overlap
  ## any other variables to be included in the model are in the fourth column and beyond
  
  names(datall)[1:3]<-c('Yobs','x','ps')
  if (ncol(datall)>3) names(datall)[4:ncol(datall)]<-paste('u',1:(ncol(datall)-3),sep='')

  ## create an overlap dataset ##
  datov<-datall[which(RO==1),]
  ## create a test overlap dataset used to predict counterfactuals with BART ##
  datovtest<-datov
  datovtest$x<-1-datovtest$x
  ## create an overlap dataset for plugging into the spline model ##
  sp_ind<-which(datov$ps>quantile(datov$ps,.05) & datov$ps<quantile(datov$ps,.95))
  datov_sp<-datov[sp_ind,]
  ## create a non-overlap dataset ##
  datno<-datall[which(RO==0),]
  ## dataset containing units in the RN with E=0 (untreated) ##
  datno0<-datall[which(RO==0 & datall$x==0),]
  ## dataset containing units in the RN with E=1 (treated) ##
  datno1<-datall[which(RO==0 & datall$x==1),]
  
  ## for everyone in the RN, find the distance from their PS to the nearest PS in the RO ##
  ROdist<-NULL
  for (i in 1:nrow(datno)){
    psi<-datno$ps[i]
    ROdist<-c(ROdist,min(abs(datov$ps-psi)))
  }

  ## initialize the dbarts sampler to implement BART in the range of overlap ##
  dbfit<-dbarts(formula=Yobs~.,data=datov,test=datovtest[,2:ncol(datovtest)])
  
  ## hyperparameters for the spline ##
  p_spl1<-8+(ncol(datall)-3)
  p_spl0<-8+(ncol(datall)-3)
  mu0<-matrix(0,nrow=p_spl1,ncol=1)
  Sigma0<-10000*diag(p_spl1)
  Sigma0_inv<-solve(Sigma0)
  a0<-1
  b0<-1
  
  ## initialize parameters etc ##
  delta<-matrix(0,nrow=nrow(datov),ncol=1)
  Y1s<-datov_sp$Yobs
  Y0s<-datov_sp$Yobs
  beta1<-matrix(0,nrow=p_spl1,ncol=1)
  sigma_spl1<-1
  beta0<-matrix(0,nrow=p_spl0,ncol=1)
  sigma_spl0<-1
  
  ## matrices to store posterior samples ##
  # delta_save<-matrix(NA,nrow=nsamp,ncol=nrow(datov))
  delta_star_save<-matrix(NA,nrow=nsamp,ncol=nrow(datall))
  # sigma_bart_save<-rep(NA,nsamp)
  # beta_save<-matrix(NA,nrow=nsamp,ncol=p_spl)
  # sigma_spl_save<-rep(NA,nsamp)
  
  ## run the sampler ##
  for (i in 1:(nburn+nsamp)){
    
    ########################################################################
    ## 1. run the dbarts sampler once and take a sample from the BART ppd ##
    ########################################################################
    
    temp<-dbfit$run(numBurnIn=0,numSamples=1)
    ppd_test<-rnorm(n=length(temp$test),mean=temp$test,sd=temp$sigma)
    
    ## form the predicted individual causal effects in the RO from the BART predictions ##
    delta[which(datov$x==1),]<-datov$Yobs[which(datov$x==1)]-ppd_test[which(datov$x==1)]
    delta[which(datov$x==0),]<-ppd_test[which(datov$x==0)]-datov$Yobs[which(datov$x==0)]
    delta_sp<-delta[sp_ind,]
    
    ## update Y1s ##
    foosp<-temp$test[sp_ind]
    Y1s[which(datov_sp$x==0)]<-foosp[which(datov_sp$x==0)]
    Y0s[which(datov_sp$x==1)]<-foosp[which(datov_sp$x==1)]
    
    ## fit spline to the estimated causal effects in the RO from BART ##
    delta_star<-rep(NA,sum(RO==0))
    bs_ps<-rcspline.eval(datov_sp$ps,knots=quantile(datov_sp$ps,probs=c(.1,.25,.5,.75,.9)),inclx = T)
    
    ## spline for treated units (E=1) in the RN ##
    if (nrow(datno1)>0){
      
      ####################################
      ## 2. run the spline sampler once ##
      ####################################
      
      bs_y1s<-rcspline.eval(Y1s,knots=quantile(Y1s,probs=c(.2,.4,.6,.8)),inclx = T)
      if (ncol(datov_sp)<=3){
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y1s))
      }
      else{
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y1s,datov_sp[,4:ncol(datov_sp)]))
      }
      
      ## sample the betas ##
      Vbeta<-solve(Sigma0_inv+((1/sigma_spl1)*t(X_spl)%*%X_spl))
      Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl1)*t(X_spl)%*%delta_sp))
      beta1<-matrix(mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl1,ncol=1)
      
      ## sample the sigma_spls ##
      a<-a0+(nrow(datov_sp)/2)
      b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta1))^2))
      sigma_spl1<-rinvgamma(n=1,shape=a,scale=b)
      
      ######################################################################
      ## 3. draw from the posterior predictive for the non-overlap region ##
      ######################################################################
      bs_ps_star<-rcspline.eval(datno1$ps,knots=quantile(datov_sp$ps,probs=c(.1,.25,.5,.75,.9)),inclx = T)
      bs_y1s_star<-rcspline.eval(datno1$Yobs,knots=quantile(Y1s,probs=c(.2,.4,.6,.8)),inclx=T)
      if (ncol(datno)<=3){
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y1s_star))
      }
      else{
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y1s_star,datno1[,4:ncol(datno1)]))
      }
      Eppd<-X_spl_star%*%beta1
      delta_star[which(datno$x==1)]<-rnorm(n=nrow(datno1),mean=Eppd,sd=sqrt(sigma_spl1+ROdist[which(datno$x==1)]*10*(max(delta)-min(delta))))
    }
    
    ## spline for untreated units (E=0) in the RN ##
    if (nrow(datno0)>0){
      
      ####################################
      ## 2. run the spline sampler once ##
      ####################################
      
      bs_y0s<-rcspline.eval(Y0s,knots=quantile(Y0s,probs=c(.2,.4,.6,.8)),inclx=T)
      if (ncol(datov_sp)<=3){
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y0s))
      }
      else{
        X_spl<-as.matrix(cbind(1,bs_ps,bs_y0s,datov_sp[,4:ncol(datov_sp)]))
      }
      
      ## sample the betas ##
      Vbeta<-solve(Sigma0_inv+((1/sigma_spl0)*t(X_spl)%*%X_spl))
      Ebeta<-Vbeta%*%(Sigma0_inv%*%mu0+((1/sigma_spl0)*t(X_spl)%*%delta_sp))
      beta0<-matrix(mvrnorm(n=1,mu=Ebeta,Sigma=Vbeta),nrow=p_spl0,ncol=1)
      
      ## sample the sigma_spls ##
      a<-a0+(nrow(datov_sp)/2)
      b<-b0+((1/2)*sum((delta_sp-(X_spl%*%beta0))^2))
      sigma_spl0<-rinvgamma(n=1,shape=a,scale=b)
      
      ######################################################################
      ## 3. draw from the posterior predictive for the non-overlap region ##
      ######################################################################
      bs_ps_star<-rcspline.eval(datno0$ps,knots=quantile(datov_sp$ps,probs=c(.1,.25,.5,.75,.9)),inclx = T)
      bs_y0s_star<-rcspline.eval(datno0$Yobs,knots=quantile(Y0s,probs=c(.2,.4,.6,.8)),inclx=T)
      if (ncol(datno)<=3){
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y0s_star))
      }
      else{
        X_spl_star<-as.matrix(cbind(1,bs_ps_star,bs_y0s_star,datno0[,4:ncol(datno0)]))
      }
      Eppd<-X_spl_star%*%beta0
      delta_star[which(datno$x==0)]<-rnorm(n=nrow(datno0),mean=Eppd,sd=sqrt(sigma_spl0+ROdist[which(datno$x==0)]*10*(max(delta)-min(delta))))
    }
    
    
    ###############################################
    ## 4. if we're past burn-in, save the output ##
    ###############################################
    if (i>nburn){
      delta_star_save[(i-nburn),which(RO==1)]<-c(delta)
      delta_star_save[(i-nburn),which(RO==0)]<-delta_star
    }
    
  }
  
  ## use the Bayesian bootstrap to get ACE posterior ##
  ace_pd<-apply(delta_star_save,1,aceBB)
  
  return(list(ace_pd,apply(delta_star_save,2,mean),apply(delta_star_save,2,quantile,probs=.025),apply(delta_star_save,2,quantile,probs=.975)))
  
}