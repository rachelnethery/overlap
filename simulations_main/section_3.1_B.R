##########################################
## THIS CODE IMPLEMENTS SIMULATION 3.1A ##
##########################################

## read in command line arguments ##
args<-commandArgs(TRUE)
wd<-args[1] #working directory#
simnum<-as.numeric(args[2]) #simulation number#
reps<-as.numeric(args[3]) #number of reps to be implemented#
N<-as.numeric(args[4]) #total N#
extr<-as.numeric(args[5]) #parameter controlling degree of non-overlap (called c in manuscript)#
a<-as.numeric(args[6]) #a parameter in non-overlap definition#
b<-as.numeric(args[7]) #b parameter in non-overlap definition#

## run simulation ##
setwd(wd)
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
source('functions/gr.R')
source('functions/bartspl.R')
source('functions/bartalone.R')
source('functions/aceBB.R')
source('functions/pw_overlap.R')

set.seed(simnum)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

save_true<-list()
save_tgr_mis<-list()
save_ugr_mis<-list()
save_tbart_mis<-list()
save_ubart_mis<-list()
save_rn_mis<-list()


simscomplete<-0
while (simscomplete<=(reps-1)){
  
  ## generate confounders (1 binary, 1 continuous) ##
  u1<-c(rbinom(n=n1,size=1,prob=.5),rbinom(n=n0,size=1,prob=.4))
  u2<-c(rnorm(n=n1,mean=2+extr,sd=2),rnorm(n=n0,mean=1,sd=1))
  
  ## construct misspecified PS estimate ##
  emod<-glm(x~u1+u2,family = binomial())
  ps_mis<-emod$fitted
  
  ## find RO and RN ## 
  RO_mis<-pw_overlap(ps=ps_mis,E=x,a=a,b=b)
  order_ps_mis<-ps_mis[order(ps_mis)]
  order_RO_mis<-RO_mis[order(ps_mis)]
  ll_mis<-max(which(order_RO_mis==1))
  
  ## skip over datasets that have no non-overlap in right tail or too much non-overlap in the left tail (RNs with >10 PS's) ##
  if (ll_mis<N & sum(1-order_RO_mis[1:ll_mis])<=10){
    
    ## ignore small non-overlap in the left tail (add all left tail to the RO) ##
    temp1<-order_ps_mis[1:ll_mis]
    ps_addRO<-temp1[which(order_RO_mis[1:ll_mis]==0)]
    RO_mis[which((ps_mis %in% ps_addRO)==1)]<-1
    
    ## true potential outcomes and causal effect ##
    Y0<-(.2*u2+.1*u2^2+1)
    Y1<-3/(1+exp(-(u2-1)*10))+(0.25*u1)-(.1*u1*u2)+.5
    ce_true<-Y1-Y0
    
    ## observed outcome data ##
    Yobs<-rep(NA,N)
    Yobs[which(x==0)]<-Y0[which(x==0)]
    Yobs[which(x==1)]<-Y1[which(x==1)]
    
    ## create trimmed and untrimmed datasets ##
    datall_mis<-data.frame(Yobs,x,ps_mis,u1,u2)
    dattr_mis<-datall_mis[which(RO_mis==1),]
    
    ######################
    ## 1. untrimmed G&R ##
    ######################
    
    ce_ugr_mis<-gr(Y=datall_mis$Yobs,trt=datall_mis$x,ps=datall_mis$ps_mis,X=cbind(datall_mis$u1,datall_mis$u2),M=500,qps=quantile(datall_mis$ps_mis,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
    keepsim_mis<-ce_ugr_mis[[1]]
    
    ## only keep this simulation and continue to run BART if the n=3 criteria in step 1 of Gutman are met ##
    if (keepsim_mis==1){
      ce_ugr_mis<-ce_ugr_mis[2:length(ce_ugr_mis)]
      
      ####################
      ## 2. trimmed G&R ##
      ####################
      
      ce_tgr_mis<-gr(Y=dattr_mis$Yobs,trt=dattr_mis$x,ps=dattr_mis$ps_mis,X=cbind(dattr_mis$u1,dattr_mis$u2),M=500,qps=quantile(dattr_mis$ps_mis,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
      ce_tgr_mis<-ce_tgr_mis[2:length(ce_tgr_mis)]
      
      #######################
      ## 3. untrimmed BART ##
      #######################
      
      testdat_mis<-datall_mis
      testdat_mis$x<-1-datall_mis$x
      ce_ubart_mis<-bartalone(xtr=datall_mis[,2:5],ytr=datall_mis[,1],xte=testdat_mis[,2:5])[c(1,4)]

      #####################
      ## 4. trimmed BART ##
      #####################
      
      testdat_mis<-dattr_mis
      testdat_mis$x<-1-dattr_mis$x
      ce_tbart_mis<-bartalone(xtr=dattr_mis[,2:5],ytr=dattr_mis[,1],xte=testdat_mis[,2:5])[c(1,4)]

      #################
      ## 5. BART+SPL ##
      #################
      
      ce_rn_mis<-bartspl(datall=datall_mis,RO=RO_mis)
      
      #################
      ## SAVE OUTPUT ##
      #################
      
      save_true<-c(save_true,list(ce_true))
      save_tgr_mis<-c(save_tgr_mis,list(ce_tgr_mis))
      save_ugr_mis<-c(save_ugr_mis,list(ce_ugr_mis))
      save_tbart_mis<-c(save_tbart_mis,list(ce_tbart_mis))
      save_ubart_mis<-c(save_ubart_mis,list(ce_ubart_mis))
      save_rn_mis<-c(save_rn_mis,list(ce_rn_mis))
      
      simscomplete<-simscomplete+1
    }
  }
}

save(save_true,save_ugr_mis,save_tgr_mis,save_ubart_mis,save_tbart_mis,save_rn_mis,
     file=paste(wd,'results/bartspline_e',extr,'_s',simnum,'.RData',sep=''))
