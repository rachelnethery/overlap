##########################################################
## THIS CODE IMPLEMENTS SIMULATION 2A IN THE SUPPLEMENT ##
##########################################################

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
source('functions/gr_binary.R')
source('functions/bartspl_bin.R')
source('functions/bartalone_binary.R')
source('functions/aceBB.R')
source('functions/pw_overlap.R')

set.seed(simnum)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

save_true<-list()
save_ugr_cor<-list()
save_ubart_cor<-list()
save_rn_cor<-list()

simscomplete<-0
while (simscomplete<=(reps-1)){
  
  ## generate confounders (1 binary, 1 continuous) ##
  u1<-c(rbinom(n=n1,size=1,prob=.5),rbinom(n=n0,size=1,prob=.4))
  u2<-c(rnorm(n=n1,mean=2+extr,sd=1.25+extr*.1),rnorm(n=n0,mean=1,sd=1))

  ## construct true PS ##
  ps_cor<-(dnorm(u2, mean = 2+extr, sd = 1.25+extr*.1, log = FALSE)*dbinom(u1,size=1,prob=.5,log=FALSE))/
    ((dnorm(u2, mean = 2+extr, sd =1.25+extr*.1, log = FALSE)*dbinom(u1,size=1,prob=.5,log=FALSE))+
       (dnorm(u2,mean=1,sd=1,log=FALSE)*dbinom(u1,size=1,prob=.4,log=FALSE)))
  
  ## find RO and RN ## 
  RO_cor<-pw_overlap(ps=ps_cor,E=x,a=a,b=b)
  order_ps_cor<-ps_cor[order(ps_cor)]
  order_RO_cor<-RO_cor[order(ps_cor)]
  ll_cor<-max(which(order_RO_cor==1))
  
  ## skip over datasets that have no non-overlap in right tail or too much non-overlap in the left tail (RNs with >10 PS's) ##
  if (ll_cor<N & sum(1-order_RO_cor[1:ll_cor])<=10){
    
    ## ignore small non-overlap in the left tail (add all left tail to the RO) ##
    temp1<-order_ps_cor[1:ll_cor]
    ps_addRO<-temp1[which(order_RO_cor[1:ll_cor]==0)]
    RO_cor[which((ps_cor %in% ps_addRO)==1)]<-1
    
    ## true potential outcomes and causal effect ##
    Y0_s<-(.2)*(u2)^3+.25*u1
    Y1_s<-(exp(.25*u2))+.5*u1*u2
    Y0_p<-exp(Y0_s)/(1+exp(Y0_s))
    Y1_p<-exp(Y1_s)/(1+exp(Y1_s))
    ce_true<-Y1_p-Y0_p
    
    ## observed outcome data ##
    Yobs<-rep(NA,N)
    Yobs[which(x==0)]<-rbinom(n=n0,size=1,prob=Y0_p[which(x==0)])
    Yobs[which(x==1)]<-rbinom(n=n1,size=1,prob=Y1_p[which(x==1)])
    
    ## create dataset ##
    datall_cor<-data.frame(Yobs,x,ps_cor,u1,u2)
    
    ######################
    ## 1. untrimmed G&R ##
    ######################

    ce_ugr_cor<-try(gr(Y=datall_cor$Yobs,trt=datall_cor$x,ps=datall_cor$ps_cor,X=cbind(datall_cor$u1,datall_cor$u2),M=500,qps=quantile(datall_cor$ps_cor,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T)))
    if ("try-error" %in% class(ce_ugr_cor)){}
    else{
      
      keepsim_cor<-ce_ugr_cor[[1]]
      
      ## only keep this simulation and continue to run BART if the n=3 criteria in step 1 of Gutman are met ##
      if (keepsim_cor==1){
        ce_ugr_cor<-ce_ugr_cor[2:length(ce_ugr_cor)]
        
        #######################
        ## 2. untrimmed BART ##
        #######################
        
        testdat_cor<-datall_cor
        testdat_cor$x<-1-datall_cor$x
        ce_ubart_cor<-bartalone(xtr=datall_cor[,2:5],ytr=datall_cor[,1],xte=testdat_cor[,2:5])[c(1,4)]
        
        #################
        ## 3. BART+SPL ##
        #################
        
        ce_rn_cor<-bartspl_bin(datall=datall_cor,RO=RO_cor)
        
        #################
        ## SAVE OUTPUT ##
        #################
        
        save_true<-c(save_true,list(ce_true))
        save_ugr_cor<-c(save_ugr_cor,list(ce_ugr_cor))
        save_ubart_cor<-c(save_ubart_cor,list(ce_ubart_cor))
        save_rn_cor<-c(save_rn_cor,list(ce_rn_cor))
        
        simscomplete<-simscomplete+1
      }
    }
  }
}

save(save_true,save_ugr_cor,save_ubart_cor,save_rn_cor,
     file=paste(wd,'results/bartspline_cor_e',extr,'_s',simnum,'.RData',sep=''))