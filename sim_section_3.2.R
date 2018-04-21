#####################################################################
## THIS CODE IMPLEMENTS THE SIMULATION IN SECTION 3.2 OF THE PAPER ##
#####################################################################

## read in command line arguments ##
args<-commandArgs(TRUE)
wd<-args[1] #working directory #
simnum<-as.numeric(args[2]) #simulation number#
reps<-as.numeric(args[3]) #number of reps to be implemented#
N<-as.numeric(args[4]) #total N#
extr<-as.numeric(args[5]) #parameter controlling degree of non-overlap#
a<-as.numeric(args[6]) #a parameter in interval-wise non-overlap#
b<-as.numeric(args[7]) #b parameter in interval-wise non-overlap#

## run simulation ##
setwd(wd)
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
library(lme4)
source('gr_binary.R')
source('rn_2spl_binary.R')
source('bartalone_binary.R')
source('aceBB.R')
source('iw_overlap.R')

set.seed(simnum)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

save_true<-list()
save_tgr<-list()
save_ugr<-list()
save_tbart<-list()
save_ubart<-list()
save_rn<-list()

simscomplete<-0
while (simscomplete<=(reps-1)){
  
  ## generate confounders (1 binary, 1 continuous) ##
  u1<-c(rbinom(n=n1,size=1,prob=.5),rbinom(n=n0,size=1,prob=.4))
  u2<-c(rnorm(n=n1,mean=1.25,sd=1),runif(n=n0,min=-1.15-extr,max=1.15+extr))
  
  ## construct propensity score ##
  emod<-glm(x~u1+u2,family = binomial())
  ps<-emod$fitted
  
  ## find RO and RN ## 
  RO<-iw_overlap(ps=ps,E=x,a=a,b=b)
  order_ps<-ps[order(ps)]
  order_RO<-RO[order(ps)]
  ll<-max(which(order_RO==1))
  
  ## skip datasets that have too much non-overlap in the left tail (RNs with >10 PS's) ##
  if (ll<N & sum(1-order_RO[1:ll])<=10){
    
    ## ignore small overlap in the left tail (add all left tail to the RO) ##
    temp1<-order_ps[1:ll]
    ps_addRO<-temp1[which(order_RO[1:ll]==0)]
    RO[which((ps %in% ps_addRO)==1)]<-1
    
    ## true potential outcomes and causal effect ##
    Y0_s<-.85*(u2-1)+(u2-1)^2
    Y1_s<-(1/(1+exp(-((u2*8)-1))))+(0.25*u1)-2
    Y0_p<-exp(Y0_s)/(1+exp(Y0_s))
    Y1_p<-exp(Y1_s)/(1+exp(Y1_s))
    ce_true<-Y1_p-Y0_p
    
    ## observed outcome data ##
    Yobs<-rep(NA,N)
    Yobs[which(x==0)]<-rbinom(n=n0,size=1,prob=Y0_p[which(x==0)])
    Yobs[which(x==1)]<-rbinom(n=n1,size=1,prob=Y1_p[which(x==1)])
    
    datall<-data.frame(Yobs,x,ps,u1,u2)
    dattr<-datall[which(RO==1),]
    
    ## fit untrimmed G&R ##
    ce_ugr<-gr(Y=datall$Yobs,trt=datall$x,ps=datall$ps,X=cbind(datall$u1,datall$u2),M=500,qps=quantile(datall$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
    keepsim<-ce_ugr[[1]]
    
    ## only keep this simulation and continue to run BART if the n=3 criteria in step 1 of Gutman are met ##
    if (keepsim==1){
      ce_ugr<-ce_ugr[2:length(ce_ugr)]
      
      ## fit trimmed G&R ##
      ce_tgr<-gr(Y=dattr$Yobs,trt=dattr$x,ps=dattr$ps,X=cbind(dattr$u1,dattr$u2),M=500,qps=quantile(dattr$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
      ce_tgr<-ce_tgr[2:length(ce_tgr)]
      
      ## fit untrimmed BART ##
      testdat<-datall
      testdat$x<-1-datall$x
      ce_ubart<-bartalone(xtr=datall[,2:5],ytr=datall[,1],xte=testdat[,2:5])[c(1,4)]
      
      ## fit trimmed BART ##
      testdat<-dattr
      testdat$x<-1-dattr$x
      ce_tbart<-bartalone(xtr=dattr[,2:5],ytr=dattr[,1],xte=testdat[,2:5])[c(1,4)]
      
      ## fit BART+SPL ##
      ce_rn<-rn_2spl(datall=datall,RO=RO)
      
      ## store results from each method ##
      save_true<-c(save_true,list(ce_true))
      save_tgr<-c(save_tgr,list(ce_tgr))
      save_ugr<-c(save_ugr,list(ce_ugr))
      save_tbart<-c(save_tbart,list(ce_tbart))
      save_ubart<-c(save_ubart,list(ce_ubart))
      save_rn<-c(save_rn,list(ce_rn))
    }
    
    simscomplete<-simscomplete+keepsim
    
  }
}

## save all results ##
save(save_true,save_ugr,save_tgr,save_ubart,save_tbart,save_rn,file=paste(wd,'results/bartspline_e',extr,'_s',simnum,'.RData',sep=''))
