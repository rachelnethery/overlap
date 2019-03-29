#################################################
## THIS CODE IMPLEMENTS THE SIMULATIONS IN 3.3 ##
#################################################

## read in command line arguments ##
args<-commandArgs(TRUE)
wd<-args[1] #working directory#
simnum<-as.numeric(args[2]) #simulation number#
reps<-as.numeric(args[3]) #number of reps to be implemented#
N<-as.numeric(args[4]) #total N#
kk<-args[5] #string "easy" for simulation 3.3A and string "hard" for simulation 3.3B#

## run simulation ##
setwd(wd)
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
source('functions/bartspl.R')
source('functions/bartalone.R')
source('functions/aceBB.R')
source('functions/pw_overlap.R')

## 3 different v,w settings ##
mu<-c(1.4,0.75,0)
sig<-c(1.4,1.2,1)
seed<-c(25,25,21)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

for (gg in 1:length(mu)){
  set.seed(seed[gg])
  
  ## generate confounder ##
  u<-c(rnorm(n=n1,mean=2.5,sd=2),rnorm(n=n0,mean=mu[gg],sd=sig[gg]))
  
  set.seed(simnum)
  
  ## 3 different RO/RN definitions ##
  for (hh in 1:3){
    if (hh==1) RO<-pw_overlap(ps=u,E=x,a=.05*(max(u)-min(u)),b=10)
    else if (hh==2) RO<-pw_overlap(ps=u,E=x,a=.1*(max(u)-min(u)),b=10)
    else if (hh==3) RO<-pw_overlap(ps=u,E=x,a=.15*(max(u)-min(u)),b=3)
    
    save_true<-list()
    save_rn<-list()
    save_bart<-list()
    
    sims_complete<-0
    while (sims_complete<=(reps-1)){
      
      ## true potential outcomes and causal effects ##
      if (kk=='easy'){
        Y0<-((u+(u^2/factorial(2)))/20)+1.5+rnorm(N,sd=.25)
        Y1<-(1/(1+exp(-(u-1))))+rnorm(N,sd=.25)
      } else if (kk=='hard'){
        Y0<-((u-(u^2/factorial(2))+(u^3/factorial(3)))/20)+1.5+rnorm(N,sd=.25)
        Y1<-(1/(1+exp(-(u-1))))+rnorm(N,sd=.25)
      }
      
      ce_true<-Y1-Y0
      
      ## observed outcome data ##
      Yobs<-rep(NA,N)
      Yobs[which(x==0)]<-Y0[which(x==0)]
      Yobs[which(x==1)]<-Y1[which(x==1)]
      
      ## construct dataset ##
      datall<-data.frame(Yobs,x,u)
      
      ## fit BART+SPL ##
      ce_rn<-try(bartspl(datall=datall,RO=RO))
      
      if ("try-error" %in% class(ce_rn)){}
      else{
        ## fit BART ##
        testdat<-datall
        testdat$x<-1-datall$x
        ce_bart<-bartalone(xtr=datall[,2:ncol(datall)],ytr=datall[,1],xte=testdat[,2:ncol(testdat)])
        
        ## save output ##
        save_true<-c(save_true,list(ce_true))
        save_rn<-c(save_rn,list(ce_rn))
        save_bart<-c(save_bart,list(ce_bart))
        
        sims_complete<-sims_complete+1
      }
    }
    save(save_true,save_rn,save_bart,file=paste(wd,'results/asim',kk,'_mu',mu[gg],'_hh',hh,'_s',simnum,'.RData',sep=''))
  }
}
  
  