#################################################
## THIS CODE IMPLEMENTS THE SIMULATIONS IN 3.2 ##
#################################################

## read in command line arguments ##
args<-commandArgs(TRUE)
wd<-args[1] #working directory#
simnum<-as.numeric(args[2]) #simulation number#
reps<-as.numeric(args[3]) #number of reps to be implemented#
N<-as.numeric(args[4]) #total N#
ncovs<-as.numeric(args[5]) #number of additional "potential confounders"#

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
source('functions/aceBB.R')
source('functions/gr.R')
source('functions/pw_overlap.R')

set.seed(simnum)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

save_true<-list()
save_gr_mis<-list()
save_rn_mis<-list()
save_rn_pso_mis<-list()

sims_complete<-0
while (sims_complete<=(reps-1)){
  
  ## construct 10 true confounders ##
  covs<-matrix(c(rbinom(n1,size=1,prob=.45),rbinom(n0,size=1,prob=.4),rnorm(n=n1,mean=2,sd=2),rnorm(n=n0,mean=1.3)),nrow=N,byrow = F)
  for (i in 1:4){
    covs<-cbind(covs,c(rbinom(n1,size=1,prob=.45),rbinom(n0,size=1,prob=.4)),c(rnorm(n=n1,mean=2,sd=2),rnorm(n=n0,mean=1.3)))
  }
  
  if (ncovs>0){
    ## make additional ncov covariates ##
    covs<-cbind(covs,matrix(rnorm(ncovs*N),nrow=N,ncol=ncovs))
  }
  
  ## estimate misspecified PS ##
  emod<-glm(x~covs,family=binomial())
  ps_mis<-emod$fitted
  
  ## find RO and RN ##
  RO_mis<-pw_overlap(ps=ps_mis,E=x,a=.1,b=7)

  if (sum(1-RO_mis)>0){
    
    ## true potential outcomes and causal effects ##
    Y0<-rowSums(.5*covs[,c(1,3,5,7,9)])+(15/(1+exp(-(covs[,2]*8-1))))+rowSums(covs[,c(4,6,8,10)])-5
    Y1<-rowSums(covs[,c(1,3,5,7,9)]-.5*covs[,c(2,4,6,8,10)])
    ce_true<-Y1-Y0
    
    ## observed outcome data ##
    Yobs<-rep(NA,N)
    Yobs[which(x==0)]<-Y0[which(x==0)]
    Yobs[which(x==1)]<-Y1[which(x==1)]
    
    ## create a dataset with PS only and a dataset with PS and all confounders ##
    datall_mis<-data.frame(Yobs,x,ps_mis,covs)
    datpso_mis<-data.frame(Yobs,x,ps_mis)
    
    ############
    ## 1. G&R ##
    ############
    
    ce_gr_mis<-gr(Y=datall_mis$Yobs,trt=datall_mis$x,ps=datall_mis$ps_mis,X=covs,M=500,qps=quantile(datall_mis$ps_mis,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
    keepsim_mis<-ce_gr_mis[[1]]
    
    ## only keep this simulation and continue to run BART if the n=3 criteria in step 1 of Gutman are met ##
    if (keepsim_mis==1){
      ce_gr_mis<-ce_gr_mis[2:length(ce_gr_mis)]

      #####################################
      ## 2. BART+SPL with all covariates ##
      #####################################
      
      ce_rn_mis<-try(bartspl(datall=datall_mis,RO=RO_mis))
      
      ##############################
      ## 3. BART+SPL with only PS ##
      ##############################
      
      ce_rn_pso_mis<-try(bartspl(datall=datpso_mis,RO=RO_mis)) 
      
      if ("try-error" %in% class(ce_rn_mis) | "try-error" %in% class(ce_rn_pso_mis)){}
      ## store results from each method ##
      else{
        save_true<-c(save_true,list(ce_true))
        save_gr_mis<-c(save_gr_mis,list(ce_gr_mis))
        save_rn_mis<-c(save_rn_mis,list(ce_rn_mis))
        save_rn_pso_mis<-c(save_rn_pso_mis,list(ce_rn_pso_mis))
        
        sims_complete<-sims_complete+1
      }
    }
  }
  
}
save(save_true,save_gr_mis,save_rn_mis,save_rn_pso_mis,
     file=paste(wd,'results/ncovs',ncovs,'_s',simnum,'.RData',sep=''))
  