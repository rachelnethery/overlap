#####################################################################
## THIS CODE IMPLEMENTS THE SIMULATION IN SECTION 3.3 OF THE PAPER ##
#####################################################################

## read in command line arguments ##
args<-commandArgs(TRUE)
wd<-args[1] #working directory #
simnum<-as.numeric(args[2]) #simulation number#
reps<-as.numeric(args[3]) #number of reps to be implemented#
N<-as.numeric(args[4]) #total N#
ncovs<-as.numeric(args[5]) #number of irrelevant predictors to add ##

## run simulation ##
setwd(wd)
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
source('rn_2spl.R')
source('aceBB.R')
source('gr.R')
source('iw_overlap.R')

set.seed(simnum)

n1<-N/2
n0<-N/2

x<-c(rep(1,n1),rep(0,n0))

save_true<-list()
save_gr<-list()
save_rn<-list()
save_rn_pso<-list()

sims_complete<-0
while (sims_complete<=(reps-1)){
  
  ## make 10 true confounders ##
  covs<-matrix(c(rbinom(n1,size=1,prob=.45),rbinom(n0,size=1,prob=.4),rnorm(n=n1,mean=1.25,sd=3),runif(n=n0,min=-2,max=2.5)),nrow=N,byrow = F)
  for (i in 1:4){
    covs<-cbind(covs,c(rbinom(n1,size=1,prob=.5),rbinom(n0,size=1,prob=.4)),c(rnorm(n=n1,mean=1.25,sd=3),runif(n=n0,min=-2,max=2.5)))
  }
  
  if (ncovs>0){
    ## make additional ncov covariates ##
    covs<-cbind(covs,matrix(rnorm(ncovs*N),nrow=N,ncol=ncovs))
  }
  
  ## estimate PS ##
  emod<-glm(x~covs,family=binomial())
  ps<-emod$fitted
  
  ## find RO and RN ##
  RO<-iw_overlap(ps=ps,E=x,a=.15,b=5)

  if (sum(RO==0)>0){
    ## true potential outcomes and causal effects ##
    Y0<-rowSums(.2*covs[,c(1,3,5,7,9)])+(1/(1+exp(-(covs[,2]*8-1))))+rowSums(covs[,c(4,6,8,10)])+5
    Y1<-rowSums(.2*covs[,c(1,3,5,7,9)]-.5*covs[,c(2,4,6,8,10)])-5
    ce_true<-Y1-Y0
    
    ## observed outcome data ##
    Yobs<-rep(NA,N)
    Yobs[which(x==0)]<-Y0[which(x==0)]
    Yobs[which(x==1)]<-Y1[which(x==1)]
    
    datall<-data.frame(Yobs,x,ps,covs)
    datpso<-data.frame(Yobs,x,ps)
    
    ## fit G&R ##
    ce_gr<-gr(Y=datall$Yobs,trt=datall$x,ps=datall$ps,X=covs,M=500,qps=quantile(datall$ps,probs=c(0,.3,.4,.5,.6,.7,1),na.rm=T))
    keepsim<-ce_gr[[1]]
    
    ## only keep this simulation and continue to run BART if the n=3 criteria in step 1 of Gutman are met ##
    if (keepsim==1){
      ce_gr<-ce_gr[2:length(ce_gr)]
      
      ## fit BART+SPL with all covariates ##
      ce_rn<-try(rn_2spl(datall=datall,RO=RO))
      
      ## fit BART+SPL with only PS ##
      ce_rn_pso<-try(rn_2spl(datall=datpso,RO=RO))  
      
      if ("try-error" %in% class(ce_rn) | "try-error" %in% class(ce_rn_pso)){}
      ## store results from each method ##
      else{
        save_true<-c(save_true,list(ce_true))
        save_gr<-c(save_gr,list(ce_gr))
        save_rn<-c(save_rn,list(ce_rn))
        save_rn_pso<-c(save_rn_pso,list(ce_rn_pso))
        
        sims_complete<-sims_complete+1
      }
    }
  }
}

## save all results ##
save(save_true,save_gr,save_rn,save_rn_pso,file=paste(wd,'results/ncovs',ncovs,'_s',simnum,'.RData',sep=''))
  