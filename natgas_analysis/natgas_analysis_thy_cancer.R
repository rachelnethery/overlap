################################################################
## thyroid cancer and natural gas compressor station analysis ##
## Rachel Nethery                                             ##
## last edit 5/9/18                                           ##
################################################################

## load packages and read in necessary functions ##
library(Hmisc)
library(MASS)
library(stats)
library(splines)
library(MCMCpack)
library(BayesTree)
library(dbarts)
library(rgdal)
library(xtable)
source('iw_overlap.R')
source('rn_2spl.R')
source('aceBB.R')
source('bartalone_noBB.R')

## read in ACS confounder data ##
acs<-read.csv('data/acs_2014endyr.csv',skip=1)

## read in thyroid cancer outcome data ##
hout<-read.csv('data/thyroid.csv')

## thyroid cancer dataset contains rows for entire states as well as each county, remove the rows that correspond to entire states ##
hout<-hout[which(hout$FIPS>1000),]

## read in health outcome data ##
hcov<-read.csv('data/socialex_health_2014.csv',skip=1)

## read in compressor station location data ##
comp<-read.csv('data/Natural_Gas_Compressor_Stations.csv',stringsAsFactors = F,na.strings = "NOT APPLICABLE")

## remove compressor stations that don't have location info ##
comp<-comp[which(is.na(comp$COUNTYFIPS)==F),]

## find out how many compressors have peak dates before 1/1/2012 ##
comp$peak_date<-unlist(lapply(strsplit(as.character(comp$OP_DATE_PE),'T'),function(x) x[1]))
comp$peak_date<-as.Date(comp$peak_date,format='%Y-%m-%d')
sum(comp$peak_date<=as.Date('2012-01-01',format='%Y-%m-%d'),na.rm=T)

## read in US county shapefile ##
shape <- readOGR(dsn = "data/cb_2014_us_county_5m", layer = "cb_2014_us_county_5m")

## get county centroids from the shapefile ##
countycent<-data.frame(as.numeric(paste(as.character(shape@data$STATEFP),as.character(shape@data$COUNTYFP),sep='')),coordinates(shape))
names(countycent)[1]<-c('FIPS')

## subset datasets to include a common set of counties ##
acs<-acs[which((acs$Geo_FIPS %in% hout$FIPS)==1 & (acs$Geo_FIPS %in% hcov$Geo_FIPS)==1),]
hout<-hout[which((hout$FIPS  %in% acs$Geo_FIPS)==1),]
hcov<-hcov[which((hcov$Geo_FIPS  %in% acs$Geo_FIPS)==1),]
countycent<-countycent[which((countycent$FIPS %in% acs$Geo_FIPS)==1),]

## align ordering across datasets ##
acs<-acs[order(acs$Geo_FIPS),]
hcov<-hcov[order(hcov$Geo_FIPS),]
hout<-hout[order(hout$FIPS),]
countycent<-countycent[order(countycent$FIPS),]

## subset to only the midwestern region of the US-- counties with centroid longitude between -90 and -110 ##
acs<-acs[which(countycent$X1<(-90) & countycent$X1>(-110)),]
hcov<-hcov[which(countycent$X1<(-90) & countycent$X1>(-110)),]
hout<-hout[which(countycent$X1<(-90) & countycent$X1>(-110)),]
countycent<-countycent[which(countycent$X1<(-90) & countycent$X1>(-110)),]

## create an empty table to store thyroid cancer analysis results ##
results_tab<-NULL

###################################################
## outcome is 2014 thyroid cancer mortality rate ##
###################################################

set.seed(4)

hout$mort2014<-unlist(lapply(strsplit(as.character(hout$Mortality.Rate..2014.),' '),function(x) as.numeric(x[1])))

## create dataset with outcome and confounders ##
test1<-data.frame(## FIPS codes ##
                  hout$FIPS,
                  ## log-transformed outcome ##
                  log(hout$mort2014),
                  ## confounders ##
                  ## PCP rate/100,000 ##
                  hcov$SE_NV003_001,
                  ## % of people without insurance ##
                  hcov$SE_T006_003,
                  ## % diabetics ##
                  hcov$SE_T009_001,
                  ## % current smokers ##
                  hcov$SE_T011_001,
                  ## % of people with limited access to healthy foods ##
                  hcov$SE_T012_001,
                  ## % obese ##
                  hcov$SE_T012_003,
                  ## food environment index ##
                  hcov$SE_T013_001,
                  ## pop dens ##
                  acs$SE_T002_002,
                  ## % male ##
                  acs$PCT_SE_T004_002,
                  ## % <55 ##
                  acs$PCT_SE_T007_002+acs$PCT_SE_T007_003+acs$PCT_SE_T007_004+acs$PCT_SE_T007_005+acs$PCT_SE_T007_006+acs$PCT_SE_T007_007+acs$PCT_SE_T007_008+acs$PCT_SE_T007_009,
                  ## % white ##
                  acs$PCT_SE_T013_002,
                  ## average household size ##
                  acs$SE_T021_001,
                  ## % bachelors degree or more ##
                  acs$PCT_SE_T150_005,
                  ## % unemployed ##
                  acs$PCT_SE_T037_003,
                  ## median household income ##
                  acs$SE_T057_001,
                  ## gini index ##
                  acs$SE_T157_001,
                  ## % owner occupied housing units ##
                  acs$PCT_SE_T094_002,
                  ## median rent as proportion of income ##
                  acs$SE_T105_001,
                  ## avg commute time ##
                  acs$SE_T147_001

)

## names of columns in the analysis dataset ##
names(test1)<-c('fips','out','pcp','insure','diabetes','smoke','healthfood','obese','foodenv','popdens','male','lt55','white','hhsize',
                'edu','unemploy','hhinc','gini','tenure','rentpinc','commute')

## add an exposure indicator for each county in the dataset, 1 if it contains a comp station, 0 otherwise ##
test1$expose<-as.numeric(test1$fips %in% comp$COUNTYFIPS)

## only keep complete cases (those with no missingness) ##
test2<-test1[complete.cases(test1),]

## fit the propensity score model using BART probit ##
bartfit<-bart(x.train = test2[,3:(ncol(test1)-1)],y.train = test2[,ncol(test1)])
bartps<-colMeans(pnorm(bartfit$yhat.train))

## exposure stratified histograms of the PS ##
hist(bartps[which(test2$expose==0)],col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,1),breaks=15,xlab="",ylab='',main='')
hist(bartps[which(test2$expose==1)],col=rgb(0.8,0.8,0.8,0.5),add=T,xlab="",ylab='',breaks=15)

## add the ps to the dataset ##
test2$ps<-bartps

## identify the RO and RN ##
RO<-iw_overlap(ps=bartps,E=test2$expose,a=.1*(max(bartps)-min(bartps)),b=10)

## overlap stratified histograms of the PS ##
hist(bartps[which(RO==1)],col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,1),breaks=10,xlab="",ylab='',main='')
hist(bartps[which(RO==0)],col=rgb(0.8,0.8,0.8,0.5),add=T,xlab="",ylab='',breaks=10)

## save histogram ##
pdf('C:/Users/Rachel/Documents/Harvard/Methods_Development/overlap/bart+spl/natgas/natgas_hist.pdf',width=5,height=5)
hist(bartps[which(test2$expose==0)],col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,.8),breaks=15,ylab='',main='',xlab="Propensity Score")
hist(bartps[which(test2$expose==1)],col=rgb(0.8,0.8,0.8,0.5),add=T,xlab="",ylab='',breaks=15)
abline(v=min(bartps[which(RO==1)]),lwd=2)
abline(v=max(bartps[which(RO==1)]),lwd=2)
legend('topright',c(expression(paste(hat(xi)[i]," | ",E[i]==1,sep='')),
                    expression(paste(hat(xi)[i]," | ",E[i]==0,sep=''))),pch=c(15,15),col=c(rgb(0.8,0.8,0.8,0.5),rgb(0.1,0.1,0.1,0.5)))

dev.off()

## estimate causal effect of exposure on thyroid cancer using trimmed BART ##
test3<-test2
test3$expose<-1-test3$expose
bartfit<-bartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
                   ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

mean(bartfit[[4]])
quantile(bartfit[[4]],c(.025,.975))

## estimate causal effect of exposure on thyroid cancer using BART+SPL ##
rnfit<-rn_2spl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],RO=RO)
mean(rnfit[[1]])
quantile(rnfit[[1]],c(.025,.975))

## add results to the table ##
results_tab<-rbind(results_tab,
                   c(sum(RO==0)/nrow(test2),mean(rnfit[[1]]),quantile(rnfit[[1]],c(.025,.975))),
                   c(sum(RO==0)/nrow(test2),mean(bartfit[[4]]),quantile(bartfit[[4]],c(.025,.975))))

##################################################################
## outcome is change in thyroid cancer mortality rate 1980-2014 ##
##################################################################

set.seed(4)

hout$changemort<-unlist(lapply(strsplit(as.character(hout$X..Change.in.Mortality.Rate..1980.2014),' '),function(x) as.numeric(x[1])))

## create dataset with outcome and confounders ##
test1<-data.frame(## FIPS codes ##
                  hout$FIPS,
                  ## outcome premature death rate ##
                  hout$changemort,
                  ## confounders ##
                  ## PCP rate/100,000 ##
                  hcov$SE_NV003_001,
                  ## % of people without insurance ##
                  hcov$SE_T006_003,
                  ## % diabetics ##
                  hcov$SE_T009_001,
                  ## % current smokers ##
                  hcov$SE_T011_001,
                  ## % of people with limited access to healthy foods ##
                  hcov$SE_T012_001,
                  ## % obese ##
                  hcov$SE_T012_003,
                  ## food environment index ##
                  hcov$SE_T013_001,
                  ## pop dens ##
                  acs$SE_T002_002,
                  ## % male ##
                  acs$PCT_SE_T004_002,
                  ## % <55 ##
                  acs$PCT_SE_T007_002+acs$PCT_SE_T007_003+acs$PCT_SE_T007_004+acs$PCT_SE_T007_005+acs$PCT_SE_T007_006+acs$PCT_SE_T007_007+acs$PCT_SE_T007_008+acs$PCT_SE_T007_009,
                  ## % white ##
                  acs$PCT_SE_T013_002,
                  ## average household size ##
                  acs$SE_T021_001,
                  ## % bachelors degree or more ##
                  acs$PCT_SE_T150_005,
                  ## % unemployed ##
                  acs$PCT_SE_T037_003,
                  ## median household income ##
                  acs$SE_T057_001,
                  ## gini index ##
                  acs$SE_T157_001,
                  ## % owner occupied housing units ##
                  acs$PCT_SE_T094_002,
                  ## median rent as proportion of income ##
                  acs$SE_T105_001,
                  ## avg commute time ##
                  acs$SE_T147_001)

## names of columns in the analysis dataset ##
names(test1)<-c('fips','out','pcp','insure','diabetes','smoke','healthfood','obese','foodenv','popdens','male','lt55','white','hhsize',
                'edu','unemploy','hhinc','gini','tenure','rentpinc','commute')

## add an exposure indicator for each county in the dataset, 1 if it contains a comp station, 0 otherwise ##
test1$expose<-as.numeric(test1$fips %in% comp$COUNTYFIPS)

## only keep complete cases (those with no missingness) ##
test2<-test1[complete.cases(test1),]

## fit the propensity score model using BART probit ##
bartfit<-bart(x.train = test2[,3:(ncol(test1)-1)],y.train = test2[,ncol(test1)])
bartps<-colMeans(pnorm(bartfit$yhat.train))

## exposure stratified histograms of the PS ##
hist(bartps[which(test2$expose==0)],col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,1),breaks=15,xlab="",ylab='',main='')
hist(bartps[which(test2$expose==1)],col=rgb(0.8,0.8,0.8,0.5),add=T,xlab="",ylab='',breaks=15)

## add the ps to the dataset ##
test2$ps<-bartps

## identify the RO and RN ##
RO<-iw_overlap(ps=bartps,E=test2$expose,a=.1*(max(bartps)-min(bartps)),b=10)

## overlap stratified histograms of the PS ##
hist(bartps[which(RO==1)],col=rgb(0.1,0.1,0.1,0.5),xlim=c(0,1),breaks=10,xlab="",ylab='',main='')
hist(bartps[which(RO==0)],col=rgb(0.8,0.8,0.8,0.5),add=T,xlab="",ylab='',breaks=10)

## estimate causal effect of exposure on thyroid cancer using trimmed BART ##
test3<-test2
test3$expose<-1-test3$expose
bartfit<-bartalone(xtr=test2[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],
                   ytr=test2[which(RO==1),2],xte=test3[which(RO==1),c(ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))])

mean(bartfit[[4]])
quantile(bartfit[[4]],c(.025,.975))

## estimate causal effect of exposure on thyroid cancer using BART+SPL ##
rnfit<-rn_2spl(datall=test2[,c(2,ncol(test2)-1,ncol(test2),3:(ncol(test2)-2))],RO=RO)
mean(rnfit[[1]])
quantile(rnfit[[1]],c(.025,.975))

## add results to the table ##
results_tab<-rbind(results_tab,
                   c(sum(RO==0)/nrow(test2),mean(rnfit[[1]]),quantile(rnfit[[1]],c(.025,.975))),
                   c(sum(RO==0)/nrow(test2),mean(bartfit[[4]]),quantile(bartfit[[4]],c(.025,.975))))

## print thyroid cancer results table ##
results_tab<-results_tab[,2:4]
colnames(results_tab)<-NULL
xtable(results_tab,digits=3)
