# Clear workspace to get rid of old junk
rm(list=ls(all=TRUE))
##SET RANDOM SEED
set.seed(123)

library (RCurl)
library(boot)
library(reshape2)
library(RColorBrewer)
library(matlib)
library(lubridate)

jagsExists = url.exists("http://sourceforge.net/projects/mcmc-jags/files/latest/download")
# Regular HTTP
if(jagsExists) {
  txt = getURL("http://sourceforge.net/projects/mcmc-jags/files/latest/download")}

#### Use jags package ####
#install.packages("rjags")
library(rjags)

#############################################
data.start<-1 #Determines which year to include; set to 1 to include full TS, to 61 to just 2008+
#############################################

#Import the data and select the right age group, get rid of columns with NAs

library(splines, quietly = TRUE)
set.seed(1)

setwd('C:/Users/DMW63/Desktop/')
##################################################
##################################################
#SECTION 1: IMPORTING AND FORMATTING TIME SERIES

in.data<-read.csv('./prelog_Brazil_processed_data.csv')

output_directory<-'./output'
dir.create(output_directory, recursive=TRUE, showWarnings = FALSE)
ds1a<-in.data[in.data$age_group==9 ,]
age_groups <- paste( unique(unlist(ds1a$age_group, use.names = FALSE)))
#filter if  column mean<5?
exclude<-c( which(colMeans(ds1a[,-c(1:3)])<5 )+3, which(is.na(colMeans(ds1a[,-c(1:3)])))+3)
ds1a<-ds1a[,-exclude ] #Position of variables to filter out


##Account for code-naming differences  "all_cause_pneu_name" gives the outcome variable (e.g. pneumonia); "noj_name" gives the name of the bronchitis/bronchiolitis
#Variable that is excluded from sensitivity analyses
#All_cause_name gives name of denominator used for  some of the analyses--can be population size, non-respiratory hospitalizations, tec 

all_cause_name <- 'ach_noj'
all_cause_pneu_name <- 'J12_18'
noj_name <- 'cJ20_J22'
date_name<-'date'

##Date variables
#When do you want the analysis to start?
data_start_date <- as.Date('2004-01-01')
data_intervention_date <- as.Date('2009-12-31') 

potential_covars<-names(ds1a)[4:ncol(ds1a)]

#When does the dataset end
data_end_date <- max(as.Date(ds1a[, date_name]))
pre_period <- c(data_start_date, data_intervention_date)  #define training period
post_period <- c(data_intervention_date + 1, data_end_date) #Define post-vaccine period
#Define evaluation period

eval_period <- c( (data_end_date %m-% months(24)), data_end_date) #take last 24m as evaluation period

ds <- ds1a
ds$date<-as.Date(ds$date)
ds <- ds[, colSums(is.na(ds)) == 0] #Delete columns with NAs
ds <- ds[match(data_start_date, ds$date):nrow(ds),]
#ds[ds == 0] <- 0.5
#ds[, 3:ncol(ds)] <- log(ds[, 3:ncol(ds)])


data_start <- match(data_start_date, ds$date)

time_points <- ds$date[data_start:nrow(ds)]

post.start.index<-which(time_points==pre_period[2]+1)

if(ncol(ds)>=4){
  covars.raw <- as.data.frame(ds[data_start:nrow(ds), 4:ncol(ds)])
  names(covars.raw)<-names(ds)[4:ncol(ds)]
}else{
  covars.raw<-as.data.frame(rep(1, nrow(ds)))
  names(covars.raw)<-'one'
}
if(ncol(ds)>4){
  covars.raw <- ds[data_start:nrow(ds), 4:ncol(ds)]
}else{
  covars.raw <- as.data.frame(ds[data_start:nrow(ds), 4])
  names(covars.raw)<-names(ds)[4]
}
month_i <- as.factor(as.numeric(format(ds[, date_name][data_start:nrow(ds)], '%m')))
spline <- setNames(as.data.frame(bs(1:nrow(covars.raw), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
year_2008 <- numeric(nrow(covars.raw))
year_2008[1:nrow(covars.raw) >= match(as.Date('2008-01-01'), ds[, date_name])] <- 1
data <- cbind.data.frame(year_2008, spline, month_i)
pandemic <- ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0))

variance.covars<-apply(as.data.frame(covars.raw[1:(post.start.index-1),]),2,var, na.rm=TRUE)
with.variation<-which(variance.covars>0)
no.variation<-which(variance.covars==0)
names.keep<-names(covars.raw)[with.variation]
if(length(no.variation)>0){covars.raw<-covin.data=ars.raw[,with.variation]} #eliminates covariates that have NO VARIABIITY IN PRE-PERIOD
covars.raw<-as.data.frame(covars.raw)
names(covars.raw)<-names.keep
covars<-covars.raw 
data.sel<-covars

covar<-data.sel
#Create Monthly dummies
m1 = rep(c(1,0,0,0,0,0,0,0,0,0,0,0), nrow(data.sel)/12) 
m2 = rep(c(0,1,0,0,0,0,0,0,0,0,0,0),nrow(data.sel)/12)
m3 = rep(c(0,0,1,0,0,0,0,0,0,0,0,0), nrow(data.sel)/12)
m4 = rep(c(0,0,0,1,0,0,0,0,0,0,0,0),nrow(data.sel)/12)
m5 = rep(c(0,0,0,0,1,0,0,0,0,0,0,0), nrow(data.sel)/12)
m6 = rep(c(0,0,0,0,0,1,0,0,0,0,0,0), nrow(data.sel)/12)
m7 = rep(c(0,0,0,0,0,0,1,0,0,0,0,0), nrow(data.sel)/12)
m8 = rep(c(0,0,0,0,0,0,0,1,0,0,0,0), nrow(data.sel)/12)
m9 = rep(c(0,0,0,0,0,0,0,0,1,0,0,0),nrow(data.sel)/12)
m10 =rep(c(0,0,0,0,0,0,0,0,0,1,0,0),nrow(data.sel)/12)
m11 = rep(c(0,0,0,0,0,0,0,0,0,0,1,0),nrow(data.sel)/12)
month.ind.matrix<-cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)

#Spline for post-PCV period only
post.index<-which(time_points==post_period[1]):nrow(data.sel)
index=1:length(time_points)
bs.test.post<-bs(index[post.index],  df=8, degree = 1)
knots<-as.vector(c(attr(bs.test.post,'Boundary.knots')[1], attr(bs.test.post,'knots')))
bs.test<-bs(index, knots=knots, degree=1)
#bs.test2<-rbind(matrix(0,nrow=(post.index[1]-1), ncol=ncol(bs.test) ), bs.test) 
bs.test2<-as.data.frame(bs.test)
names(bs.test2)<-paste0('spl', 1:9)


predictors<-as.matrix(covars)
predictors <-cbind(predictors, ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0)))
dimnames(predictors)[[2]][ncol(predictors)]<-'pandemic'

#outcome.pre<-outcome
#outcome.pre[post.start.index:length(outcome)]<-NA

#COMBINE MONTHLY DUMMIES AND COVARIATES AND PANDEMIC INTO SINGLE DATAFRAME
covar.matrix<-cbind.data.frame(month.ind.matrix,pandemic,bs.test2,predictors)
covar.lab<-dimnames(covar.matrix)[[2]]
log_offset=log(rep(200000, nrow(predictors)))

outcome=ds[data_start:nrow(ds),3]
time<-1:nrow(data.sel)
time_post<-(time[post.start.index:length(outcome)]-post.start.index+1)/100

data.fit<-cbind.data.frame(outcome, covar.matrix)
data.fit$outcome<-as.integer(data.fit$outcome)
one<-rep(1,times=nrow(data.fit))

month.mat<-data.fit[,c('m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11')]
spl.mat<-data.fit[,names(bs.test2)]
yr2008<-as.numeric(time_points>=as.Date("2008-01-01"))

spl.increment<-spl.mat[2,1]-spl.mat[1,1] #what is increment of change for spline 1
data.fit$t<-seq(from=0, length.out=nrow(spl.mat), by=spl.increment )

  
#########################################
##########################################
##JAGS
###################################
jcode_g5 <- "model
{
    #OUTCOME MODEL
  for (i in 1:N)     {
  log.lambda[i]<- (int+ log_ach_noj[i]+ beta[1]*month.mat[i,1]+beta[2]*month.mat[i,2]+beta[3]*month.mat[i,3]+beta[4]*month.mat[i,4]+beta[5]*month.mat[i,5]
                        +beta[6]*month.mat[i,6]+beta[7]*month.mat[i,7]
                        +beta[8]*month.mat[i,8]+beta[9]*month.mat[i,9]+beta[10]*month.mat[i,10]+beta[11]*month.mat[i,11] 
                        +beta[12]*pandemic[i]
                        +eps[1]*spl.mat[i,1] +eps[2]*spl.mat[i,2] +eps[3]*spl.mat[i,3]+eps[4]*spl.mat[i,4]
                        +eps[5]*spl.mat[i,5]+eps[6]*spl.mat[i,6]+eps[7]*spl.mat[i,7]+eps[8]*spl.mat[i,8]
                        +eps[9]*spl.mat[i,9])
  lambda[i]<-exp(log.lambda[i])
 log.lambda.novax[i]<- int+ log_ach_noj[i]+ beta[1]*month.mat[i,1]+beta[2]*month.mat[i,2]+beta[3]*month.mat[i,3]+beta[4]*month.mat[i,4]+beta[5]*month.mat[i,5]
                        +beta[6]*month.mat[i,6]+beta[7]*month.mat[i,7]
                        +beta[8]*month.mat[i,8]+beta[9]*month.mat[i,9]+beta[10]*month.mat[i,10]+beta[11]*month.mat[i,11] 
                        +beta[12]*pandemic[i]
                        +eps[1]*t.index[i]
  log.rr.vax[i]<- log.lambda[i] - log.lambda.novax[i]

  outcome[i] ~ dnegbin(p1[i], r1) # Observation variation
  p1[i]<-r1/(r1+lambda[i])

  }
  r1 ~ dunif(0, 50)


#prior for main model
int~dnorm(0,1e-4)
for(j in 1:12){  beta[j]~dnorm(0,1e-4)} #monthly covars and pandemic

#autoregressive structure for spline coefficients to impose smoothness
for(k in 1:2){eps[k]~dnorm(0, 1e-4)} #pre-vax trend
for(k in 3:9){eps[k]~dnorm(eps[k-1],tau.spl1)}  

    #precision of Random effect for accouning for overdispersion
    sd.spl1 ~ dunif(0,100)
    var.spl1<- sd.spl1^2	
    tau.spl1 <- 1/var.spl1

}"
  
jdat5 <- list(N=nrow(covar.matrix), log_ach_noj=log(data.fit$ach_noj), outcome=data.fit$outcome,t.index=data.fit$t, month.mat=month.mat,spl.mat=spl.mat, pandemic=data.fit$pandemic)
jmod5 <- jags.model(textConnection(jcode_g5), data=jdat5, n.chains=2, n.adapt=1000)
update(jmod5,5000)
jpos5 <- coda.samples(jmod5, c('beta', 'delta', "log.lambda",'log.lambda.novax','log.rr.vax'),	 n.iter=10000, thin=5) 
summary.pred<-summary(jpos5)
#plot(jpos5,ask=TRUE)

pred.ipd.quantiles2 <-as.data.frame(summary.pred[2])
fitted<-exp(pred.ipd.quantiles2[grepl("log.lambda[", rownames(pred.ipd.quantiles2), fixed=TRUE), c(1,3,5)])
counterfact<-exp(pred.ipd.quantiles2[grepl("log.lambda.novax", rownames(pred.ipd.quantiles2), fixed=TRUE), c(1,3,5)])
rr.vax<-exp(pred.ipd.quantiles2[grepl("log.rr.vax", rownames(pred.ipd.quantiles2), fixed=TRUE), c(1,3,5)])

pdf("plots.pdf")
#PLOT PREDICTED VS OBSERVED
matplot(counterfact, type='l',lty=c(2,1,2), col='lightgray', ylim=c(0,max(fitted)*1.2) )
points(fitted[,2], type="l", col='red')
points(data.fit$outcome, type="p")
title("observed vs expected")


#Plot log rate ratio trajectory
matplot(rr.vax, type='l',lty=c(2,1,2), col='lightgray', ylim=c(0.5,1.5))
abline(h=1)
title("Rate ratio trajectory")
dev.off()


