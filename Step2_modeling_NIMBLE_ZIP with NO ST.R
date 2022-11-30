################################################
# Step2: Winbugs modeling
################################################

setwd("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/data/")

memory.limit(size=100000)
library(nimble)
library(coda)
library(dplyr)
library(rgdal)
library(spdep)

#########################
# 1. Load dataset for 2016-2019
#########################
adj<-c(read.csv("Korea_2018_adj_FINAL.csv"))$x
num<-c(read.csv("Korea_2018_num_FINAL.csv"))$x
sumNumNeigh=sum(num)

hepA.dat<-read.csv("hepA_2011_2020.csv") # 250*525 (sido, sigungu, sigungu_cd, X201101, X201102,...)
hepA<-NULL
for(year in 2016:2019){
  temp.dat<-hepA.dat %>% dplyr::select(starts_with(paste("X",year,sep="")))
  if(year==2016) {hepA<-temp.dat}
  if(year!=2016) {hepA<-cbind(hepA,temp.dat)}
}
N.site<-nrow(hepA)
N.time<-ncol(hepA)
years<-c(rep(1,53),rep(2,52),rep(3,52),rep(4,52))
y<-as.matrix(hepA)

##################################
# covariates
##################################
pop.dat<-read.csv("covariates/population_2011_2020.csv") # 250*393
pop.total.id<-paste("X",2016:2019,"년_계_총인구수",sep="")
pop<-as.matrix(pop.dat[,pop.total.id])

popden.dat<-read.csv("covariates/popdensity_2011_2020.csv") 
popden<-round(as.matrix(popden.dat[,paste("X",2016:2019,sep="")]),3)

pop.male.dat<-cbind(pop.dat[,1:3],pop.dat[,paste("X",2016:2019,"년_남_총인구수",sep="")])
male<-round(as.matrix(pop.male.dat[,4:7]/pop*100),3)

#Age0019.rate.dat<-(pop.dat[,paste("X",2016:2019,"년_계_0.9세",sep="")] + pop.dat[,paste("X",2016:2019,"년_계_10.19세",sep="")])/pop*100
Age3049.rate.dat<-(pop.dat[,paste("X",2016:2019,"년_계_30.39세",sep="")] + pop.dat[,paste("X",2016:2019,"년_계_40.49세",sep="")])/pop*100

#Age2039.rate.dat<-(pop.dat[,paste("X",2016:2019,"년_계_20.29세",sep="")] + pop.dat[,paste("X",2016:2019,"년_계_30.39세",sep="")])/pop*100
#Age4059.rate.dat<-(pop.dat[,paste("X",2016:2019,"년_계_40.49세",sep="")] + pop.dat[,paste("X",2016:2019,"년_계_50.59세",sep="")])/pop*100
#Age2059.rate.dat<-Age2039.rate.dat + Age4059.rate.dat

age<-round(as.matrix(Age3049.rate.dat),3)

foreigner.dat<-read.csv("covariates/foreigner_cases_2011_2020.csv")
fore<-as.matrix(foreigner.dat[,paste("reigsterd_foreigner_",2016:2019,sep="")])  # total foreigner 
fore<-round(log(fore),3)

edu<-c(read.csv("covariates/HighEducation_2015.csv")[,"Highedu_rate"])
# income (2016-2019)
income.raw.dat<-read.csv("covariates/earned_income_tax_2016_2019.csv")
income.dat<-income.raw.dat[,paste("X",2016:2019,".급여총계.금액..백만원.",sep="")]/income.raw.dat[,paste("X",2016:2019,".급여총계.인원..명.",sep="")]
income<-round(as.matrix(income.dat),3)

doctor.dat<-read.csv("covariates/doctors_2011_2020.csv")
doctor<-as.matrix(doctor.dat[,paste("doctors_",2016:2019,sep="")])

sewage<-round(as.matrix(read.csv("covariates/public_sewage_2011_2019_v2.csv")[,9:12]),3)
water<-round(as.matrix(read.csv("covariates/water_supply_2011_2019_v2.csv")[,9:12]),3)
birthrate<-round(as.matrix(read.csv("covariates/birthrate_2011_2020_v2.csv")[,9:12]),3)

temp.ave.dat<-read.csv("covariates/temp_ave_2011_2020.csv")
prec.dat<-read.csv("covariates/prec_2011_2020.csv")
humid.dat<-read.csv("covariates/humi_2011_2020.csv")
temp<-NULL; precip<-NULL; humid<-NULL
for(year in 2016:2019){
  temp.ave<-temp.ave.dat %>% dplyr::select(starts_with(paste("X",year,sep="")))
  temp.prec<-prec.dat %>% dplyr::select(starts_with(paste("X",year,sep="")))
  temp.humid<-humid.dat %>% dplyr::select(starts_with(paste("X",year,sep="")))
  if(year==2016) {temp<-temp.ave; precip<-temp.prec; humid<-temp.humid}
  if(year!=2016) {temp<-cbind(temp,temp.ave); precip<-cbind(precip,temp.prec); humid<-cbind(humid,temp.humid)}
}

temp<-round(as.matrix(temp),3)
precip<-round(as.matrix(precip),3)
humid<-round(as.matrix(humid),3)


#######################################
# 2. NIMBLE - working on!! 
# 2-1. No ST ZIP model
#######################################
model.ID<-"ZIP with NO ST_v2"


# the model
Model <- nimbleCode({ 
  for(i in 1:N.site) {
    for(j in 1:N.time) {
       y[i,j] ~ dpois(mu[i,j])
       mu[i,j] <- (1-zr[i,j])*lambda[i,j]
       zr[i,j] ~ dbern(p0[i,j])
       log(lambda[i,j]) <- beta0 + beta[1] * age[i,years[j]] + beta[2]*fore[i,years[j]] + beta[3]*edu[i]+beta[4]*income[i,years[j]] + beta[5]*sewage[i,years[j]] + beta[6]*water[i,years[j]]+beta[7]*birthrate[i,years[j]] + beta[8]*male[i,years[j]] + beta[9]*doctor[i,years[j]] + alpha[1] * temp[i,j]+ alpha[2]* precip[i,j]+alpha[3]*humid[i,j]+log(popden[i,years[j]])+ ST[i,j] 
       logit_p0[i,j] <- gamma0 + gamma[1] * age[i,years[j]] + gamma[2]*fore[i,years[j]] + gamma[3]*edu[i]+gamma[4]*income[i,years[j]] + gamma[5]*sewage[i,years[j]] + gamma[6]*water[i,years[j]]+gamma[7]*birthrate[i,years[j]]  + gamma[8]*male[i,years[j]] + gamma[9]*doctor[i,years[j]] + delta[1] * temp[i,j]+ delta[2]* precip[i,j]+delta[3]*humid[i,j]
       p0[i,j] <- expit(logit_p0[i,j])
    
     ST[i,j]~dnorm(0,sd=10)
  }}


beta0~dnorm(0, sd=10)
gamma0~dnorm(0, sd=10)

  for(k in 1:9) {
    beta[k]~dnorm(0, sd=10)
    gamma[k]~dnorm(0, sd=10)
  }
  for(k in 1:3) {
    alpha[k]~dnorm(0, sd=10)
    delta[k]~dnorm(0, sd=10)
  }
})

nimbleData <- list(y=y,age=age,fore=fore,edu=edu,income=income,sewage=sewage,water=water,birthrate=birthrate,male=male,doctor=doctor,temp=temp,precip=precip,humid=humid,popden=popden)
Consts <- list(N.site=250,N.time=209,years=years)
inits <- list(beta0=-0.1,gamma0=-0.1,beta=rep(0,9),alpha=rep(0,3),gamma=rep(0,9),delta=rep(0,3),ST=matrix(0,nrow=N.site,ncol=N.time),zr=structure(.Data=rep(0, N.site*N.time), .Dim=c(N.site,N.time)))

nimbleModel <- nimbleModel(code = Model, name = "nimbleModel",constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("beta0","gamma0","beta","alpha","gamma","delta","mu")) 

nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(nimbleModel)
CnbMCMC <- compileNimble(nbMCMC, project = nimbleModel)

# Sampling
T <- Sys.time()
results <- runMCMC(CnbMCMC,niter=50000,nburnin=0,inits=inits,nchains=1,thin=50,samplesAsCodaMCMC = T) # check for convergence (2/3/11:30 AM - start) : 1.2 days
Sys.time() - T

# file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/%s_summary.csv",model.ID)

T <- Sys.time()
for (i in 1:100){ 
  results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
  write.csv(as.matrix(results), file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
  rm(results)
}
Sys.time() - T
# for 23 hours

## results summary
est1 <- read.csv(file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/Output_%s_%d.csv",model.ID,1))

est <- data.frame(matrix(NA, 5000, ncol(est1)))
for (i in 1:100) {
  est[(1+50*(i-1)):(50*i),] <- read.csv(file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/Output_%s_%d.csv",model.ID,i))
}
colnames(est) <- colnames(read.csv(file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/Output_%s_%d.csv",model.ID,i)))

summary <- data.frame("mean"=rep(NA, ncol(est1)), "med"=rep(NA, ncol(est1)), "95%CI_L"=rep(NA, ncol(est1)), "95%CI_U"=rep(NA, ncol(est1)), "sd"=rep(NA, ncol(est1)))
summary$mean <- round(apply(est, 2, mean),4)
for (i in 1:ncol(est1)) {
  summary$med[i] <- round(quantile(est[,i], probs=0.5),4)
  summary$X95.CI_L[i] <- round(quantile(est[,i], probs=0.025),4)
  summary$X95.CI_U[i] <- round(quantile(est[,i], probs=0.975),4)
  summary$sd[i] <- round(sd(est[,i]),3)
}
rownames(summary) <- colnames(est)
View(summary)
write.csv(summary, file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/%s_summary.csv",model.ID))

mu.ind <- which(colnames(est)=="mu.1..1.")
mu <- est[,mu.ind:(mu.ind+(N.site*N.time)-1)]

y.pred <- matrix(NA, 5000, (N.time*N.site))
for (i in 1:(N.time*N.site)) {
  for (j in 1:5000) {
    y.pred[j,i] <- rpois(1, mu[j,i])
  }
}

y.pred.mean <- apply(y.pred, 2, mean)
y.pred.mean.m<-matrix(y.pred.mean,nrow=N.site,ncol=N.time)
write.csv(y.pred.mean.m, file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/%s_ypred.csv",model.ID),row.names=FALSE)

y <- unlist(y)
(MAE <- mean(abs(y-y.pred.mean))) # MAE

(MSPE <- mean((y-y.pred.mean)^2)) # MSPE

# DIC (2Dbar-D(est))
Dbar.samp <- rep(NA, 5000)
for (i in 1:5000) {Dbar.samp[i] <- -2*sum(-mu[i,]+y*log(mu[i,]+0.01)-lfactorial(y))} # 666022: clustering w/o spatial RE
Dbar <- mean(Dbar.samp)

mu.est <- apply(mu, 2, mean)
D.est <- -2*sum(-mu.est+y*log(mu.est+0.01)-lfactorial(y))
(DIC <- 2*Dbar - D.est)
(pD <- Dbar-D.est)

measure.values<-c(MAE,MSPE,Dbar,DIC,pD)
names(measure.values)<-c("MAE","MSPE","Dbar","DIC","pD")
write.csv(measure.values, file=sprintf("D:/SynologyDrive/Research/Collaboration/Prof_Kim_Jeong/HepA modeling/output/%s_measures.csv",model.ID))





