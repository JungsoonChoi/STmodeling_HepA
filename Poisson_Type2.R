################################################
# Modeling code for Poisson with Type 2
# in the paper "Investigating the spatio-temporal variation of Hepatitis A in Korea using a Bayesian model"
#
# Jungsoon Choi
################################################

library(nimble)
library(coda)
library(dplyr)
library(rgdal)
library(spdep)

#########################
# 1. INPUT: 
#
# N.site: number of areas
# N.time: number of time points
# years: yearly index corresponding of time point
#
# y (N.site*N.time): cases
# age, fore, income,sewage,water,birthrate,male,doctor (N.site * number of years): covariates 
# edu: N.site
# temp, precip, humid (N.site * N.time): water covariate  
# adj, num, simNumNeigh: spatial information
#########################
# Data should be given


#######################################
# 2. Poisson with Type2 : NIMBLE code
#######################################
model.ID<-"Poi with Type2"

# the model
Model <- nimbleCode({ 
  for(i in 1:N.site) {
    for(j in 1:N.time) {
      y[i,j] ~ dpois(mu[i,j])
      log(mu[i,j]) <- beta0 + beta[1] * age[i,years[j]] + beta[2]*fore[i,years[j]] + beta[3]*edu[i]+beta[4]*income[i,years[j]] + beta[5]*sewage[i,years[j]] + beta[6]*water[i,years[j]] + beta[7]*birthrate[i,years[j]] + beta[8]*male[i,years[j]] + beta[9]*doctor[i,years[j]]+ alpha[1] * temp[i,j]+ alpha[2]* precip[i,j]+alpha[3]*humid[i,j]+log(popden[i,years[j]])+ ST[i,j] 

      ST[i,j]<-u[i]+v[i]+t[j]+eta[j]
    }}
  
  for (k in 1:sumNumNeigh) {weights[k] <- 1}
  v[1:N.site] ~ dcar_normal(adj[1:sumNumNeigh], weights[1:sumNumNeigh], num[1:N.site], tau.v, zero_mean = 1)
  tau.v <- 1/pow(sigma.v, 2)
  sigma.v ~ dunif(0,10)
  
  t[1]~dnorm(0,sd=sigma.t)
  for(j in 2:N.time){
    t[j]~dnorm(t[j-1],sd=sigma.t)
  }
  sigma.t~ dunif(0,10)
  
  for(i in 1:N.site) {
    u[i]~dnorm(0, sd=sigma.u)
  }
  sigma.u ~ dunif(0,10)
  
  for(j in 1:N.time) {
    eta[j]~dnorm(0, sd=sigma.eta)
  }
  sigma.eta ~ dunif(0,10)
  beta0~dnorm(0, sd=10)
  
  for(k in 1:9) {
    beta[k]~dnorm(0, sd=10)
    }
  for(k in 1:3) {
    alpha[k]~dnorm(0, sd=10)
    }
})

nimbleData <- list(y=y,age=age,popden=popden, fore=fore,edu=edu,income=income,sewage=sewage,water=water,birthrate=birthrate,male=male,doctor=doctor,temp=temp,precip=precip,humid=humid)
Consts <- list(N.site=250,N.time=209,years=years,sumNumNeigh=sumNumNeigh,adj=adj,num=num)
inits <- list(beta0=-0.1,beta=rep(0,9),alpha=rep(0,3),
              u=rep(0,N.site),v=rep(0,N.site),t=rep(0,N.time),eta=rep(0,N.time),sigma.v=0.1,sigma.u=0.1,sigma.t=0.1,sigma.eta=0.1)
nimbleModel <- nimbleModel(code = Model, name = "nimbleModel",constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("beta0","beta","alpha","sigma.u","sigma.v","sigma.t","sigma.eta","u","v","t","eta","mu")) 
nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(nimbleModel)
CnbMCMC <- compileNimble(nbMCMC, project = nimbleModel)

# Sampling for burn-in
results <- runMCMC(CnbMCMC,niter=50000,nburnin=0,inits=inits,nchains=1,thin=50,samplesAsCodaMCMC = T) 

# sampling
i=1
results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
est <- data.frame(matrix(NA, 5000, ncol(as.matrix(results))))
est[(1+50*(i-1)):(50*i),] <-as.matrix(results)
rm(results)

for (i in 2:100){ 
  results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
  write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
  est[(1+50*(i-1)):(50*i),] <-as.matrix(results)
  rm(results)
}

colnames(est) <- colnames(read.csv(file=sprintf("D:/Output_%s_%d.csv",model.ID,1)))

summary <- data.frame("mean"=rep(NA, ncol(est)), "med"=rep(NA, ncol(est)), "95%CI_L"=rep(NA, ncol(est)), "95%CI_U"=rep(NA, ncol(est)), "sd"=rep(NA, ncol(est)))
summary$mean <- round(apply(est, 2, mean),4)
for (i in 1:ncol(est)) {
  summary$med[i] <- round(quantile(est[,i], probs=0.5),4)
  summary$X95.CI_L[i] <- round(quantile(est[,i], probs=0.025),4)
  summary$X95.CI_U[i] <- round(quantile(est[,i], probs=0.975),4)
  summary$sd[i] <- round(sd(est[,i]),3)
}
rownames(summary) <- colnames(est)
write.csv(summary, file=sprintf("D:/%s_summary.csv",model.ID))


