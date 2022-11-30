################################################
# Modeling code for ZIP with Type 1 & Type 3
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
# 2. ZIP with Type1 : NIMBLE code
#######################################
model.ID<-"ZIP with Type1"

# the model
Model <- nimbleCode({ 
  for(i in 1:N.site) {
    for(j in 1:N.time) {
      y[i,j] ~ dpois(mu[i,j])
      mu[i,j] <- (1-zr[i,j])*lambda[i,j]
      zr[i,j] ~ dbern(p0[i,j])
      log(lambda[i,j]) <- beta0 + beta[1] * age[i,years[j]] + beta[2]*fore[i,years[j]] + beta[3]*edu[i]+beta[4]*income[i,years[j]] + beta[5]*sewage[i,years[j]] + beta[6]*water[i,years[j]]+beta[7]*birthrate[i,years[j]] + beta[8]*male[i,years[j]] + beta[9]*doctor[i,years[j]] + alpha[1] * temp[i,j]+ alpha[2]* precip[i,j]+alpha[3]*humid[i,j]+log(popden[i,years[j]]) 
      logit_p0[i,j] <- gamma0 + gamma[1] * age[i,years[j]] + gamma[2]*fore[i,years[j]] + gamma[3]*edu[i]+gamma[4]*income[i,years[j]] + gamma[5]*sewage[i,years[j]] + gamma[6]*water[i,years[j]]+gamma[7]*birthrate[i,years[j]]  + gamma[8]*male[i,years[j]] + gamma[9]*doctor[i,years[j]] + delta[1] * temp[i,j]+ delta[2]* precip[i,j]+delta[3]*humid[i,j]
      p0[i,j] <- expit(logit_p0[i,j])
      
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


nimbleData <- list(y=y,age=age,popden=popden, fore=fore,edu=edu,income=income,sewage=sewage,water=water,birthrate=birthrate,male=male,doctor=doctor,temp=temp,precip=precip,humid=humid)
Consts <- list(N.site=250,N.time=209,years=years)
inits <- list(beta0=-0.1,gamma0=-0.1,beta=rep(0,9),alpha=rep(0,3),gamma=rep(0,9),delta=rep(0,3),zr=structure(.Data=rep(0, N.site*N.time), .Dim=c(N.site,N.time)))

nimbleModel <- nimbleModel(code = Model, name = "nimbleModel",constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("beta0","gamma0","beta","alpha","gamma","delta","mu")) 

nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(nimbleModel)
CnbMCMC <- compileNimble(nbMCMC, project = nimbleModel)

# Sampling for burn-in
results <- runMCMC(CnbMCMC,niter=50000,nburnin=0,inits=inits,nchains=1,thin=50,samplesAsCodaMCMC = T) # check for convergence (4/15 2PM - start) : 1.2 days

# sampling
i=1
results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
est <- data.frame(matrix(NA, 5000, ncol(as.matrix(results))))
est[(1+50*(i-1)):(50*i),] <-as.matrix(results)

for (i in 2:100){ 
  results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
  write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
  est[(1+50*(i-1)):(50*i),] <-as.matrix(results)
  rm(results)
}
colnames(est) <- colnames(read.csv(file=sprintf("D:/Output_%s_%d.csv",model.ID,1)))

# Save summary 
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

#####################################################
# 3. Using residual, fit ST components
#####################################################
summary.dat<-read.csv(file=sprintf("D:/%s_summary.csv",model.ID))
mu.ind <- which(summary.dat$X=="mu.1..1.")
mu.est<- matrix(summary.dat$mean[mu.ind:(mu.ind+(N.site*N.time)-1)],nrow=N.site,ncol=N.time)   # N.site * N.time
res<-log(y+0.1)-log(mu.est)

model.ID<-"ZIP with Type1_residual"

# the model
Model <- nimbleCode({ 
  for(i in 1:N.site) {
    for(j in 1:N.time) {
      res[i,j] ~ dnorm(ST[i,j],sd=sigma.res)
      ST[i,j]<-alpha+u[i]+v[i]+t[j]+eta[j]
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
  sigma.res~dunif(0,10)
  alpha~dnorm(0, sd=10)
  
})

nimbleData <- list(res=res)
Consts <- list(N.site=250,N.time=209,sumNumNeigh=sumNumNeigh,adj=adj,num=num)
inits <- list(alpha=0.1,u=rep(0,N.site),v=rep(0,N.site),t=rep(0,N.time),eta=rep(0,N.time),sigma.v=0.1,sigma.u=0.1,sigma.t=0.1,sigma.eta=0.1,sigma.res=0.1)

nimbleModel <- nimbleModel(code = Model, name = "nimbleModel",constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("alpha","sigma.u","sigma.v","sigma.t","sigma.eta","u","v","t","eta","ST")) 

nbMCMC <- buildMCMC(MCMCconfig)
Cnb <- compileNimble(nimbleModel)
CnbMCMC <- compileNimble(nbMCMC, project = nimbleModel)

# burn-in
results <- runMCMC(CnbMCMC,niter=50000,nburnin=0,inits=inits,nchains=1,thin=50,samplesAsCodaMCMC = T) 

# sampling
i=1
results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
est <- data.frame(matrix(NA, 5000, ncol(as.matrix(results))))
est[(1+50*(i-1)):(50*i),] <-as.matrix(results)

for (i in 2:100){ 
  results <- runMCMC(CnbMCMC,niter=2500,nchains=1,thin=50,samplesAsCodaMCMC = T)
  write.csv(as.matrix(results), file=sprintf("D:/Output_%s_%d.csv",model.ID,i),row.names=FALSE)
  est[(1+50*(i-1)):(50*i),] <-as.matrix(results)
  rm(results)
}
colnames(est) <- colnames(read.csv(file=sprintf("D:/Output_%s_%d.csv",model.ID,1)))

# Save summary 
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




#####################################################
# 4. ZIP_Type3
#####################################################
model.ID<-"ZIP with Type1_residual"
summary.dat<-read.csv(file=sprintf("D:/%s_summary.csv",model.ID))

ST.ind <- which(summary.dat$X=="ST.1..1.")
ST.est<- matrix(summary.dat$mean[ST.ind:(ST.ind+(N.site*N.time)-1)],nrow=N.site,ncol=N.time)   # N.site * N.time
alpha.ind <- which(summary.dat$X=="alpha")
alpha.est<-summary.dat$mean[alpha.ind]
ST.new<-ST.est-alpha.est

model.ID<-"ZIP with Type3"

# the model
Model <- nimbleCode({ 
  for(i in 1:N.site) {
    for(j in 1:N.time) {
      y[i,j] ~ dpois(mu[i,j])
      mu[i,j] <- (1-zr[i,j])*lambda[i,j]
      zr[i,j] ~ dbern(p0[i,j])
      log(lambda[i,j]) <- beta0 + beta[1] * age[i,years[j]] + beta[2]*fore[i,years[j]] + beta[3]*edu[i]+beta[4]*income[i,years[j]] + beta[5]*sewage[i,years[j]] + beta[6]*water[i,years[j]]+beta[7]*birthrate[i,years[j]] + beta[8]*male[i,years[j]] + beta[9]*doctor[i,years[j]] + alpha[1] * temp[i,j]+ alpha[2]* precip[i,j]+alpha[3]*humid[i,j]+log(popden[i,years[j]])+ST[i,j]+err[i,j]   
      logit_p0[i,j] <- gamma0 + gamma[1] * age[i,years[j]] + gamma[2]*fore[i,years[j]] + gamma[3]*edu[i]+gamma[4]*income[i,years[j]] + gamma[5]*sewage[i,years[j]] + gamma[6]*water[i,years[j]]+gamma[7]*birthrate[i,years[j]]  + gamma[8]*male[i,years[j]] + gamma[9]*doctor[i,years[j]] + delta[1] * temp[i,j]+ delta[2]* precip[i,j]+delta[3]*humid[i,j]
      p0[i,j] <- expit(logit_p0[i,j])
      err[i,j]~dnorm(0, sd=sigma.err)
      
    }}
  
  sigma.err~dunif(0,10)
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

nimbleData <- list(y=y,age=age,fore=fore,edu=edu,income=income,sewage=sewage,water=water,birthrate=birthrate,male=male,doctor=doctor,temp=temp,precip=precip,humid=humid,popden=popden,ST=ST.new)
Consts <- list(N.site=250,N.time=209,years=years)
inits <- list(beta0=-0.1,gamma0=-0.1,beta=rep(0,9),alpha=rep(0,3),gamma=rep(0,9),delta=rep(0,3),sigma.err=0.1, err=matrix(0,nrow=N.site,ncol=N.time), zr=structure(.Data=rep(0, N.site*N.time), .Dim=c(N.site,N.time)))

nimbleModel <- nimbleModel(code = Model, name = "nimbleModel",constants = Consts, data = nimbleData, inits = inits)
MCMCconfig <- configureMCMC(nimbleModel,monitors=c("beta0","gamma0","beta","alpha","gamma","delta","sigma.err","err","mu")) 

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
