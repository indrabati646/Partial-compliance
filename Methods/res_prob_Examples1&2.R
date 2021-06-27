#Gives the Stage-1 response probability for Examples 1 and 2
# We use the 'nimble' package in R which lets run BUGS model from R

source("bnp_model_Examples1&2.R")

library(nimble)


response_prob <- function(data,result){
  data <- result$dt
  
  indy1 <- which(data[,"a1"]==1 & data[,"s"]==1)
  indy2 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1)
  indy3 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1)
  indy4 <- which(data[,"a1"]==-1 & data[,"s"]==1)
  indy5 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1)
  indy6 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1)
  
  data[indy1,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD1
  data[indy2,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD2
  data[indy3,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD3
  data[indy4,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD4
  data[indy5,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD5
  data[indy6,c("D1","D2","D3","D4")] <- result$ME[[length(result$ME)]]$DD6

data0 <- data[data[,"a1"]==1,]
n0 <- nrow(data0)
a1 <- data0[,"a1"]; frac <- data0[,"s"]
frac[frac==-1] <- 0; D1 <- data0[,"D1"]; D2 <- data0[,"D2"]
D3 <- data0[,"D3"]
X1 <- data0[,"X1"]; X2<- data0[,"X2"]; X3 <- data0[,"X3"]

code<- nimbleCode({
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  alpha~ dgamma(1,1)
  w[1:L] <- stick_breaking(v[1:(L-1)])
  for(i in 1:n) {
    z[i] ~ dcat(w[1:L])
    logit(p[i])<- alpha1[z[i]]+beta11[z[i]]*D1[i]+beta12[z[i]]*D2[i]+
      beta3[z[i]]*X1[i]+beta4[z[i]]*X2[i]+beta5[z[i]]*X3[i]
    frac[i]~ dbern(p[i])
  }
  for(i in 1:L){
    alpha1[i]~dnorm(1,.5)
    beta11[i]~dnorm(1,.5)
    beta12[i]~dnorm(1,.5)
    beta3[i]~dnorm(1,.5)
    beta4[i]~dnorm(1,.5)
    beta5[i] ~dnorm(1,.5)
  }
  
})
#D11<- D1[frac==1]; D12 <- D1[frac==-1]; n1 <- length(D11)
#n2<- length(D12)
dat1 <- list(frac=frac,D1=D1,D2=D2,X1=X1,X2=X2,X3=X3)
constants <- list(n=n0, L=5)
inits <- list(alpha1=rnorm(5,1,1),beta11=rnorm(5,1,1),beta12=rnorm(5,1,0.5),
              z = sample(1:4, size = constants$n, replace = TRUE),
              beta3=rnorm(5,1,1),beta4=rnorm(5,1,1),beta5=rnorm(5,1,1),
              alpha=1,v  = rbeta(constants$L-1, 1, 1))
pump <- nimbleModel(code = code, name = "pump", 
                    constants = constants,
                    data = dat1, inits = inits)
cpump <- compileNimble(pump, showCompilerOutput = T)
pumpConf <- configureMCMC(pump, print = TRUE)
pumpConf$addMonitors(c("alpha1", "beta11","beta12","beta3","beta4","beta5",
                       "alpha","z","v","w"))
pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump)
niter <- 100+length(result$ME)
set.seed(1)
samples <- runMCMC(CpumpMCMC, niter = niter,nburnin=100)


data1 <- data[data[,"a1"]==-1,]
n1 <- nrow(data1)
a1 <- data1[,"a1"]; frac <- data1[,"s"]
frac[frac==-1] <- 0; D1 <- data1[,"D1"]; D2 <- data1[,"D2"]
X1 <- data1[,"X1"]; X2<- data1[,"X2"]; X3 <- data1[,"X3"]

code<- nimbleCode({
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  for(i in 1:L){
    alpha1[i]~dnorm(1,0.5)
    beta21[i]~dnorm(1,0.5)
    beta22[i]~dnorm(1,0.5)
    beta3[i]~dnorm(1,0.5)
    beta4[i]~dnorm(1,0.5)
    beta5[i] ~dnorm(1,0.5)
  }
  alpha~ dgamma(1,1)
  w[1:L] <- stick_breaking(v[1:(L-1)])
  for(i in 1:n) {
    z[i] ~ dcat(w[1:L])
    logit(p[i])<- alpha1[z[i]]+beta21[z[i]]*D1[i]+beta22[z[i]]*D2[i]+
      beta3[z[i]]*X1[i]+beta4[z[i]]*X2[i]+beta5[z[i]]*X3[i]
    frac[i]~ dbern(p[i])
  }
})
#D11<- D1[frac==1]; D12 <- D1[frac==-1]; n1 <- length(D11)
#n2<- length(D12)
dat1 <- list(frac=frac,D1=D1,D2=D2,X1=X1,X2=X2,X3=X3)
constants <- list(n=n1, L=5)
inits <- list(alpha1=rnorm(5,1,1),beta21=rnorm(5,1,1),beta22=rnorm(5,1,1),
              z = sample(1:4, size = constants$n, replace = TRUE),
              beta3=rnorm(5,1,1),beta4=rnorm(5,1,1),beta5=rnorm(5,1,1),
              alpha=1,v  = rbeta(constants$L, 1, 1))
pump <- nimbleModel(code = code, name = "pump", 
                    constants = constants,
                    data = dat1, inits = inits)
cpump <- compileNimble(pump, showCompilerOutput = T)
pumpConf <- configureMCMC(pump, print = TRUE)
pumpConf$addMonitors(c("alpha1", "beta21","beta22","beta3","beta4","beta5",
                       "alpha","z","v","w"))
pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump)
niter <- 100+length(result$ME)
set.seed(1)
samples1 <- runMCMC(CpumpMCMC, niter = niter,nburnin=100)
return(list(samples=samples,samples1=samples1))
}

