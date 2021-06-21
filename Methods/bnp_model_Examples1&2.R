#Function for fitting the nonparametric Bayesian model for a single data of the type in Examples 1 and 2. 
# Examples 1 and 2 satisfy the assumption D3=D4.
# Some ideas have been borrowed from Kim et. al (2019).

library(dirichletprocess)
library(mnormt)
library(rootSolve)
library(matrixcalc)
library(condMVNorm)
library(KScorrect)
library(MASS)
library(mvtnorm)
library(truncnorm)

source("metropolis_For_D.R")
source("metropolis_For_R.R")

bnp_mod <- function(data){
  
  #Pre-processing the data-frame for analysis
  data <- matrix(unlist(data),nrow=nrow(data),ncol=ncol(data))
  colnames(data) <- c("a1","a2","s","D1","D2","D3","y","y1","X1","X2","X3")
  D4 <- 0
  data <- cbind(data,D4)
  data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1),
       "D4"] <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1),
                     "D3"]
  
  
  indy1 <- which(data[,"a1"]==1 & data[,"s"]==1)
  indy2 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1)
  indy3 <- which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1)
  indy4 <- which(data[,"a1"]==-1 & data[,"s"]==1)
  indy5 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1)
  indy6 <- which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1)
  
  data[c(indy4,indy5,indy6),"D1"] <- 0
  data[c(indy1,indy2,indy3), "D2"] <- 0
  data[c(indy1,indy3,indy4,indy5,indy6),"D3"] <- 0
  data[c(indy1,indy2,indy3,indy4,indy6),"D4"] <- 0
  
  
  x1 <- data[,"X1"]; x2 <- data[,"X2"]; x3 <- data[,"X3"]
  x <- cbind(x1,x2,x3)
  a1<-data[,"a1"]
  
  y1 <- data[which(data[,"a1"]==1 & data[,"s"]==1),]
  y2 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1),]
  y3 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1),]
  y4 <- data[which(data[,"a1"]==-1 & data[,"s"]==1),]
  y5 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1),]
  y6 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1),]
  
  xx1 <- y1[,c("X1","X2","X3")]
  xx2 <- y2[,c("X1","X2","X3")]
  xx3 <- y3[,c("X1","X2","X3")]
  xx4 <- y4[,c("X1","X2","X3")]
  xx5 <- y5[,c("X1","X2","X3")]
  xx6 <- y6[,c("X1","X2","X3")]
  
  aa1<-y1[,"a1"]; aa2<-y2[,"a1"]; aa3<-y3[,"a1"]
  aa4<-y4[,"a1"]; aa5<-y5[,"a1"]; aa6<-y6[,"a1"]
  
  n1 <- nrow(y1); n2 <- nrow(y2); n3 <- nrow(y3); n4 <- nrow(y4) 
  n5 <- nrow(y5); n6 <- nrow(y6)
  n <- n1+n2+n3+n4+n5+n6 
  y <- data[,"y"]
  yy1 <- y1[,"y"]; yy2 <- y2[,"y"]; yy3 <- y3[,"y"]; yy4 <- y4[,"y"]
  yy5 <- y5[,"y"]; yy6 <- y6[,"y"]
  
  DD1 <- y1[,c("D1","D2","D3","D4")]
  DD2 <- y2[,c("D1","D2","D3","D4")]
  DD3 <- y3[,c("D1","D2","D3","D4")]
  DD4 <- y4[,c("D1","D2","D3","D4")]
  DD5 <- y5[,c("D1","D2","D3","D4")]
  DD6 <- y6[,c("D1","D2","D3","D4")]
  DD1[,2:4] <- cbind(rtruncnorm(n1,mean(DD1[,1]),0.1,a=0,b=1),0,0)
  DD2[,c(2,4)] <- cbind(rtruncnorm(n2,mean(DD2[,1]),0.1,a=0,b=1),0)
  DD3[,2:4] <- cbind(rtruncnorm(n3,mean(DD3[,1]),0.1,a=0,b=1),
                     rtruncnorm(n3,mean(DD3[,1]),0.1,a=0,b=1),0)
  DD4[,c(1,3,4)] <-cbind(rtruncnorm(n4,mean(DD4[,2]),0.1,a=0,b=1),0,
                         rtruncnorm(n4,mean(DD4[,2]),0.1,a=0,b=1))
  DD5[,c(1,3)] <- cbind(rtruncnorm(n5,mean(DD5[,2]),0.1,a=0,b=1),0)
  DD6[,c(1,3,4)] <- cbind(rtruncnorm(n6,mean(DD6[,2]),0.1,a=0,b=1),0,
                          rtruncnorm(n6,mean(DD6[,2]),0.1,a=0,b=1))
  D <- rbind(DD1,DD2,DD3,DD4,DD5,DD6)
  
  K <-3; mcmc <- 10000 ; P<-2
  dim.cov <- 2;dimx <- dim.cov+1
  y_1<-y_2<-y_3<-y_4<-y_5<-y_6<- list()
  R<- diag(4)
  
  #Initial values for the parameters of the marginals of observed compliances
  para.dd1 <- matrix(nrow=mcmc,ncol=K+K+3+K+5+1+n)
  para.dd1[1,] <-para.dd1[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                                                                    rbind(xx1,xx2,xx3)[,2]+rbind(xx1,xx2,xx3)[,3]))[1],K),
                                   coefficients(lm(D[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                                                     rbind(xx1,xx2,xx3)[,2]+rbind(xx1,xx2,xx3)[,3]))[-1],
                                   rep(0.1,K),1,1,1,mean(D[1:(n1+n2+n3),1]),1,1, 
                                   sample(1:K,n, replace=TRUE))
  
  
  para.dd2 <- matrix(nrow=mcmc,ncol=K+K+3+K+5+1+n)
  para.dd2[1,] <-para.dd2[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                                                                    rbind(xx4,xx5,xx6)[,2]+rbind(xx4,xx5,xx6)[,3]))[1],K),
                                   coefficients(lm(D[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                                                     rbind(xx4,xx5,xx6)[,2]+rbind(xx4,xx5,xx6)[,3]))[-1],
                                   rep(0.1,K),1,1,1,mean(D[(n1+n2+n3+1):n,2]),1,1, 
                                   sample(1:K,n, replace=TRUE))
  
  
  para.dd3 <- matrix(nrow=mcmc,ncol=K+K+4+K+5+1+n)
  para.dd3[1,] <-para.dd3[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[(n1+1):(n1+n2),3]~xx2[,1]+xx2[,2]+xx2[,3]
                                                                  +D[(n1+1):(n1+n2),1]))[1],K),
                                   coefficients(lm(D[(n1+1):(n1+n2),3]~xx2[,1]+xx2[,2]+xx2[,3]
                                                   +D[(n1+1):(n1+n2),1]))[-1],
                                   rep(0.1,K),1,1,1,mean(DD2[,3]),1,1, 
                                   sample(1:K, n, replace=TRUE))
  
  para.dd4 <- matrix(nrow=mcmc,ncol=K+K+4+K+5+1+n)
  para.dd4[1,] <-para.dd4[2,] <- c(rep(1/K,K),rep(coefficients(lm(D[c(n1+n2+n3+n4+1):(n-n6),4]~xx5[,1]+xx5[,2]+xx5[,3]
                                                                  +D[(n1+n2+n3+n4+1):(n-n6),2]))[1],K),
                                   coefficients(lm(D[(n1+n2+n3+n4+1):(n-n6),4]~xx5[,1]+xx5[,2]+xx5[,3]
                                                   +D[(n1+n2+n3+n4+1):(n-n6),2]))[-1],
                                   rep(0.1,K),1,1,1,mean(DD5[,4]),1,1, 
                                   sample(1:K, n, replace=TRUE))
  
  # Copula correlation matrix
  para.C<-matrix(nrow = mcmc, ncol = 12)
  para.C[1,] <- para.C[2,] <- c(rep(NA,6),rep(0,6))
  
  ind1 <- (K+1):(2*K) # intercept
  ind3<- (2*K+1):(2*K+3)
  ind2 <- (2*K+4):(3*K+3) # variance
  ind4 <- 3*K+4 # alpha
  ind5 <- ind4+1 # alpha_beta0
  ind6 <- ind5+1 # alpha_sigma
  ind7 <- ind6+1 # mu_beta0
  ind8 <- ind7+1 # sigma_beta0
  
  ind9 <- length(para.dd1[1,])
  
  ind13<- (2*K+1):(2*K+4)
  ind12 <- (2*K+5):(3*K+4) # variance
  ind14 <- 3*K+5 # alpha
  ind15 <- ind14+1 # alpha_beta0
  ind16 <- ind15+1 # alpha_sigma
  ind17 <- ind16+1 # mu_beta0
  ind18 <- ind17+1 # sigma_beta0
  
  ind19 <- length(para.dd3[1,])
  
  
  cov.d11 <- diag(vcov(lm(D[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                            rbind(xx1,xx2,xx3)[,2]+rbind(xx1,xx2,xx3)[,3]))[1,1],K)
  cov.d12 <- diag(vcov(lm(D[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                            rbind(xx4,xx5,xx6)[,2]+rbind(xx4,xx5,xx6)[,3]))[1,1],K)
  cov.d13 <- diag(vcov(lm(D[(n1+1):(n1+n2),3]~xx2[,1]+xx2[,2]+xx2[,3]+D[(n1+1):(n1+n2),1]))[1,1],K)
  cov.d14 <- diag(vcov(lm(D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),4]~xx5[,1]+
                            xx5[,2]+xx5[,3]+D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),2]))[1,1],K)
  
  cov2.d11 <- diag(diag(vcov(lm(D[1:(n1+n2+n3),1]~rbind(xx1,xx2,xx3)[,1]+
                                  rbind(xx1,xx2,xx3)[,2]+rbind(xx1,xx2,xx3)[,3]))[-1,-1]),3)
  cov2.d12 <- diag(diag(vcov(lm(D[(n1+n2+n3+1):n,2]~rbind(xx4,xx5,xx6)[,1]+
                                  rbind(xx4,xx5,xx6)[,2]+rbind(xx4,xx5,xx6)[,3]))[-1,-1]),3)
  cov2.d13 <- diag(diag(vcov(lm(D[(n1+1):(n1+n2),3]~xx2[,1]+xx2[,2]+xx2[,3]+D[(n1+1):(n1+n2),1]))[-1,-1]),4)
  cov2.d14 <- diag(diag(vcov(lm(D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),4]~xx5[,1]+
                                  xx5[,2]+xx5[,3]+D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),2]))[-1,-1]),4)
  
  
  ME <- list()
  h <- D; h1 <- DD1; h2<- DD2
  h3<-DD3; h4<-DD4; h5<- DD5; h6<- DD6
  
  w1 <- cbind(yy1,h1[,1:2],xx1); wbar1 <- apply(w1,2,mean)
  w2 <- cbind(yy2,h2[,1:3],xx2); wbar2 <- apply(w2,2,mean)
  w3 <- cbind(yy3,h3[,1:3],xx3); wbar3 <- apply(w3,2,mean)
  w4 <- cbind(yy4,h4[,1:2],xx4); wbar4 <- apply(w4,2,mean)
  w5 <- cbind(yy5,h5[,c(1,2,4)],xx5); wbar5 <- apply(w5,2,mean)
  w6 <- cbind(yy6,h6[,c(1,2,4)],xx6); wbar6 <- apply(w6,2,mean)
  
  wcov1 <- cov(w1); wcov2 <- cov(w2); wcov3 <- cov(w3)
  wcov4 <- cov(w4); wcov5 <- cov(w5); wcov6 <- cov(w6)
  
  #Fitting DPM for the joint distribution of (Y,H,X0)
  g0priors1 <- list(mu0=wbar1,Lambda=10*diag(6),kappa0=6,nu=6)
  dp1 <- Fit(DirichletProcessMvnormal(y=w1,g0Priors = g0priors1),100)
  
  g0priors2 <- list(mu0=wbar2,Lambda=10*diag(7),kappa0=7,nu=7)
  dp2 <- Fit(DirichletProcessMvnormal(y=w2,g0Priors = g0priors2),100)
  
  g0priors3 <- list(mu0=wbar3,Lambda=10*diag(7),kappa0=7,nu=7)
  dp3 <- Fit(DirichletProcessMvnormal(y=w3,g0Priors = g0priors3),100)
  
  g0priors4 <- list(mu0=wbar4,Lambda=10*diag(6),kappa0=6,nu=6)
  dp4 <- Fit(DirichletProcessMvnormal(y=w4,g0Priors = g0priors4),100)
  
  g0priors5 <- list(mu0=wbar5,Lambda=10*diag(7),kappa0=7,nu=7)
  dp5 <- Fit(DirichletProcessMvnormal(y=w5,g0Priors = g0priors5),100)
  
  g0priors6 <- list(mu0=wbar6,Lambda=10*diag(7),kappa0=7,nu=7)
  dp6 <- Fit(DirichletProcessMvnormal(y=w6,g0Priors = g0priors6),100)
  
  y_1<-y_2<-y_3<-y_4<-y_5<-y_6 <- list()
  acc1<-acc2 <-acc3<-acc4<-acc5<-acc6<-0
  cc <-0
  
  for(t in 3:mcmc){
    w1 <-cbind(yy1,h1[,1:2],xx1)
    wbar1 <- as.vector(apply(w1,2,mean))
    wcov1 <- cov(w1)
    g0priors1 <- list(mu0=wbar1,Lambda=10*diag(6),kappa0=6,nu=6)
    dp11 <- Fit(dp1,1)
    dp1 <- dp11
    
    w2 <-cbind(yy2,h2[,1:3],xx2)
    wbar2 <- as.vector(apply(w2,2,mean))
    wcov2 <- cov(w2)
    g0priors2 <- list(mu0=wbar2,Lambda=10*diag(7),kappa0=7,nu=7)
    dp12 <- Fit(dp2,1)
    dp2 <- dp12
    
    w3 <-cbind(yy3,h3[,1:3],xx3)
    wbar3 <- as.vector(apply(w3,2,mean))
    wcov3 <- cov(w3)
    g0priors3 <- list(mu0=wbar3,Lambda=10*diag(7),kappa0=7,nu=7)
    dp13 <- Fit(dp3,1)
    dp3 <- dp13
    
    w4 <-cbind(yy4,h4[,1:2],xx4)
    wbar4 <- as.vector(apply(w4,2,mean))
    wcov4 <- cov(w4)
    g0priors4 <- list(mu0=wbar4,Lambda=10*diag(6),kappa0=6,nu=6)
    dp14 <- Fit(dp4,1)
    dp4 <- dp14
    
    w5 <-cbind(yy5,h5[,c(1,2,4)],xx5)
    wbar5 <- as.vector(apply(w5,2,mean))
    wcov5 <- cov(w5)
    g0priors5 <- list(mu0=wbar5,Lambda=10*diag(7),kappa0=7,nu=7)
    dp15 <- Fit(dp5,1)
    dp5 <- dp15
    
    w6 <-cbind(yy6,h6[,c(1,2,4)],xx6)
    wbar6 <- as.vector(apply(w6,2,mean))
    wcov6 <- cov(w6)
    g0priors6 <- list(mu0=wbar6,Lambda=10*diag(7),kappa0=7,nu=7)
    dp16 <- Fit(dp6,1)
    dp6 <- dp16
    
    
    pi1 <- table(dp1$clusterLabels)
    Sigma1 <- list()
    Mu1 <- list()
    for(i in 1:length(pi1)){
      Mu1[[i]] <- dp1$clusterParameters$mu[,,i]
      Sigma1[[i]] <- dp1$clusterParameters$sig[,,i]
    }
    pi1 <- round(table(dp1$clusterLabels)/length(yy1),2)
    
    pi2 <- table(dp2$clusterLabels)
    Sigma2 <- list()
    Mu2 <- list()
    for(i in 1:length(pi2)){
      Mu2[[i]] <- dp2$clusterParameters$mu[,,i]
      Sigma2[[i]] <- dp2$clusterParameters$sig[,,i]
    }
    pi2 <- round(table(dp2$clusterLabels)/length(yy2),2)
    
    pi3 <- table(dp3$clusterLabels)
    Sigma3 <- list()
    Mu3 <- list()
    for(i in 1:length(pi3)){
      Mu3[[i]] <- dp3$clusterParameters$mu[,,i]
      Sigma3[[i]] <- dp3$clusterParameters$sig[,,i]
    }
    pi3 <- round(table(dp3$clusterLabels)/length(yy3),2)
    
    pi4 <- table(dp4$clusterLabels)
    Sigma4 <- list()
    Mu4 <- list()
    for(i in 1:length(pi4)){
      Mu4[[i]] <- dp4$clusterParameters$mu[,,i]
      Sigma4[[i]] <- dp4$clusterParameters$sig[,,i]
    }
    pi4 <- round(table(dp4$clusterLabels)/length(yy4),2)
    
    pi5 <- table(dp5$clusterLabels)
    Sigma5 <- list()
    Mu5 <- list()
    for(i in 1:length(pi5)){
      Mu5[[i]] <- dp5$clusterParameters$mu[,,i]
      Sigma5[[i]] <- dp5$clusterParameters$sig[,,i]
    }
    pi5 <- round(table(dp5$clusterLabels)/length(yy5),2)
    
    pi6 <- table(dp6$clusterLabels)
    Sigma6 <- list()
    Mu6 <- list()
    for(i in 1:length(pi6)){
      Mu6[[i]] <- dp6$clusterParameters$mu[,,i]
      Sigma6[[i]] <- dp6$clusterParameters$sig[,,i]
    }
    pi6 <- round(table(dp6$clusterLabels)/length(yy6),2)
    
    
    
    SEQ <- seq(1, mcmc, by=5)
    
    CondMVN<-function(mean, sigma, dependent.ind, given.ind, X.given, 
                      check.sigma = TRUE)
    {
      B <- sigma[dependent.ind, dependent.ind]
      C <- sigma[dependent.ind, given.ind, drop = FALSE]
      D <- sigma[given.ind, given.ind]
      CDinv <- C %*% solve(D)
      cMu <- mean[dependent.ind] + CDinv %*% t(X.given - 
                                                 matrix(mean[given.ind],
                                                        nrow=dim(X.given)[1], ncol=dim(X.given)[2],byrow=TRUE))
      list(condMean = cMu)
    }
    if(t %in% SEQ){
      cc <- cc + 1
      
      den<-sapply(1:length(pi1), function(i) dmnorm(cbind(h1[,1:2],xx1), 
                                                    Mu1[[i]][-c(1)], 
                                                    Sigma1[[i]][-c(1),
                                                                -c(1)])*pi1[i])/
        rowSums(sapply(1:length(pi1), function(i) dmnorm(cbind(h1[,1:2],xx1), 
                                                         Mu1[[i]][-c(1)], 
                                                         Sigma1[[i]][-c(1),-c(1)])*pi1[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_1[[cc]]<-rowSums(sapply(1:length(pi1), function(i) CondMVN(mean=Mu1[[i]] , sigma=Sigma1[[i]],
                                                                   dependent=c(1), given=c(2:6),X.given=
                                                                     cbind(h1[,1:2],xx1))$condMean)*(den1))
      
      den<-sapply(1:length(pi2), function(i) dmnorm(cbind(h2[,1:3],xx2), 
                                                    Mu2[[i]][-c(1)], 
                                                    Sigma2[[i]][-c(1),
                                                                -c(1)])*pi2[i])/
        rowSums(sapply(1:length(pi2), function(i) dmnorm(cbind(h2[,1:3],xx2), 
                                                         Mu2[[i]][-c(1)], 
                                                         Sigma2[[i]][-c(1),-c(1)])*pi2[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_2[[cc]]<-rowSums(sapply(1:length(pi2), function(i) CondMVN(mean=Mu2[[i]] , sigma=Sigma2[[i]],
                                                                   dependent=c(1), given=c(2:7),X.given=
                                                                     cbind(h2[,1:3],xx2))$condMean)*(den1))
      
      den<-sapply(1:length(pi3), function(i) dmnorm(cbind(h3[,1:3],xx3), 
                                                    Mu3[[i]][-c(1)], 
                                                    Sigma3[[i]][-c(1),
                                                                -c(1)])*pi3[i])/
        rowSums(sapply(1:length(pi3), function(i) dmnorm(cbind(h3[,1:3],xx3), 
                                                         Mu3[[i]][-c(1)], 
                                                         Sigma3[[i]][-c(1),-c(1)])*pi3[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_3[[cc]]<-rowSums(sapply(1:length(pi3), function(i) CondMVN(mean=Mu3[[i]] , sigma=Sigma3[[i]],
                                                                   dependent=c(1), given=c(2:7),X.given=
                                                                     cbind(h3[,1:3],xx3))$condMean)*(den1))
      
      den<-sapply(1:length(pi4), function(i) dmnorm(cbind(h4[,1:2],xx4), 
                                                    Mu4[[i]][-c(1)], 
                                                    Sigma4[[i]][-c(1),
                                                                -c(1)])*pi4[i])/
        rowSums(sapply(1:length(pi4), function(i) dmnorm(cbind(h4[,1:2],xx4), 
                                                         Mu4[[i]][-c(1)], 
                                                         Sigma4[[i]][-c(1),-c(1)])*pi4[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_4[[cc]]<-rowSums(sapply(1:length(pi4), function(i) CondMVN(mean=Mu4[[i]] , sigma=Sigma4[[i]],
                                                                   dependent=c(1), given=c(2:6),X.given=
                                                                     cbind(h4[,1:2],xx4))$condMean)*(den1))
      
      den<-sapply(1:length(pi5), function(i) dmnorm(cbind(h5[,c(1,2,4)],xx5), 
                                                    Mu5[[i]][-c(1)], 
                                                    Sigma5[[i]][-c(1),
                                                                -c(1)])*pi5[i])/
        rowSums(sapply(1:length(pi5), function(i) dmnorm(cbind(h5[,c(1,2,4)],xx5), 
                                                         Mu5[[i]][-c(1)], 
                                                         Sigma5[[i]][-c(1),-c(1)])*pi5[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_5[[cc]]<-rowSums(sapply(1:length(pi5), function(i) CondMVN(mean=Mu5[[i]] , sigma=Sigma5[[i]],
                                                                   dependent=c(1), given=c(2:7),X.given=
                                                                     cbind(h5[,c(1,2,4)],xx5))$condMean)*(den1))
      
      den<-sapply(1:length(pi6), function(i) dmnorm(cbind(h6[,c(1,2,4)],xx6), 
                                                    Mu6[[i]][-c(1)], 
                                                    Sigma6[[i]][-c(1),
                                                                -c(1)])*pi6[i])/
        rowSums(sapply(1:length(pi6), function(i) dmnorm(cbind(h6[,c(1,2,4)],xx6), 
                                                         Mu6[[i]][-c(1)], 
                                                         Sigma6[[i]][-c(1),-c(1)])*pi6[i]))
      den1<-ifelse(is.na(den), 0, den)
      y_6[[cc]]<-rowSums(sapply(1:length(pi6), function(i) CondMVN(mean=Mu6[[i]] , sigma=Sigma6[[i]],
                                                                   dependent=c(1), given=c(2:7),X.given=
                                                                     cbind(h6[,c(1,2,4)],xx6))$condMean)*(den1))
      
      
      
      ME[[cc]] <- list(DD1=DD1,DD2=DD2,DD3=DD3,DD4=DD4,DD5=DD5,DD6=DD6,Mu1=Mu1,
                       Mu2=Mu2,Mu3=Mu3,Mu4=Mu4,Mu5=Mu5,Mu6=Mu6,Sigma1=Sigma1,
                       Sigma2=Sigma2,Sigma3=Sigma3,Sigma4=Sigma4,Sigma5=Sigma5,
                       Sigma6=Sigma6,pi1=pi1,pi2=pi2,pi3=pi3,pi4=pi4,pi5=pi5,pi6
                       =pi6,para.C=para.C,h1=h1,h2=h2,h3=h3,h4=h4,h5=h5,h6=h6)
    }
    if(t > 200){
      cov.d11 <- cov(para.dd1[3:(t-1),(K+1):(2*K)])
    }
    
    h <- ifelse(is.nan(h), qnorm(0.1^5), h)
    para.dd1[t,] <- metropolis(h=h, Y=D[,1], X=rbind(xx1,xx2,xx3,xx4,xx5,xx6),
                               R=R, w_pre=para.dd1[t-1,1:K],
                               beta0_pre=para.dd1[t-1,ind1], 
                               beta_pre=para.dd1[t-1,ind3],
                               sigma_pre=para.dd1[t-1,ind2],
                               alpha_pre=para.dd1[t-1,ind4], 
                               alpha_beta0_pre=para.dd1[t-1,ind5], 
                               alpha_sigma_pre=para.dd1[t-1,ind6],
                               mu_beta0_pre=para.dd1[t-1,ind7], 
                               sigma_beta0_pre=para.dd1[t-1,ind8], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=1, 
                               cov1=cov.d11,cov2=cov2.d11,gamma=para.dd1[t-1,(ind9-n)], 
                               Z=para.dd1[t-1,(ind9-n+1):ind9])
    
    Zm1 <- para.dd1[t,(ind9-n+1):ind9][1:(n1+n2+n3)]
    h[1:(n1+n2+n3),1] <- qnorm(pmin(1-0.1^5,pmax(0.1^5, ptruncnorm(D[1:(n1+n2+n3),1],a=0,b=1,
                                                                   para.dd1[t,ind1][Zm1]+
                                                                     rbind(xx1,xx2,xx3)%*%as.matrix(para.dd1[t,ind3]),
                                                                   sqrt(para.dd1[t,ind2])[Zm1]))))
    
    
    
    if(t > 200){
      cov.d12 <- cov(para.dd2[3:(t-1),(K+1):(2*K)])
    }
    h <- ifelse(is.nan(h), qnorm(0.1^5), h)
    para.dd2[t,] <- metropolis(h=h, Y=D[,2], 
                               X=rbind(xx1,xx2,xx3,xx4,xx5,xx6),R=R, 
                               w_pre=para.dd2[t-1,1:K],
                               beta0_pre=para.dd2[t-1,ind1], 
                               beta_pre=para.dd2[t-1,ind3],
                               sigma_pre=para.dd2[t-1,ind2],
                               alpha_pre=para.dd2[t-1,ind4], 
                               alpha_beta0_pre=para.dd2[t-1,ind5], 
                               alpha_sigma_pre=para.dd2[t-1,ind6],
                               mu_beta0_pre=para.dd2[t-1,ind7], 
                               sigma_beta0_pre=para.dd2[t-1,ind8], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=2, 
                               cov1=cov.d12,cov2=cov2.d12,gamma=para.dd2[t-1,(ind9-n)], 
                               Z=para.dd2[t-1,(ind9-n+1):ind9])
    
    Zm2 <- para.dd2[t,(ind9-n+1):ind9][(n1+n2+n3+1):n]
    h[(n1+n2+n3+1):n,2] <- qnorm(pmin(1-0.1^5,pmax(0.1^5, 
                                                   ptruncnorm(D[(n1+n2+n3+1):n,2],a=0,b=1,
                                                              para.dd2[t,ind1][Zm2]+
                                                                rbind(xx4,xx5,xx6)%*%as.matrix(para.dd2[t,ind3]),
                                                              sqrt(para.dd2[t,ind2])[Zm2]))))
    
    
    
    if(t > 200){
      cov.d13 <- cov(para.dd3[3:(t-1),(K+1):(2*K)])
    }
    h <- ifelse(is.nan(h), qnorm(0.1^5), h)
    para.dd3[t,] <- metropolis(h=h, Y=D[,3], X=cbind(rbind(xx1,xx2,xx3,xx4,xx5,xx6),D[,1]),R=R, 
                               w_pre=para.dd3[t-1,1:K],
                               beta0_pre=para.dd3[t-1,ind1], 
                               beta_pre=para.dd3[t-1,ind13],
                               sigma_pre=para.dd3[t-1,ind12],
                               alpha_pre=para.dd3[t-1,ind14], 
                               alpha_beta0_pre=para.dd3[t-1,ind15], 
                               alpha_sigma_pre=para.dd3[t-1,ind16],
                               mu_beta0_pre=para.dd3[t-1,ind17], 
                               sigma_beta0_pre=para.dd3[t-1,ind18], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=3, 
                               cov1=cov.d13,cov2=cov2.d13,gamma=para.dd3[t-1,(ind19-n)], 
                               Z=para.dd3[t-1,(ind19-n+1):ind19])
    
    Zm3 <- para.dd3[t,(ind19-n+1):ind19][(n1+1):(n1+n2)]
    h[(n1+1):(n1+n2),3] <- qnorm(pmin(1-0.1^5,pmax(0.1^5, 
                                                   ptruncnorm(D[(n1+1):(n1+n2),3],a=0,b=1,
                                                              para.dd3[t,ind1][Zm3]+
                                                                cbind(xx2,D[(n1+1):(n1+n2),1])%*%as.matrix(para.dd3[t,ind13]),
                                                              sqrt(para.dd3[t,ind12])[Zm3]))))
    if(t > 200){
      cov.d14 <- cov(para.dd4[3:(t-1),(K+1):(2*K)])
    }
    h <- ifelse(is.nan(h), qnorm(0.1^5), h)
    para.dd4[t,] <- metropolis(h=h, Y=D[,4], 
                               X=cbind(rbind(xx1,xx2,xx3,xx4,xx5,xx6),D[,2]),R=R, 
                               w_pre=para.dd4[t-1,1:K],
                               beta0_pre=para.dd4[t-1,ind1], 
                               beta_pre=para.dd4[t-1,ind13],
                               sigma_pre=para.dd4[t-1,ind12],
                               alpha_pre=para.dd4[t-1,ind14], 
                               alpha_beta0_pre=para.dd4[t-1,ind15], 
                               alpha_sigma_pre=para.dd4[t-1,ind16],
                               mu_beta0_pre=para.dd4[t-1,ind17], 
                               sigma_beta0_pre=para.dd4[t-1,ind18], K=K, 
                               eps=0.1,del1=15, del2=10, del3=20, index=4, 
                               cov1=cov.d14,cov2=cov2.d14,gamma=para.dd4[t-1,(ind19-n)], 
                               Z=para.dd4[t-1,(ind19-n+1):ind19])
    
    Zm4 <- para.dd4[t,(ind19-n+1):ind19][(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]
    h[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),4] <- qnorm(pmin(1-0.1^5,pmax(0.1^5, 
                                                                     ptruncnorm(D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),4],a=0,b=1,
                                                                                para.dd4[t,ind1][Zm4]+
                                                                                  cbind(xx5,D[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),2])%*%as.matrix(para.dd4[t,ind13]),
                                                                                sqrt(para.dd4[t,ind12])[Zm4]))))
    
    para.C[t,] <- metropolisC1(h=h, rho=para.C[t-1,7:12])
    prop1 <- para.C[t,7:12]
    
    R <- matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],
                  1,prop1[6],prop1[3],prop1[5],prop1[6],1),4,4,byrow=TRUE)
    
    
    
    h1_prop <- h1
    DD1_prop <- DD1
    
    MNORM <- function(n = 1, mean = rep(0, d), varcov){
      sqrt.varcov <-  chol(varcov)
      d <- ncol(sqrt.varcov)
      return(drop(mean + t(matrix(rnorm(n * d), d, n)) %*% sqrt.varcov))
    }
    
    DD1_prop[,2] <- rtruncnorm(n1,DD1[,2],sd(DD1[,2])/5,a=0,b=1)
    clus1 <- para.dd2[t,(ind9-n1+1):ind9][1:n1]
    h1_prop[,2] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD1_prop[,2],a=0,b=1,
                                                     para.dd2[t,ind1][clus1]+
                                                       xx1%*%as.matrix(para.dd2[t,ind3]),
                                                     sqrt(para.dd2[t,ind2])[clus1]))))
    
    #clus2 <- para.dd3[t,(ind9-n+1):ind9][1:n1]
    #h1_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #                                     pnorm(DD1_prop[,3],para.dd3[t,ind1][clus2],
    #                                           sqrt(para.dd3[t,ind2])[clus2]))))
    
    
    rat1 <- dmnorm(h1_prop[,1:2],rep(0,2), R[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD1_prop[,2],a=0,b=1,para.dd2[t,ind1][clus1]+
                                    xx1%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus1])))+
      #log(pmax(0.1^100,dnorm(DD1_prop[,3])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi1), function(i) dmnorm(cbind(yy1,h1_prop[,1:2],xx1), Mu1[[i]], Sigma1[[i]])*pi1[i]))/
                 rowSums(sapply(1:length(pi1), function(i) dmnorm(cbind(h1_prop[,1:2],xx1), Mu1[[i]][-1], Sigma1[[i]][-c(1),-c(1)])*pi1[i]))))+
      sapply(1:n1, function(c) log(dtruncnorm(DD1[c,2],DD1_prop[c,2], 
                                              sd(DD1[,2])/5,a=0,b=1)))
    
    rat2 <- dmnorm(h1[,1:2],rep(0,2), R[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD1[,2],a=0,b=1,para.dd2[t,ind1][clus1]+
                                    xx1%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus1])))+
      #log(pmax(0.1^100,dnorm(DD1[,3])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi1), function(i) dmnorm(cbind(yy1,h1[,1:2],xx1), Mu1[[i]], Sigma1[[i]])*pi1[i]))/
                 rowSums(sapply(1:length(pi1), function(i) dmnorm(cbind(h1[,1:2],xx1), Mu1[[i]][-c(1)], Sigma1[[i]][-c(1),-c(1)])*pi1[i]))))+
      sapply(1:n1, function(c) log(dtruncnorm(DD1_prop[c,2],DD1[c,2], 
                                              sd(DD1[,2])/5,a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc1 <- rep(0,n1)
    
    for(ii in 1:n1){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD1_prop[ii,2] <- DD1[ii,2]
        h1_prop[ii,2] <- h1[ii,2]
      }else{
        DD1[ii,2] <- DD1_prop[ii,2]
        h1[ii,2] <- h1_prop[ii,2]
        acc1[ii] <- 1
      }
    }
    
    h2_prop <- h2
    DD2_prop <- DD2
    
    DD2_prop[,2] <- rtruncnorm(n2,DD2[,2],sd(DD2[,2])/5,a=0,b=1)
    clus3 <- para.dd2[t,(ind9-n+1):ind9][(n1+1):(n1+n2)]
    #clus13 <- para.dd3[t,(ind19-n+1):ind19][(n1+1):(n1+n2)]
    h2_prop[,2] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD2_prop[,2],a=0,b=1,
                                                     para.dd2[t,ind1][clus3]+
                                                       xx2%*%as.matrix(para.dd2[t,ind3]),
                                                     sqrt(para.dd2[t,ind2])[clus3]))))
    
    
    rat1 <- dmnorm(h2_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2_prop[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus3]+
                                    xx2%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus3])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi2), function(i) dmnorm(cbind(yy2,h2_prop[,1:3],xx2), Mu2[[i]], Sigma2[[i]])*pi2[i]))/
                 rowSums(sapply(1:length(pi2), function(i) dmnorm(cbind(h2_prop[,1:3],xx2), Mu2[[i]][-1], Sigma2[[i]][-c(1),-c(1)])*pi2[i]))))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2[c,2],DD2_prop[c,2], 
                                              sd(DD2[,2])/5,a=0,b=1)))
    
    rat2 <- dmnorm(h2,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD2[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus3]+
                                    xx2%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus3])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi2), function(i) dmnorm(cbind(yy2,h2[,1:3],xx2), Mu2[[i]], Sigma2[[i]])*pi2[i]))/
                 rowSums(sapply(1:length(pi2), function(i) dmnorm(cbind(h2[,1:3],xx2), Mu2[[i]][-c(1)], Sigma2[[i]][-c(1),-c(1)])*pi2[i]))))+
      sapply(1:n2, function(c) log(dtruncnorm(DD2_prop[c,2],DD2[c,2], 
                                              sd(DD2[,2])/5,a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc2 <- rep(0,n2)
    
    for(ii in 1:n2){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD2_prop[ii,2] <- DD2[ii,2]
        h2_prop[ii,2] <- h2[ii,2]
      }else{
        DD2[ii,2] <- DD2_prop[ii,2]
        h2[ii,2] <- h2_prop[ii,2]
        acc2[ii] <- 1
      }
    }
    
    h3_prop <- h3
    DD3_prop <- DD3
    
    DD3_prop[,2:3] <- MNORM(n3,DD3[,2:3],cov(DD3[,2:3])/5)
    clus4 <- para.dd2[t,(ind9-n+1):ind9][(n1+n2+1):(n1+n2+n3)]
    h3_prop[,2] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD3_prop[,2],a=0,b=1,
                                                     para.dd2[t,ind1][clus4]+
                                                       xx3%*%as.matrix(para.dd2[t,ind3]),
                                                     sqrt(para.dd2[t,ind2])[clus4]))))
    
    clus5 <- para.dd3[t,(ind19-n+1):ind19][(n1+n2+1):(n1+n2+n3)]
    h3_prop[,3] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD3_prop[,3],a=0,b=1,
                                                     para.dd3[t,ind1][clus5]+
                                                       cbind(xx3,DD3[,1])%*%as.matrix(para.dd3[t,ind13]),
                                                     sqrt(para.dd3[t,ind12])[clus5]))))
    #h3_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #pnorm(DD3_prop[,3],para.dd3[t,ind1][clus5],
    #sqrt(para.dd3[t,ind2])[clus5]))))
    
    rat1 <- dmnorm(h3_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3_prop[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus4]+
                                    xx3%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus4])))+
      log(pmax(0.1^100,dtruncnorm(DD3_prop[,3],a=0,b=1,
                                  para.dd3[t,ind1][clus5]+
                                    cbind(xx3,DD3[,1])%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus5])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi3), function(i) dmnorm(cbind(yy3,h3_prop[,1:3],xx3), Mu3[[i]], Sigma3[[i]])*pi3[i]))/
                 rowSums(sapply(1:length(pi3), function(i) dmnorm(cbind(h3_prop[,1:3],xx3), Mu3[[i]][-1], Sigma3[[i]][-c(1),-c(1)])*pi3[i]))))+
      sapply(1:n3, function(c) dmnorm(DD3[c,2:3],DD3_prop[c,2:3], 
                                      cov(DD3[,2:3])/5, log=TRUE))
    
    rat2 <- dmnorm(h3,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD3[,2],a=0,b=1,
                                  para.dd2[t,ind1][clus4]+
                                    xx3%*%as.matrix(para.dd2[t,ind3]),
                                  sqrt(para.dd2[t,ind2])[clus4])))+
      log(pmax(0.1^100,dtruncnorm(DD3[,3],a=0,b=1,
                                  para.dd3[t,ind1][clus5]+
                                    cbind(xx3,DD3[,1])%*%as.matrix(para.dd3[t,ind13]),
                                  sqrt(para.dd3[t,ind12])[clus5])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi3), function(i) dmnorm(cbind(yy3,h3[,1:3],xx3), Mu3[[i]], Sigma3[[i]])*pi3[i]))/
                 rowSums(sapply(1:length(pi3), function(i) dmnorm(cbind(h3[,1:3],xx3), Mu3[[i]][-c(1)], Sigma3[[i]][-c(1),-c(1)])*pi3[i]))))+
      sapply(1:n3, function(c) dmnorm(DD3_prop[c,2:3],DD3[c,2:3], 
                                      cov(DD3[,2:3])/5, log=TRUE))
    
    rat <- rat1 - rat2
    
    acc3 <- rep(0,n3)
    
    for(ii in 1:n3){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD3_prop[ii,2:3] <- DD3[ii,2:3]
        h3_prop[ii,2:3] <- h3[ii,2:3]
      }else{
        DD3[ii,2:3] <- DD3_prop[ii,2:3]
        h3[ii,2:3] <- h3_prop[ii,2:3]
        acc3[ii] <- 1
      }
    }
    h4_prop<-h4; DD4_prop<-DD4
    
    DD4_prop[,1] <- rtruncnorm(n4,DD4[,1],sd(DD4[,1])/5,a=0,b=1)
    clus6 <- para.dd1[t,(ind9-n+1):ind9][(1+n1+n2+n3):(n1+n2+n3+n4)]
    h4_prop[,1] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD4_prop[,1],a=0,b=1,
                                                     para.dd1[t,ind1][clus6]+
                                                       xx4%*%as.matrix(para.dd1[t,ind3]),
                                                     sqrt(para.dd1[t,ind2])[clus6]))))
    
    #clus7 <- para.dd3[t,(ind9-n+1):ind9][1:n1][(1+n1+n2+n3):(n1+n2+n3+n4)]
    #h4_prop[,3] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #                                   pnorm(DD4_prop[,3],para.dd3[t,ind1][clus7],
    #                                        sqrt(para.dd3[t,ind2])[clus7]))))
    
    
    rat1 <- dmnorm(h4_prop[,1:2],c(0,0), R[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD4_prop[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus6]+
                                    xx4%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus6])))+
      #log(pmax(0.1^100,dnorm(DD4_prop[,3])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi4), function(i) dmnorm(cbind(yy4,h4_prop[,1:2],xx4), Mu4[[i]], Sigma4[[i]])*pi4[i]))/
                 rowSums(sapply(1:length(pi4), function(i) dmnorm(cbind(h4_prop[,1:2],xx4), Mu4[[i]][-1], Sigma4[[i]][-c(1),-c(1)])*pi4[i]))))+
      sapply(1:n4, function(c) log(dtruncnorm(DD4[c,1],DD4_prop[c,1], 
                                              sd(DD4[,1])/5, a=0,b=1)))
    
    rat2 <- dmnorm(h4[,1:2],c(0,0), R[1:2,1:2], log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD4[,1],a=0,b=1,para.dd1[t,ind1][clus6]+
                                    xx4%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus6])))+
      #log(pmax(0.1^100,dnorm(DD4[,3])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi4), function(i) dmnorm(cbind(yy4,h4[,1:2],xx4), Mu4[[i]], Sigma4[[i]])*pi4[i]))/
                 rowSums(sapply(1:length(pi4), function(i) dmnorm(cbind(h4[,1:2],xx4), Mu4[[i]][-c(1)], Sigma4[[i]][-c(1),-c(1)])*pi4[i]))))+
      sapply(1:n4, function(c) log(dtruncnorm(DD4_prop[c,1],DD4[c,1], 
                                              sd(DD4[,1])/5,a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc4 <- rep(0,n4)
    
    for(ii in 1:n4){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD4_prop[ii,1] <- DD4[ii,1]
        h4_prop[ii,1] <- h4[ii,1]
      }else{
        DD4[ii,1] <- DD4_prop[ii,1]
        h4[ii,1] <- h4_prop[ii,1]
        acc4[ii] <- 1
      }
    }
    
    h5_prop <- h5
    DD5_prop <- DD5
    
    DD5_prop[,1] <- rtruncnorm(n5,DD5[,1],sd(DD5[,1])/5,a=0,b=1)
    clus8 <- para.dd1[t,(ind9-n+1):ind9][(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]
    #clus18 <- para.dd4[t,(ind19-n+1):ind19][(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5)]
    h5_prop[,1] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD5_prop[,1],a=0,b=1,para.dd1[t,ind1][clus8]+
                                                       xx5%*%as.matrix(para.dd1[t,ind3]),
                                                     sqrt(para.dd1[t,ind2])[clus8]))))
    # h5_prop[,4] <-qnorm(pmin(1-0.1^3,pmax(0.1^100,
    #                                 pnorm(DD5[,4],para.dd4[t,ind1][clus18]+
    #                                        cbind(xx5,DD5[,2])%*%as.matrix(para.dd4[t,ind13]),
    #                                     sqrt(para.dd4[t,ind12])[clus18]))))
    
    
    rat1 <- dmnorm(h5_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5_prop[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus8]+
                                    xx5%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus8])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi5), function(i) dmnorm(cbind(yy5,h5_prop[,c(1,2,4)],xx5), Mu5[[i]], Sigma5[[i]])*pi5[i]))/
                 rowSums(sapply(1:length(pi5), function(i) dmnorm(cbind(h5_prop[,c(1,2,4)],xx5), Mu5[[i]][-1], Sigma5[[i]][-c(1),-c(1)])*pi5[i]))))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5[c,1],DD5_prop[c,1], 
                                              sd(DD5[,1]),a=0,b=1)))
    
    rat2 <- dmnorm(h5,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD5[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus8]+
                                    xx5%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus8])))+
      log(pmax(0.1^200,rowSums(sapply(1:length(pi5), function(i) dmnorm(cbind(yy5,h5[,c(1,2,4)],xx5), Mu5[[i]], Sigma5[[i]])*pi5[i]))/
                 rowSums(sapply(1:length(pi5), function(i) dmnorm(cbind(h5[,c(1,2,4)],xx5), Mu5[[i]][-c(1)], Sigma5[[i]][-c(1),-c(1)])*pi5[i]))))+
      sapply(1:n5, function(c) log(dtruncnorm(DD5_prop[c,1],DD5[c,1], 
                                              sd(DD5[,1]),a=0,b=1)))
    
    rat <- rat1 - rat2
    
    acc5 <- rep(0,n5)
    
    for(ii in 1:n5){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD5_prop[ii,1] <- DD5[ii,1]
        h5_prop[ii,1] <- h5[ii,1]
      }else{
        DD5[ii,1] <- DD5_prop[ii,1]
        h5[ii,1] <- h5_prop[ii,1]
        acc5[ii] <- 1
      }
    }
    
    h6_prop <- h6
    DD6_prop <- DD6
    
    DD6_prop[,c(1,4)] <- MNORM(n6,DD6[,c(1,4)],cov(DD6[,c(1,4)])/5)
    clus9 <- para.dd1[t,(ind9-n+1):ind9][(n1+n2+n3+n4+n5+1):n]
    h6_prop[,1] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD6_prop[,1],a=0,b=1,
                                                     para.dd1[t,ind1][clus9]+
                                                       xx6%*%as.matrix(para.dd1[t,ind3]),
                                                     sqrt(para.dd1[t,ind2])[clus9]))))
    
    clus10 <- para.dd4[t,(ind19-n+1):ind19][(n1+n2+n3+n4+n5+1):n]
    h6_prop[,4] <-qnorm(pmin(1-0.1^5,pmax(0.1^100,
                                          ptruncnorm(DD6_prop[,4],a=0,b=1,
                                                     para.dd4[t,ind1][clus10]+
                                                       cbind(xx6,DD6[,2])%*%as.matrix(para.dd4[t,ind13]),
                                                     sqrt(para.dd4[t,ind12])[clus10]))))
    
    rat1 <- dmnorm(h6_prop,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6_prop[,1],a=0,b=1,
                                  para.dd1[t,ind1][clus9]+
                                    xx6%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus9])))+
      log(pmax(0.1^100,dtruncnorm(DD6_prop[,4],a=0,b=1,
                                  para.dd4[t,ind1][clus10]+
                                    cbind(xx6,DD6[,2])%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus10])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi6), function(i) dmnorm(cbind(yy6,h6_prop[,c(1,2,4)],xx6), Mu6[[i]], Sigma6[[i]])*pi6[i]))/
                 rowSums(sapply(1:length(pi6), function(i) dmnorm(cbind(h6_prop[,c(1,2,4)],xx6), Mu6[[i]][-1], Sigma6[[i]][-c(1),-c(1)])*pi6[i]))))+
      sapply(1:n6, function(c) dmnorm(DD6[c,c(1,4)],DD6_prop[c,c(1,4)], 
                                      cov(DD6[,c(1,4)])/5, log=TRUE))
    
    rat2 <- dmnorm(h6,rep(0,4), R, log=TRUE)+
      log(pmax(0.1^100,dtruncnorm(DD6[,1],a=0,b=1,para.dd1[t,ind1][clus9]+
                                    xx6%*%as.matrix(para.dd1[t,ind3]),
                                  sqrt(para.dd1[t,ind2])[clus9])))+
      log(pmax(0.1^100,dtruncnorm(DD6[,4],a=0,b=1,para.dd4[t,ind1][clus10]+
                                    cbind(xx6,DD6[,2])%*%as.matrix(para.dd4[t,ind13]),
                                  sqrt(para.dd4[t,ind12])[clus10])))+
      log(pmax(0.1^300,rowSums(sapply(1:length(pi6), function(i) dmnorm(cbind(yy6,h6[,c(1,2,4)],xx6), Mu6[[i]], Sigma6[[i]])*pi6[i]))/
                 rowSums(sapply(1:length(pi6), function(i) dmnorm(cbind(h6[,c(1,2,4)],xx6), Mu6[[i]][-c(1)], Sigma6[[i]][-c(1),-c(1)])*pi6[i]))))+
      sapply(1:n6, function(c) dmnorm(DD6_prop[c,c(1,4)],DD6[c,c(1,4)], 
                                      cov(DD6[,c(1,4)])/5, log=TRUE))
    
    rat <- rat1 - rat2
    
    acc6 <- rep(0,n6)
    
    for(ii in 1:n6){    
      if(log(runif(1))>rat[ii] | is.na(rat[ii])){
        DD6_prop[ii,c(1,4)] <- DD6[ii,c(1,4)]
        h6_prop[ii,c(1,4)] <- h6[ii,c(1,4)]
      }else{
        DD6[ii,c(1,4)] <- DD6_prop[ii,c(1,4)]
        h6[ii,c(1,4)] <- h6_prop[ii,c(1,4)]
        acc6[ii] <- 1
      }
    }
  }
  return(list(y_1=y_1,y_2=y_2,y_3=y_3,y_4=y_4,y_5=y_5,y_6=y_6,ME=ME,dt=data))
}
