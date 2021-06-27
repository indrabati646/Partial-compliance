# Obtaining the EDTRs for Example 3

source("bnp_model_Example3.R")
source("Simulate_Example3.R")
source("res_prob_Example3.R")

edtr_fun <- function(data){
  results_list <- bnp_mod(data)
  data <- matrix(unlist(data),nrow=nrow(data),ncol=ncol(data))
  colnames(data) <- c("a1","a2","s","D1","D2","D3","D4","y","y1","X1","X2","X3")
  
  
  PS_75 <- c(0.75,0.75,0.75,0.75,median(data[,"X1"]),median(data[,"X2"]))
  
  
  
  PS_50 <- c(0.50,0.50,0.50,0.50,median(data[,"X1"]),median(data[,"X2"]))
  
  
  PS_25 <- c(0.25,0.25,0.25,0.25,median(data[,"X1"]),median(data[,"X2"]))
  PS_00 <- c(1,1,1,1,median(data[,"X1"]),median(data[,"X2"]))
  
  Mu_seq1 <- unlist(results_list$ME[[1]]$Mu1)
  Mu_seq2 <- unlist(results_list$ME[[1]]$Mu2)
  Mu_seq3 <- unlist(results_list$ME[[1]]$Mu3)
  Mu_seq4 <- unlist(results_list$ME[[1]]$Mu4)
  Mu_seq5 <- unlist(results_list$ME[[1]]$Mu5)
  Mu_seq6 <- unlist(results_list$ME[[1]]$Mu6)
  
  for(i in 2:length(results_list$ME)){
    Mu_seq1 <- cbind(Mu_seq1,apply(simplify2array(results_list$ME[[i]]$Mu1),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi1)))
    Mu_seq2 <- cbind(Mu_seq2,apply(simplify2array(results_list$ME[[i]]$Mu2),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi2)))
    Mu_seq3 <- cbind(Mu_seq3,apply(simplify2array(results_list$ME[[i]]$Mu3),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi3)))
    Mu_seq4 <- cbind(Mu_seq4,apply(simplify2array(results_list$ME[[i]]$Mu4),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi4)))
    Mu_seq5 <- cbind(Mu_seq5,apply(simplify2array(results_list$ME[[i]]$Mu5),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi5)))
    Mu_seq6 <- cbind(Mu_seq6,apply(simplify2array(results_list$ME[[i]]$Mu6),
                                   c(1,2),function(x) sum(x*results_list$ME[[i]]$pi6)))
    
  }
  
  
  Mu_seq1 <- rowMeans(Mu_seq1); Mu_seq2 <- rowMeans(Mu_seq2)
  Mu_seq3 <- rowMeans(Mu_seq3); Mu_seq4 <- rowMeans(Mu_seq4)
  Mu_seq5 <- rowMeans(Mu_seq5); Mu_seq6 <- rowMeans(Mu_seq6)
  
  Sig_seq1 <- matrix(unlist(results_list$ME[[1]]$Sigma1),ncol=5,
                     byrow=T)/length(results_list$ME)
  Sig_seq2 <- matrix(unlist(results_list$ME[[1]]$Sigma2),ncol=6,
                     byrow=T)/length(results_list$ME)
  Sig_seq3 <- matrix(unlist(results_list$ME[[1]]$Sigma3),ncol=6,
                     byrow=T)/length(results_list$ME)
  Sig_seq4 <-matrix(unlist(results_list$ME[[1]]$Sigma4),ncol=5,
                    byrow=T)/length(results_list$ME)
  Sig_seq5 <- matrix(unlist(results_list$ME[[1]]$Sigma5),ncol=6,
                     byrow=T)/length(results_list$ME)
  Sig_seq6 <- matrix(unlist(results_list$ME[[1]]$Sigma6),ncol=6,
                     byrow=T)/length(results_list$ME)
  
  for(i in 2:length(results_list$ME)){
    Sig_seq1 <- Sig_seq1+apply(simplify2array(results_list$ME[[i]]$Sigma1),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi1))/length(results_list$ME)
    Sig_seq2 <- Sig_seq2+apply(simplify2array(results_list$ME[[i]]$Sigma2),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi2))/length(results_list$ME)
    Sig_seq3 <- Sig_seq3+apply(simplify2array(results_list$ME[[i]]$Sigma3),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi3))/length(results_list$ME)
    Sig_seq4 <- Sig_seq4+apply(simplify2array(results_list$ME[[i]]$Sigma4),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi4))/length(results_list$ME)
    Sig_seq5 <- Sig_seq5+apply(simplify2array(results_list$ME[[i]]$Sigma5),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi5))/length(results_list$ME)
    Sig_seq6 <- Sig_seq6+apply(simplify2array(results_list$ME[[i]]$Sigma6),
                               c(1,2),function(x) sum(x*results_list$ME[[i]]$pi6))/length(results_list$ME)
    
  }
  res_prob <- response_prob(sim_1[[1]],results_list)
  a1ind1 <- which(colnames(res_prob$samples)=="alpha1[1]")
  a1ind2 <- which(colnames(res_prob$samples)=="alpha1[5]")
  b11ind1 <- which(colnames(res_prob$samples)=="beta11[1]")
  b11ind2 <- which(colnames(res_prob$samples)=="beta11[5]")
  b12ind1 <- which(colnames(res_prob$samples)=="beta12[1]")
  b12ind2 <- which(colnames(res_prob$samples)=="beta12[5]")
  b13ind1 <- which(colnames(res_prob$samples)=="beta3[1]")
  b13ind2 <- which(colnames(res_prob$samples)=="beta3[5]")
  b14ind1 <- which(colnames(res_prob$samples)=="beta4[1]")
  b14ind2 <- which(colnames(res_prob$samples)=="beta4[5]")

  
  w1ind1 <- which(colnames(res_prob$samples)=="w[1]")
  w1ind2 <- which(colnames(res_prob$samples)=="w[5]")
  
  
  a2ind1 <- which(colnames(res_prob$samples1)=="alpha1[1]")
  a2ind2 <- which(colnames(res_prob$samples1)=="alpha1[5]")
  b21ind1 <- which(colnames(res_prob$samples1)=="beta21[1]")
  b21ind2 <- which(colnames(res_prob$samples1)=="beta21[5]")
  b22ind1 <- which(colnames(res_prob$samples1)=="beta22[1]")
  b22ind2 <- which(colnames(res_prob$samples1)=="beta22[5]")
  b23ind1 <- which(colnames(res_prob$samples1)=="beta3[1]")
  b23ind2 <- which(colnames(res_prob$samples1)=="beta3[5]")
  b24ind1 <- which(colnames(res_prob$samples1)=="beta4[1]")
  b24ind2 <- which(colnames(res_prob$samples1)=="beta4[5]")
  
  w2ind1 <- which(colnames(res_prob$samples1)=="w[1]")
  w2ind2 <- which(colnames(res_prob$samples1)=="w[5]")
  
  
  
  
  y_seq1 <- condMVN(Mu_seq1,Sig_seq1,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_75[-c(3,4)])
  y_seq2 <- condMVN(Mu_seq2,Sig_seq2,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_75[-4])
  y_seq3 <- condMVN(Mu_seq3,Sig_seq3,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_75[-4])
  y_seq4 <- condMVN(Mu_seq4,Sig_seq4,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_75[-c(3,4)])
  y_seq5 <- condMVN(Mu_seq5,Sig_seq5,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_75[-3])
  y_seq6 <- condMVN(Mu_seq6,Sig_seq6,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_75[-3])
  
  out1 <- rnorm(length(results_list$ME),y_seq1$condMean,sqrt(y_seq1$condVar)/
                  length(results_list$ME[[1]]$DD1[,1]))
  out2 <- rnorm(length(results_list$ME),y_seq2$condMean,sqrt(y_seq2$condVar)/
                  length(results_list$ME[[1]]$DD2[,1]))
  out3 <- rnorm(length(results_list$ME),y_seq3$condMean,sqrt(y_seq3$condVar)/
                  length(results_list$ME[[1]]$DD3[,1]))
  out4 <- rnorm(length(results_list$ME),y_seq4$condMean,sqrt(y_seq4$condVar)/
                  length(results_list$ME[[1]]$DD4[,1]))
  out5 <- rnorm(length(results_list$ME),y_seq5$condMean,sqrt(y_seq5$condVar)/
                  length(results_list$ME[[1]]$DD5[,1]))
  out6 <- rnorm(length(results_list$ME),y_seq6$condMean,sqrt(y_seq6$condVar)/
                  length(results_list$ME[[1]]$DD6[,1]))
  
  prob1 <- matrix(NA,length(results_list$ME),5)
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob1[i,j] <- res_prob$samples[i,w1ind1+j-1]*(res_prob$samples[i,a1ind1+j-1]+
                                                   res_prob$samples[i,b11ind1+j-1]*PS_75[1]+
                                                   res_prob$samples[i,b12ind1+j-1]*PS_75[2]+
                                                   res_prob$samples[i,b13ind1+j-1]*PS_75[5]+
                                                   res_prob$samples[i,b14ind1+j-1]*PS_75[6])
    }
  }
  
  g1_75 <- sapply(rowMeans(prob1),expit)
  
  prob2 <- matrix(NA,length(results_list$ME),5)
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob2[i,j] <- res_prob$samples1[i,w2ind1+j-1]*(res_prob$samples1[i,a2ind1+j-1]+
                                                     res_prob$samples1[i,b21ind1+j-1]*PS_75[1]+
                                                     res_prob$samples1[i,b22ind1+j-1]*PS_75[2]+
                                                     res_prob$samples1[i,b23ind1+j-1]*PS_75[5]+
                                                     res_prob$samples1[i,b24ind1+j-1]*PS_75[6])
    }
  }
  
  g2_75 <- sapply(rowMeans(prob2),expit)
  
  theta_1_75 <- out1*g1_75+out2*(1-g1_75)
  theta_2_75 <- out1*g1_75+out3*(1-g1_75)
  theta_3_75 <- out4*g2_75+out5*(1-g2_75)
  theta_4_75 <- out4*g2_75+out6*(1-g2_75)
  
  #Compute means and standard errors for each embedded DTR
  theta_mean_75 <- cbind(theta_1_75,theta_2_75,theta_3_75,theta_4_75)
  colnames(theta_mean_75) <- c("EDTR 1","EDTR 2", 'EDTR 3',"EDTR 4")
  EDTR_75_mean <- apply(theta_mean_75,2,mean)
  EDTR_75_sd <- apply(theta_mean_75,2,sd)
  diff1<- -theta_mean_75-max(-EDTR_75_mean)
  quan1<- apply(diff1,2,function(x) quantile(x,1-0.05/3))
  bst_75 <- which(quan1<0)
  
  
  y_seq1 <- condMVN(Mu_seq1,Sig_seq1,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_50[-c(3,4)])
  y_seq2 <- condMVN(Mu_seq2,Sig_seq2,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_50[-4])
  y_seq3 <- condMVN(Mu_seq3,Sig_seq3,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_50[-4])
  y_seq4 <- condMVN(Mu_seq4,Sig_seq4,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_50[-c(3,4)])
  y_seq5 <- condMVN(Mu_seq5,Sig_seq5,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_50[-3])
  y_seq6 <- condMVN(Mu_seq6,Sig_seq6,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_50[-3])
  
  out1 <- rnorm(length(results_list$ME),y_seq1$condMean,sqrt(y_seq1$condVar)/
                  length(results_list$ME[[1]]$DD1[,1]))
  out2 <- rnorm(length(results_list$ME),y_seq2$condMean,sqrt(y_seq2$condVar)/
                  length(results_list$ME[[1]]$DD2[,1]))
  out3 <- rnorm(length(results_list$ME),y_seq3$condMean,sqrt(y_seq3$condVar)/
                  length(results_list$ME[[1]]$DD3[,1]))
  out4 <- rnorm(length(results_list$ME),y_seq4$condMean,sqrt(y_seq4$condVar)/
                  length(results_list$ME[[1]]$DD4[,1]))
  out5 <- rnorm(length(results_list$ME),y_seq5$condMean,sqrt(y_seq5$condVar)/
                  length(results_list$ME[[1]]$DD5[,1]))
  out6 <- rnorm(length(results_list$ME),y_seq6$condMean,sqrt(y_seq6$condVar)/
                  length(results_list$ME[[1]]$DD6[,1])) 
  
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob1[i,j] <- res_prob$samples[i,w1ind1+j-1]*(res_prob$samples[i,a1ind1+j-1]+
                                                   res_prob$samples[i,b11ind1+j-1]*PS_50[1]+
                                                   res_prob$samples[i,b12ind1+j-1]*PS_50[2]+
                                                   res_prob$samples[i,b13ind1+j-1]*PS_50[5]+
                                                   res_prob$samples[i,b14ind1+j-1]*PS_50[6])
    }
  }
  
  g1_50 <- sapply(rowMeans(prob1),expit)
  
  prob2 <- matrix(NA,length(results_list$ME),5)
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob1[i,j] <- res_prob$samples1[i,w2ind1+j-1]*(res_prob$samples1[i,a2ind1+j-1]+
                                                     res_prob$samples1[i,b21ind1+j-1]*PS_50[1]+
                                                     res_prob$samples1[i,b22ind1+j-1]*PS_50[2]+
                                                     res_prob$samples1[i,b23ind1+j-1]*PS_50[5]+
                                                     res_prob$samples1[i,b24ind1+j-1]*PS_50[6])
    }
  }
  
  g2_50 <- sapply(rowMeans(prob2),expit)
  
  theta_1_50 <- (out1*g1_50+out2*(1-g1_50))
  theta_2_50 <- (out1*g1_50+out3*(1-g1_50))
  theta_3_50 <- (out4*g2_50+out5*(1-g2_50))
  theta_4_50 <- (out4*g2_50+out6*(1-g2_50))
  
  #Compute means and standard errors for each embedded DTR
  theta_mean_50 <- cbind(theta_1_50,theta_2_50,theta_3_50,theta_4_50)
  colnames(theta_mean_50) <- c("EDTR 1","EDTR 2", 'EDTR 3',"EDTR 4")
  EDTR_50_mean <- apply(theta_mean_50,2,mean)
  EDTR_50_sd <- apply(theta_mean_50,2,sd)
  diff2<- -theta_mean_50-max(-EDTR_50_mean)
  quan2 <- apply(diff2,2,function(x) quantile(x,1-0.05/3))
  bst_50 <- which(quan2<0)
  
  
  y_seq1 <- condMVN(Mu_seq1,Sig_seq1,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_25[-c(3,4)])
  y_seq2 <- condMVN(Mu_seq2,Sig_seq2,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_25[-4])
  y_seq3 <- condMVN(Mu_seq3,Sig_seq3,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_25[-4])
  y_seq4 <- condMVN(Mu_seq4,Sig_seq4,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_25[-c(3,4)])
  y_seq5 <- condMVN(Mu_seq5,Sig_seq5,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_25[-3])
  y_seq6 <- condMVN(Mu_seq6,Sig_seq6,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_25[-3])
  
  out1 <- rnorm(length(results_list$ME),y_seq1$condMean,sqrt(y_seq1$condVar)/
                  length(results_list$ME[[1]]$DD1[,1]))
  out2 <- rnorm(length(results_list$ME),y_seq2$condMean,sqrt(y_seq2$condVar)/
                  length(results_list$ME[[1]]$DD2[,1]))
  out3 <- rnorm(length(results_list$ME),y_seq3$condMean,sqrt(y_seq3$condVar)/
                  length(results_list$ME[[1]]$DD3[,1]))
  out4 <- rnorm(length(results_list$ME),y_seq4$condMean,sqrt(y_seq4$condVar)/
                  length(results_list$ME[[1]]$DD4[,1]))
  out5 <- rnorm(length(results_list$ME),y_seq5$condMean,sqrt(y_seq5$condVar)/
                  length(results_list$ME[[1]]$DD5[,1]))
  out6 <- rnorm(length(results_list$ME),y_seq6$condMean,sqrt(y_seq6$condVar)/
                  length(results_list$ME[[1]]$DD6[,1])) 
  
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob1[i,j] <- res_prob$samples[i,w1ind1+j-1]*(res_prob$samples[i,a1ind1+j-1]+
                                                   res_prob$samples[i,b11ind1+j-1]*PS_25[1]+
                                                   res_prob$samples[i,b12ind1+j-1]*PS_25[2]+
                                                   res_prob$samples[i,b13ind1+j-1]*PS_25[5]+
                                                   res_prob$samples[i,b14ind1+j-1]*PS_25[6])
    }
  }
  
  g1_25 <- sapply(rowMeans(prob1),expit)
  
  prob2 <- matrix(NA,length(results_list$ME),5)
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob2[i,j] <- res_prob$samples1[i,w2ind1+j-1]*(res_prob$samples1[i,a2ind1+j-1]+
                                                     res_prob$samples1[i,b21ind1+j-1]*PS_25[1]+
                                                     res_prob$samples1[i,b22ind1+j-1]*PS_25[2]+
                                                     res_prob$samples1[i,b23ind1+j-1]*PS_25[5]+
                                                     res_prob$samples1[i,b24ind1+j-1]*PS_25[6])
    }
  }
  
  g2_25 <- sapply(rowMeans(prob2),expit)
  
  theta_1_25 <- out1*g1_25+out2*(1-g1_25)
  theta_2_25 <- out1*g1_25+out3*(1-g1_25)
  theta_3_25 <- out4*g2_25+out5*(1-g2_25)
  theta_4_25 <- out4*g2_25+out6*(1-g2_25)
  
  #Compute means and standard errors for each embedded DTR
  theta_mean_25 <- cbind(theta_1_25,theta_2_25,theta_3_25,theta_4_25)
  colnames(theta_mean_25) <- c("EDTR 1","EDTR 2", 'EDTR 3',"EDTR 4")
  EDTR_25_mean <- apply(theta_mean_25,2,mean)
  EDTR_25_sd <- apply(theta_mean_25,2,sd)
  diff3 <- -theta_mean_25-max(-EDTR_25_mean)
  quan3 <- apply(diff3,2,function(x) quantile(x,1-0.05/3))
  bst_25 <- which(quan3<0)
  
  
  y_seq1 <- condMVN(Mu_seq1,Sig_seq1,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_00[-c(3,4)])
  y_seq2 <- condMVN(Mu_seq2,Sig_seq2,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_00[-4])
  y_seq3 <- condMVN(Mu_seq3,Sig_seq3,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_00[-4])
  y_seq4 <- condMVN(Mu_seq4,Sig_seq4,dependent.ind=c(1),
                    given.ind=c(2:5),X.given = PS_00[-c(3,4)])
  y_seq5 <- condMVN(Mu_seq5,Sig_seq5,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_00[-3])
  y_seq6 <- condMVN(Mu_seq6,Sig_seq6,dependent.ind=c(1),given.ind=c(2:6),
                    X.given = PS_00[-3])
  
  out1 <- rnorm(length(results_list$ME),y_seq1$condMean,sqrt(y_seq1$condVar)/
                  length(results_list$ME[[1]]$DD1[,1]))
  out2 <- rnorm(length(results_list$ME),y_seq2$condMean,sqrt(y_seq2$condVar)/
                  length(results_list$ME[[1]]$DD2[,1]))
  out3 <- rnorm(length(results_list$ME),y_seq3$condMean,sqrt(y_seq3$condVar)/
                  length(results_list$ME[[1]]$DD3[,1]))
  out4 <- rnorm(length(results_list$ME),y_seq4$condMean,sqrt(y_seq4$condVar)/
                  length(results_list$ME[[1]]$DD4[,1]))
  out5 <- rnorm(length(results_list$ME),y_seq5$condMean,sqrt(y_seq5$condVar)/
                  length(results_list$ME[[1]]$DD5[,1]))
  out6 <- rnorm(length(results_list$ME),y_seq6$condMean,sqrt(y_seq6$condVar)/
                  length(results_list$ME[[1]]$DD6[,1])) 
  
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob1[i,j] <- res_prob$samples[i,w1ind1+j-1]*(res_prob$samples[i,a1ind1+j-1]+
                                                   res_prob$samples[i,b11ind1+j-1]*PS_00[1]+
                                                   res_prob$samples[i,b12ind1+j-1]*PS_00[2]+
                                                   res_prob$samples[i,b13ind1+j-1]*PS_00[5]+
                                                   res_prob$samples[i,b14ind1+j-1]*PS_00[6])
    }
  }
  
  g1_00 <- sapply(rowMeans(prob1),expit)
  
  for(i in 1:length(results_list$ME)){
    for(j in 1:5){
      prob2[i,j] <- res_prob$samples1[i,w2ind1+j-1]*(res_prob$samples1[i,a2ind1+j-1]+
                                                     res_prob$samples1[i,b21ind1+j-1]*PS_00[1]+
                                                     res_prob$samples1[i,b22ind1+j-1]*PS_00[2]+
                                                     res_prob$samples1[i,b23ind1+j-1]*PS_00[5]+
                                                     res_prob$samples1[i,b24ind1+j-1]*PS_00[6])
    }
  }
  
  g2_00 <- sapply(rowMeans(prob2),expit)
  
  
  theta_1_00 <- out1*g1_00+out2*(1-g1_00)
  theta_2_00 <- out1*g1_00+out3*(1-g1_00)
  theta_3_00 <- out4*g2_00+out5*(1-g2_00)
  theta_4_00 <- out4*g2_00+out6*(1-g2_00)
  
  #Compute means and standard errors for each embedded DTR
  theta_mean_00 <- cbind(theta_1_00,theta_2_00,theta_3_00,theta_4_00)
  colnames(theta_mean_00) <- c("EDTR 1","EDTR 2", 'EDTR 3',"EDTR 4")
  EDTR_00_mean <- apply(theta_mean_00,2,mean)
  EDTR_00_sd <- apply(theta_mean_00,2,sd)
  diff4 <- -theta_mean_00-max(-EDTR_00_mean)
  quan4 <- apply(diff4,2,function(x) quantile(x,1-0.05/3))
  bst_00 <- which(quan4<0)
  return(list(EDTR_00_mean,EDTR_75_mean,EDTR_50_mean,EDTR_25_mean,bst_00,bst_75,
              bst_50,bst_25))
}

#Running in parallel
library(doParallel)  
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores, type="FORK")  
result <- parLapply(cl, sim_1, EDTR_fun)  
stopCluster(cl)


#True EDTRs
true_dtr_00_1 <- 0.5*(0.7+0.6*exp(1+1)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                  0.2+0.7*1+0.7*1-0.9*median(data[,"X1"])+0.3*median(data[,"X2"])+0.7*median(data[,"X3"])+0.9*1)
true_dtr_00_2 <- 0.5*(0.7+0.6*exp(1+1)+0.8*median(dt[,"X1"])-0.2*median(dt[,"X2"])+
                        0.2+0.2+0.6*1+0.7*1+0.9*median(data[,"X1"])+0.2*median(data[,"X2"])+0.6*median(data[,"X3"])+0.8*1)
true_dtr_00_3 <- 0.5*(0.7+0.6*1+0.6*1+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.5*1+0.6*1+0.7*log(1+1)-0.5*median(data[,"X2"])+median(data[,"X3"]))
true_dtr_00_4 <- 0.5*(0.7+0.6*1+0.6*1+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.8*1+0.7*1+0.3*1-0.5*median(data[,"X2"])+0.9*median(data[,"X3"]))

true_dtr_75_1 <- 0.5*(0.7+0.6*exp(1+0.75)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.2+0.7*0.75+0.7*0.75-0.9*median(data[,"X1"])+0.3*median(data[,"X2"])+0.7*median(data[,"X3"])+0.9*0.75)
true_dtr_75_2 <- 0.5*(0.7+0.6*exp(1+0.75)+0.8*median(dt[,"X1"])-0.2*median(dt[,"X2"])+
                        0.2+0.2+0.6*0.75+0.7*0.75+0.9*median(data[,"X1"])+0.2*median(data[,"X2"])+0.6*median(data[,"X3"])+0.8*0.75)
true_dtr_75_3 <- 0.5*(0.7+0.6*0.75+0.6*0.75+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.5*0.75+0.6*0.75+0.7*log(1+0.75)-0.5*median(data[,"X2"])+median(data[,"X3"]))
true_dtr_75_4 <- 0.5*(0.7+0.6*0.75+0.6*0.75+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.8*0.75+0.7*0.75+0.3*0.75-0.5*median(data[,"X2"])+0.9*median(data[,"X3"]))

true_dtr_75_1 <- 0.5*(0.7+0.6*exp(1+0.75)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.2+0.7*0.75+0.7*0.75-0.9*median(data[,"X1"])+0.3*median(data[,"X2"])+0.7*median(data[,"X3"])+0.9*0.75)
true_dtr_75_2 <- 0.5*(0.7+0.6*exp(1+0.75)+0.8*median(dt[,"X1"])-0.2*median(dt[,"X2"])+
                        0.2+0.2+0.6*0.75+0.7*0.75+0.9*median(data[,"X1"])+0.2*median(data[,"X2"])+0.6*median(data[,"X3"])+0.8*0.75)
true_dtr_75_3 <- 0.5*(0.7+0.6*0.75+0.6*0.75+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.5*0.75+0.6*0.75+0.7*log(1+0.75)-0.5*median(data[,"X2"])+median(data[,"X3"]))
true_dtr_75_4 <- 0.5*(0.7+0.6*0.75+0.6*0.75+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.8*0.75+0.7*0.75+0.3*0.75-0.5*median(data[,"X2"])+0.9*median(data[,"X3"]))


true_dtr_50_1 <- 0.5*(0.7+0.6*exp(1+0.50)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.2+0.7*0.50+0.7*0.50-0.9*median(data[,"X1"])+0.3*median(data[,"X2"])+0.7*median(data[,"X3"])+0.9*0.50)
true_dtr_50_2 <- 0.5*(0.7+0.6*exp(1+0.50)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.2+0.2+0.6*0.5+0.7*0.5+0.9*median(data[,"X1"])+0.2*median(data[,"X2"])+0.6*median(data[,"X3"])+0.8*0.5)
true_dtr_50_3 <- 0.5*(0.7+0.6*0.5+0.6*0.5+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.5*0.5+0.6*0.5+0.7*log(1+0.5)-0.5*median(data[,"X2"])+median(data[,"X3"]))
true_dtr_50_4 <- 0.5*(0.7+0.6*0.5+0.6*0.5+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.8*0.5+0.7*0.5+0.3*0.5-0.5*median(data[,"X2"])+0.9*median(data[,"X3"]))

true_dtr_25_1 <- 0.5*(0.7+0.6*exp(1+0.25)+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.2+0.7*0.25+0.7*0.25-0.9*median(data[,"X1"])+0.3*median(data[,"X2"])+0.7*median(data[,"X3"])+0.9*0.25)
true_dtr_25_2 <- 0.5*(0.7+0.6*exp(1+0.25)+0.8*median(data[,"X1"])-0.2*median(dt[,"X2"])+
                        0.2+0.2+0.6*0.25+0.7*0.25+0.9*median(data[,"X1"])+0.2*median(data[,"X2"])+0.6*median(data[,"X3"])+0.8*0.25)
true_dtr_25_3 <- 0.5*(0.7+0.6*0.25+0.6*0.25+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.5*0.25+0.6*0.25+0.7*log(1+0.25)-0.5*median(data[,"X2"])+median(data[,"X3"]))
true_dtr_25_4 <- 0.5*(0.7+0.6*0.25+0.6*0.25+0.8*median(data[,"X1"])-0.2*median(data[,"X2"])+
                        0.3+0.8*0.25+0.7*0.25+0.3*0.25-0.5*median(data[,"X2"])+0.9*median(data[,"X3"]))



