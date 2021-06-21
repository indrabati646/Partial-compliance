source("/Users/indrabatibhattacharya/Documents/biom_code/bnp_model_Examples1&2.R")
source("/Users/indrabatibhattacharya/Documents/biom_code/Simulate_Examples1&2.R")

results_list <- lapply(sim_1,bnp_mod)
bias_1 <- bias_2 <- bias_3 <- bias_4 <- bias_5 <- bias_6 <- c()
sd_1 <- sd_2 <- sd_3 <- sd_4 <- sd_5 <- sd_6 <- c()

for(i in 1:length(sim_1)){
  data<- sim_1[[i]]
  data <- matrix(unlist(data),nrow=nrow(data),ncol=ncol(data))
  colnames(data) <- c("a1","a2","s","D1","D2","D3","y","y1","X1","X2","X3")
  
  y1 <- data[which(data[,"a1"]==1 & data[,"s"]==1),]
  y2 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==1),]
  y3 <- data[which(data[,"a1"]==1 & data[,"s"]==0 & data[,"a2"]==-1),]
  y4 <- data[which(data[,"a1"]==-1 & data[,"s"]==1),]
  y5 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==1),]
  y6 <- data[which(data[,"a1"]==-1 & data[,"s"]==0 & data[,"a2"]==-1),]
  x1 <- data[,"X1"]; x2 <- data[,"X2"]; x3 <- data[,"X3"]
  x <- cbind(x1,x2,x3)
  n1 <- nrow(y1); n2 <- nrow(y2); n3 <- nrow(y3); n4 <- nrow(y4) 
  n5 <- nrow(y5); n6 <- nrow(y6)
  n <- n1+n2+n3+n4+n5+n6 
  y <- data[,"y1"]
  yy1 <- y1[,"y1"]; yy2 <- y2[,"y1"]; yy3 <- y3[,"y1"]; yy4 <- y4[,"y1"]
  yy5 <- y5[,"y1"]; yy6 <- y6[,"y1"]
  ymat <- matrix(unlist(results_list[[i]][[1]]),ncol=length(yy1),byrow=T)
  ymean1 <- apply(ymat,2,mean)
  yvar1 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_1[i]<- mean(yy1-ymean1)
  sd_1[i]<- mean(sqrt(yvar1))
  
  ymat <- matrix(unlist(results_list[[i]][[2]]),ncol=length(yy2),byrow=T)
  ymean2 <- apply(ymat,2,mean)
  yvar2 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_2[i]<- mean(yy2-ymean2)
  sd_2[i]<- mean(sqrt(yvar2))
  
  ymat <- matrix(unlist(results_list[[i]][[3]]),ncol=length(yy3),byrow=T)
  ymean3 <- apply(ymat,2,mean)
  yvar3 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_3[i] <- mean(yy3-ymean3)
  sd_3[i] <- mean(sqrt(yvar3))
  
  ymat <- matrix(unlist(results_list[[i]][[4]]),ncol=length(yy4),byrow=T)
  ymean4 <- apply(ymat,2,mean)
  yvar4 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_4[i] <- mean(yy4-ymean4)
  sd_4[i] <-mean(sqrt(yvar4))
  
  
  ymat <- matrix(unlist(results_list[[i]][[5]]),ncol=length(yy5),byrow=T)
  ymean5 <- apply(ymat,2,mean)
  yvar5 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_5[i] <- mean(yy5-ymean5)
  sd_5[i] <- mean(sqrt(yvar5))
  
  ymat <- matrix(unlist(results_list[[i]][[6]]),ncol=length(yy6),byrow=T)
  ymean6 <- apply(ymat,2,mean)
  yvar6 <- apply(ymat,2,function(x) var(x)*(length(x)-1)/length(x))
  bias_6[i] <- mean(yy6-ymean6)
  sd_6[i] <- mean(sqrt(yvar6))
}

bias_1_final <- mean(bias_1); sd_1_final <- mean(sd_1)
bias_2_final <- mean(bias_2); sd_2_final <- mean(sd_2)
bias_3_final <- mean(bias_3); sd_3_final <- mean(sd_3)
bias_4_final <- mean(bias_4); sd_4_final <- mean(sd_4)
bias_5_final <- mean(bias_5); sd_5_final <- mean(sd_5)
bias_6_final <- mean(bias_6); sd_6_final <- mean(sd_6)

