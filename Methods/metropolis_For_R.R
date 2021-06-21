#Metropolis algorithm for the copula correlation matrix R


metropolisC1 = function(h,rho){
  acc1 <- rep(NA,6)
  prop1 <- rho 
  
  mem <- function(prop1){
    mempty<-matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],
                     1,prop1[6],prop1[3],prop1[5],prop1[6],1),4,4,byrow=TRUE)
    return(mempty)
  }
  for(q in 1:6){
    prop1[q] <- 1
    RHO1 <- det(mem(prop1))
    prop1[q] <- 0
    RHO0 <- det(mem(prop1))
    prop1[q] <- -1
    RHO_1 <- det(mem(prop1))
    Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
    interval <- multiroot(Fun, start=c(-1.1,1.1), maxiter = 200)$root  # Interval satisfying positive definite matrix restriction
    
     prop1[q] <- runif(1, rho[q]-0.3,rho[q]+0.3 )      
    
    # Proposed R
    COR <-  mem(prop1)
    
    # Current R
    RHO <-  mem(rho)
    
    # Acceptance ratio
    if(is.positive.definite(COR)){        
      rat <- -(nrow(h))/2*log(det(COR))-
        0.5*sum(diag(h%*%(solve(COR))%*%t(h)))+
        dunif(prop1[q],min(interval),max(interval),log=TRUE)+
        dunif(rho[q],prop1[q]-0.3,prop1[q]+0.3,log=TRUE)+
        (nrow(h))/2*log(det(RHO))+0.5*sum(diag(h%*%(solve(RHO))%*%t(h)))-
        dunif(rho[q],min(interval),max(interval),log=TRUE)-
        dunif(prop1[q],rho[q]-0.3,rho[q]+0.3,log=TRUE)
      if(is.na(rat)){
        prop1[q] <- rho[q]
        acc1[q] <- 0
      }else{
        if(log(runif(1))>rat) {
          prop1[q] <- rho[q]
          acc1[q] <- 0
        }else{
          rho[q] <- prop1[q]}
      }
    }else{
      prop1[q] <- rho[q]
      acc1[q] <- 0
    }    
  }
  return(c(acc1,prop1))
}

