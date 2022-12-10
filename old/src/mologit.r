model mologit {
for (i in 1:n){    
    Y[i] ~ dcat(P[i,])
    P[i,1] <- max(min(1 - Q[i,1],1),0)
    
    for (i.cut in 2:n.cut){
      P[i,i.cut] <- Q[i,i.cut-1] - Q[i,i.cut]
    }
    P[i,n.cut+1] <- max(min(Q[i,n.cut],1),0)
    
  # random effects!!
    for (i.cut in 1:n.cut){
      logit(Q[i,i.cut]) <- Z[i,i.cut]
      Z[i,i.cut] <- b1*age[i] + b2*edu[i] + b3*inc[i] - C[yr[i],i.cut] 
    }
  }
  
  # PRIORS: 
  b1 ~ dnorm(0,.0001)
  b2 ~ dnorm(0,.0001)
  b3 ~ dnorm(0,.0001)
  
  for (i.year in 1:n.year){
    C[i.year,1] ~ dnorm(0,.0001)I(,C[i.year,2])
    C[i.year,2] ~ dnorm(0,.0001)I(C[i.year,1],C[i.year,3])
    C[i.year,3] ~ dnorm(0,.0001)I(C[i.year,2],)
  }

}
