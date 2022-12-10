###############################################################################################################
## GROWTH MODELS FOR PIT-TAGGED JUVENILE SRSS CHINOOK 
###############################################################################################################

data
{
  for(y in 1:N.years)
  {
    ones[y] <- 1
  }
}

model
{
  T1 <- 0
  T2 <- 300
  
  #------------------------------------------------------------------------------
  # HYPER-PRIORS
  #------------------------------------------------------------------------------
  
  # Hyper-means and hyper-SDs of growth curve parameters
  mu.lL1 ~ dnorm(0,1e-6)           # length at reference age T1 that "anchors" the growth curve
  sigma.lL1 ~ dunif(0,10)
  mu.lL2 ~ dnorm(0,1e-6)T(mu.lL1,) # length at reference age T2 that "anchors" the growth curve
  sigma.lL2 ~ dunif(0,10)
  
  mu.lb ~ dnorm(0,1e-6)T(0,)       # shape parameter
  sigma.lb ~ dunif(0,10)
  
  sigma ~ dunif(0,10)              # residual SD
  
  #------------------------------------------------------------------------------
  # RANDOM EFFECTS FOR EACH YEAR
  #------------------------------------------------------------------------------
  
  for(y in 1:N.years)
  {
    L1[y] ~ dlnorm(mu.lL1, pow(sigma.lL1,-2))
    L20[y] ~ dlnorm(mu.lL2, pow(sigma.lL2,-2))T(L1[y],)
    L2[y] <- max(L1[y], L20[y])    # this SHOULD be redundant, as test[1:N.years] confirms...but it isn't
    test[y] <- L2[y] == L20[y]
    #     ones[y] ~ dbern(test[y])
    b[y] ~ dlnorm(mu.lb, pow(sigma.lb,-2))T(0,)
  }
  
  #------------------------------------------------------------------------------
  # LIKELIHOOD OF OBSERVATIONS
  #------------------------------------------------------------------------------
  
  # Mark-recapture observations of initial and final length and elapsed time
  for(i in 1:length(LR))
  {
    # Growth curve parameter values for this smolt year
    # (avoids clunky expression for fitted value)
    L1i[i] <- L1[smolt.year.MR[i]]
    L2i[i] <- L2[smolt.year.MR[i]]
    bi[i] <- b[smolt.year.MR[i]]
    
    # Predicted length at recapture
    LR.fit[i] <- (LM[i]^bi[i] + (L2i[i]^bi[i] - L1i[i]^bi[i])*dt[i]/(T2 - T1))^(1/bi[i])
    
    # Likelihood
    LR[i] ~ dlnorm(log(LR.fit[i]), pow(sigma,-2))
  }
  
  # Measurements of individuals that were only observed once
  for(j in 1:length(L))
  {
    # Growth curve parameter values for M-R fish in this smolt year
    # (avoids clunky expression for fitted value)
    L1j[j] <- L1[smolt.year.noMR[j]]
    L2j[j] <- L2[smolt.year.noMR[j]]
    bj[j] <- b[smolt.year.noMR[j]]
    
    # Predicted length when observed
    L.fit[j] <- (L1j[j]^bj[j] + (L2j[j]^bj[j] - L1j[j]^bj[j])*(tt[j] - T1)/(T2 - T1))^(1/bj[j])
    
    # Likelihood
    L[j] ~ dlnorm(log(L.fit[j]), pow(sigma,-2))
  } 
}

