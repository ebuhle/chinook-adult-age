###############################################################################################################
## HIERARCHICAL MODEL OF SALMON ADULT AGE FOR PIT-TAG DATA
## 
## This model uses observations of adult age-at-return of fish that were tagged as outmigrating juveniles
## to make inferences about the time-series properties of the unobserved "states" (conditional probabilities
## of returning at a given age, conditioned on surviving to adulthood). The state vector in year y consists
## of log ratio-transformed age-at-return probabilities lp[1:M], where M = number of ocean ages - 1.
## The multivariate process model includes two "levels" corresponding to different time scales, (1) interannual
## (i.e., across outmigrant cohorts) and (2) intra-annual (within-cohort variation among individual fish):
##
##   lp[i,1:M] = lp.yr[y,1:M] + lp.ind[i,1:M].
##
## The cohort-level submodel is a seemingly unrelated regression (SUR) with AR(1) errors. Dropping the state
## indexing for simplicity:
##
##   lp.yr[y,] = lp.mu[y,] + Phi %*% (lp.yr[y-1,] - lp.mu[y-1,]) + v[y,]
##
## where Phi is a diagonal matrix of AR(1) coefficients (possibly zero) and 
## v[y,] ~ MVN(0, Sigma), with Sigma a (possibly diagonal) MxM process error covariance matrix.
##
## The annual mean vector includes an intercept lp0 and any covariate effects given by the matrix X (these
## might include, e.g., PDO or SST) with a numX x M matrix of regression coefficients a, 
## some of which may be set to zero so that different state components depend on different covariates:
##
##   lp.mu[y,] = lp0 + X[y,1:numX] %*% a.
##
## Finally, the individual-level submodel is given by another SUR, with covariate matrix Z (these might
## include, e.g., julian day or size at migration) and a numZ x M matrix of regression coefficients b,
## some of which may be set to zero so that different state components depend on different covariates:
##
##   lp.ind[i,] = Z[i,1:numZ] %*% b.
##
## (NOTE: We may eventually allow the b-elements to be time-varying, but for now they are constant.)
##
## The data consist of observations of individual adult age, which we model as multinomial given the probability
## vector p (i.e., the inverse log ratio of lp):
##
##   age[i] ~ Multinom(1,p).
##
##
## DATA:
##   X          year x variable matrix of cohort-level covariates. If no covars, use a column matrix of zeros.
##   Xindx      variable x state 0/1 matrix. The regression model for state j includes X[,k] iff Xindx[k,j]==1.
##   Z          individual x variable matrix of indiv-level covariates. If no covars, use a column matrix of zeros.
##   Zindx      variable x state 0/1 matrix. The regression model for state j includes Z[,k] iff Zindx[k,j]==1.
##   IID        0/1 scalar indicating whether to include autocorrelation in process model (IID==0) or not (IID==1)
##   year       vector identifying the migration year of each individual, shifted so the first year is 1
##   age        vector of individual observations of adult ocean age, shifted if necessary so youngest age is 1
###############################################################################################################


data
{  
  Y <- max(year)       # number of outmigration years in dataset
  M <- max(age) - 1    # dimension of state vector
    
  dimX <- dim(X)       # year x variable matrix of covariates
  numX <- dimX[2]      # number of cohort-level (interannual) covariates 

  dimZ <- dim(Z)       # individual x variable matrix of covariates
  numZ <- dimZ[2]      # number of individual-level (intra-annual) covariates
}


model
{
  #------------------------------------------------------------------------------
  # PRIORS
  #------------------------------------------------------------------------------
  
  # Regression coefficients for cohort- and individual-level process model
  for(j in 1:M)
  {
    lp0[j] ~ dnorm(0,1e-6)          # intercept vector
    
    # cohort-level covariate slopes
    for(k in 1:numX)
    {
      a0[k,j] ~ dnorm(0,1e-6)       # numX x M matrix of coefs
      a[k,j] <- a0[k,j]*Xindx[k,j]  # set to zero if specified via Xindx
    }
    
    phi0[j] ~ dunif(-1,1)           # AR(1) coefs
    phi[j] <- phi0[j]*(1 - IID)     # no autocorrelation if IID==1
    
    # individual-level covariate slopes
    for(k in 1:numZ)
    {
      b0[k,j] ~ dnorm(0,1e-6)       # numZ x M matrix of coefs (CURRENTLY TIME-INVARIANT; MAKE RANDOM?)
      b[k,j] <- b0[k,j]*Zindx[k,j]  # set to zero if specified via Zindx
    }
  }
  
  ### CHECK SENSITIVITY TO PRIOR SCALE
  ### CAN invSigma BE MADE DIAGONAL? JUST MULTIPLYING ELEMENT-WISE BY IDENTITY MATRIX DOESN'T WORK
  ### B/C JAGS CANNOT FIND APPROPRIATE SAMPLER FOR WISHART NODE IN THAT CASE

  # construct prior covariance scale matrix
  for(i in 1:M)
  {
    for(j in 1:M)
    {
      scaleSigma[i,j] <- 1*(i==j)
    }
  }
  
  invSigma ~ dwish(scaleSigma,M+1)  # precision (inverse covariance) matrix for cohort-level process model

  #---------------------------------------------------------------------------------------
  # COHORT-LEVEL PROCESS MODEL
  # Loop over outmigration years, evaluate likelihood of annual states 
  # (alr-transformed age-at-return probabilities) under SUR-AR(1) process model. 
  # Note that initial state distribution must be specified separately to initialize 
  # the process model.
  #---------------------------------------------------------------------------------------
  
  lp.mu[1,1:M] <- lp0 + X[1,1:numX] %*% a         # expected alr(p) for year 1
  lp.yr[1,1:M] ~ dmnorm(lp.mu[1,1:M], invSigma)   # alr(p) for year 1
  lp.err[1,1:M] <- lp.yr[1,1:M] - lp.mu[1,1:M]    # alr(p) "residual" for year 1
  
  for(j in 1:M)                                   # inverse-alr transform for convenience
  {
    elp.yr[1,j] <- exp(lp.yr[1,j])
  }
  elp.yr[1,M+1] <- 1
  p.yr[1,1:(M+1)] <- elp.yr[1,1:(M+1)]/sum(elp.yr[1,1:(M+1)])      # age-at-return for year y    
  
  for(y in 2:Y)
  {
    lp.mu[y,1:M] <- lp0 + X[y,1:numX] %*% a       # expected alr(p) for year y
    lp.yr[y,1:M] ~ dmnorm(lp.mu[y,1:M] + phi*lp.err[y-1,1:M], invSigma)   # alr(p) for year y
    lp.err[y,1:M] <- lp.yr[y,1:M] - lp.mu[y,1:M]  # alr(p) "residual" for year y
    
    for(j in 1:M)                                 # inverse-alr transform for convenience
    {
      elp.yr[y,j] <- exp(lp.yr[y,j])
    }
    elp.yr[y,M+1] <- 1
    p.yr[y,1:(M+1)] <- elp.yr[y,1:(M+1)]/sum(elp.yr[y,1:(M+1)])     # age-at-return for year y    
  }
  
  
  #---------------------------------------------------------------------------------------
  # INDIVIDUAL-LEVEL PROCESS MODEL AND LIKELIHOOD
  # Loop over all individual fish. Adjust "expected" state vector for the cohort (i.e.,
  # smolt year) by individual-level covariate effects. Evaluate multinomial likelihood
  # of observed ocean age at adult return.
  #---------------------------------------------------------------------------------------
  
  for(i in 1:length(age))
  {
    lp[i,1:M] <- lp.yr[year[i],1:M] + Z[i,1:numZ] %*% b  # alr(p) for individual i
    
    for(j in 1:M)  # loop needed because "You cannot vectorize link functions. Sorry."
    {
      elp[i,j] <- exp(lp[i,j])
    }
    elp[i,M+1] <- 1
    p[i,1:(M+1)] <- elp[i,1:(M+1)]/sum(elp[i,1:(M+1)])   # age-at-return for individual i
    age[i] ~ dcat(p[i,1:(M+1)])                          # likelihood for individual i
    cc[i] <- p[i,age[i]]==max(p[i,1:(M+1)])              # correct classification score (0/1)
  }
  
  ccr <- mean(cc)                                        # overall correct classification rate
  
  # Got all that?
}





