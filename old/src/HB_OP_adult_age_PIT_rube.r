###############################################################################################################
## HIERARCHICAL MODEL OF SALMON ADULT AGE FOR PIT-TAG DATA
## 
## This model uses observations of adult age-at-return of fish that were tagged as outmigrating juveniles
## to make inferences about the time-series properties of the unobserved "states" (conditional probabilities
## of returning at a given age, conditioned on surviving to adulthood). The state vector in year y consists
## of age-at-return probabilities p[1:M], where M = number of ocean ages. These are represented as a partition 
## of an underlying continuous latent variable alpha. Following the standard "latent variable" formulation of the
## ordered probit model, alpha is normally distributed conditional on the predictors, and the partition is achieved
## by a vector of cutpoints tau[1:(M-1)].
##
## The process model for alpha includes two "levels" corresponding to different time scales, (1) interannual
## (i.e., across outmigrant cohorts) and (2) intra-annual (within-cohort variation among individual fish):
##
##   alpha[i] = alpha.yr[y] + alpha.ind[i].
##
## The cohort-level submodel is a linear regression with AR(1) errors:
##
##   alpha.yr[y] = alpha.mu[y] + rho * (alpha.yr[y-1] - alpha.mu[y-1]) + v[y]
##
## where rho is an AR(1) coefficient (possibly zero) and v[y] ~ N(0, sigma.yr).
##
## The annual mean includes any covariate effects given by the matrix X (these
## might include, e.g., PDO or SST) with a vector of regression coefficients a:
##
##   alpha.mu[y] = a0.int + X[y,1:numX] %*% a.
##
## Finally, the individual-level submodel is given by another linear predictor, with covariate matrix Z
## (e.g., julian day or size at migration) and a vector of regression coefficients b:
##
##   alpha.ind[i] = Z[i,1:numZ] %*% b + e[i]
##
## where e[i] ~ N(0, sigma).
##
## (NOTE: We may eventually allow the b-elements to be time-varying, but for now they are constant.)
##
## As implemented here, the model is overparameterized; the intercept of alpha, the full set of M-1 cutpoints,
## and the "probit" SD sigma are not jointly identifiable. This redundant parameterization greatly improves mixing,
## but we transform to the set of identified parameters (i.e., the ones that should be monitored) as follows,
## where the "0" suffix denotes the redundant, top-level parameterization:
## 
##   a = (a0)/sigma
##   b = (b0)/sigma
##   tau = (tau0 - a0.int)/sigma
##
## The data consist of observations of individual adult age, which we model as multinomial given the prob vector p:
##
##   age[i] ~ Multinom(1,p).
##
##
## DATA:
##   X          year x variable matrix of cohort-level covariates. If no covars, use a column matrix of zeros.
##   Z          individual x variable matrix of indiv-level covariates. If no covars, use a column matrix of zeros.
##   IID        0/1 scalar indicating whether to include autocorrelation in process model (IID==0) or not (IID==1)
##   year       vector identifying the migration year of each individual, shifted so the first year is 1
##   age        vector of individual observations of adult ocean age, shifted if necessary so youngest age is 1
###############################################################################################################


data
{  
  Y <- max(year)       # number of outmigration years in dataset
  M <- max(age)        # dimension of multinomial probability vector
    
#   dimX <- dim(X)       # year x variable matrix of covariates
#   numX <- dimX[2]      # number of cohort-level (interannual) covariates 
# 
  dimZ <- dim(Z)       # individual x variable matrix of covariates
  numZ <- dimZ[2]      # number of individual-level (intra-annual) covariates
}


model
{
  #------------------------------------------------------------------------------
  # PRIORS
  #------------------------------------------------------------------------------
  
  a0.int ~ dnorm(0,1e-6)    # intercept (redundant parameterization)
  
  # Cohort-level regression coefficients
#   for(k in 1:numX)
#   {
#     a0[k] ~ dnorm(0,1e-6)
#   }
  
  FOR(a0, LC.alpha0.mu, fuck, ? ~ dnorm(0,1e-6))    # cohort-level regression coefs (rube FOR syntax)
  a <- a0/sigma  # map to identified parameters
  
  # Individual-level regression coefficients
  # (CURRENTLY TIME-INVARIANT; MAKE RANDOM?)
  for(k in 1:numZ)
  {
    b0[k] ~ dnorm(0,1e-6)
  }
  
#   FOR(b0, LC.alpha0.ind, , ? ~ dnorm(0,1e-6))  # individual-level regression coefs (rube FOR syntax)
  b <- b0/sigma  # map to identified parameters
  
  rho0 ~ dunif(-1,1)           # AR(1) coefs
  rho <- rho0*(1 - IID)        # no autocorrelation if IID==1
  sigma.yr0 ~ dunif(0,10)      # process error SD
  sigma.yr <- sigma.yr0/sigma  # map to identified parameters
  sigma ~ dunif(0,10)          # "residual" or individual-level error SD (redundant parameterization)
  
  # Cutpoints
  for(j in 1:(M-1))
  {
    tau.vals[j] ~ dnorm(0,1e-6)
  }
  tau0 <- sort(tau.vals)
  tau <- (tau0 - a0.int)/sigma  # map to identified parameters
  
  #---------------------------------------------------------------------------------------
  # COHORT-LEVEL PROCESS MODEL
  # Loop over outmigration years, evaluate likelihood of annual state (latent age indicator) 
  # under the process model. Note that initial state distribution must be specified 
  # separately to initialize the process model.
  #---------------------------------------------------------------------------------------
  
  alpha0.mu[1] <- a0.int + sum(X[1,1:length(a0)]*a0)     # expected alpha0 for year 1
  alpha0.mu[1] <- LC(a0, LC.alpha0.mu, , 1)              # expected alpha (rube LC syntax)
  alpha0.yr[1] ~ dnorm(alpha0.mu[1], pow(sigma.yr0,-2))  # alpha0 for year 1
  alpha0.err[1] <- alpha0.yr[1] - alpha0.mu[1]           # alpha0 "residual" for year 1
  
  # partition alpha0 using cutpoints
  F.yr[1,1] <- 0
  for(j in 1:(M-1)) 
  {
    F.yr[1,j+1] <- phi((tau0[j] - alpha0.yr[1])/sigma)
    p.yr[1,j] <- max(min(F.yr[1,j+1] - F.yr[1,j], 1), 1e-10)
  }
  p.yr[1,M] <- max(min(1 - F.yr[1,M], 1), 1e-10)

  for(y in 2:Y)
  {
    alpha0.mu[y] <- a0.int + sum(X[y,1:length(a0)]*a0) # expected alpha0 for year y
#     alpha0.mu[y] <- LC(a0, LC.alpha0.mu, , y)        # expected alpha (rube LC syntax)
    alpha0.yr[y] ~ dnorm(alpha0.mu[y] + rho*alpha0.err[y-1], pow(sigma.yr0,-2))   # alpha for year y
    alpha0.err[y] <- alpha0.yr[y] - alpha0.mu[y]     # alpha0 "residual" for year y
    
    # partition alpha0 using cutpoints
    F.yr[y,1] <- 0
    for(j in 1:(M-1))                                   
    {
      F.yr[y,j+1] <- phi((tau0[j] - alpha0.yr[y])/sigma)
      p.yr[y,j] <- max(min(F.yr[y,j+1] - F.yr[y,j], 1), 1e-10)
    }
    p.yr[y,M] <- max(min(1 - F.yr[y,M], 1), 1e-10)
  }
  
  alpha.yr <- (alpha0.yr - a0.int)/sigma             # rescaled (identified) alpha
  
  #---------------------------------------------------------------------------------------
  # INDIVIDUAL-LEVEL PROCESS MODEL AND LIKELIHOOD
  # Loop over all individual fish. Adjust "expected" state vector for the cohort (i.e.,
  # smolt year) by individual-level covariate effects. Evaluate multinomial likelihood
  # of observed ocean age at adult return.
  #---------------------------------------------------------------------------------------
  
  for(i in 1:length(age))
  {
#     alpha0.ind[i] <- LC(b0, LC.alpha0.ind, , i)            # individual effects (rube LC syntax)
#     alpha0[i] <- alpha0.yr[year[i]] + alpha0.ind[i]        # alpha0 for individual i
    alpha0[i] <- alpha0.yr[year[i]] + sum(Z[i,1:length(b0)]*b0)  # alpha0 for individual i
    
    # partition alpha0 using cutpoints
    F[i,1] <- 0
    for(j in 1:(M-1))                                   
    {
      F[i,j+1] <- phi((tau0[j] - alpha0[i])/sigma)
      p[i,j] <- max(min(F[i,j+1] - F[i,j], 1), 1e-10)
    }
    p[i,M] <- max(min(1 - F[i,M], 1), 1e-10)

    age[i] ~ dcat(p[i,1:M])                                # likelihood for individual i
    cc[i] <- p[i,age[i]]==max(p[i,1:M])                    # correct classification score (0/1)
  }
  ccr <- mean(cc)                                          # overall correct classification rate
  
  # Got all that?
}

