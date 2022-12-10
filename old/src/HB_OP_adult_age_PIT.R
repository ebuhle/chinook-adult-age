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
##   alpha.ind[i] = Z[i,1:numZ] %*% b[y,1:numZ] + e[i]
##
## where e[i] ~ N(0, sigma).
##
## The b coefficients are time-varying, so each individual-level effect varies from year to year
## following an uncorrelated random walk (essentially a DLM for the latent state):
##
##   b[1,k] ~ N(0,1e-6)
##   b[y,k] ~ N(b[y-1,k], sigma.b[k]).
##
## As implemented here, the model is overparameterized; the intercept of alpha, the full set of M-1 cutpoints,
## and the "probit" SD sigma are not jointly identifiable. This redundant parameterization greatly improves mixing,
## but we transform to the set of identified parameters (i.e., the ones that should be monitored) as follows,
## where the "0" suffix denotes the redundant, top-level parameterization:
## 
##   a = a0/sigma
##   sigma.b = sigma.b0/sigma
##   b = b0/sigma
##   sigma.yr = sigma0.yr/sigma
##   tau = (tau0 - a0.int)/sigma
##
## The data consist of observations of individual adult age, which we model as multinomial given the prob vector p:
##
##   age[i] ~ Multinom(1,p).
##
## DATA:
##   X          year x variable data frame of cohort-level covariates. If no covars, use a column of zeros.
##   Z          individual x variable data frame of indiv-level covariates. If no covars, use a column of zeros.
##   year       vector identifying the migration year of each individual, shifted so the first year is 1
##   age        vector of individual observations of adult ocean age, shifted if necessary so youngest age is 1
##
## NOTE: This code uses the rube package to facilitate conditional run-time substitution of various
##  model options. When rube() is called with this model file as its modelFile argument, the
##  following optional arguments must be specified:
##
##  cases       Character vector with up to two elements. One specifies whether or not to estimate an
##                AR(1) coef for alpha.yr ("AR1" or NULL) or not ("IID") while the other specifies 
##                whether indiv-level covariates are included ("Zfull" or NULL) or not ("Z0").
##  varList     A named list with two elements, LC.yr and LC.ind. Each one is a character string
##                resembling the RHS of a standard R formula which specifies the cohort- and indiv-level
##                covars to include, respectively. 
##                YOU SHOULD ALWAYS INCLUDE "- 1" in LC.ind! 
##                To exclude either set of covars from the model entirely, (1) omit the corresponding LC
##                element of varList, (2) use "X0" or "Z0" in cases, and (3) modify the inits 
##                and parameters.to.save.
###############################################################################################################


data
{  
  Y <- max(year)       # number of outmigration years in dataset
  M <- max(age)        # dimension of multinomial probability vector
  #   dimX <- dim(X)       # year x variable matrix of covariates
  #   numX <- dimX[2]      # number of cohort-level (interannual) covariates 
  #   
  #   dimZ <- dim(Z)       # individual x variable matrix of covariates
  #   numZ <- dimZ[2]      # number of individual-level (intra-annual) covariates
}


model
{
  #------------------------------------------------------------------------------
  # PRIORS
  #------------------------------------------------------------------------------
  
  #   a0.int ~ dnorm(0,1e-6)       # intercept (redundant parameterization)
  #   
  #   # Cohort-level regression coefficients
  #   for(k in 1:numX)
  #   {
  #     a0[k] ~ dnorm(0,1e-6)
  #   }
  #   a <- a0/sigma                # map to identified parameters
  
  # Cohort-level regression coefficients
  IFCASE(!X0)
  FOR(a., LC.yr, 0, ? ~ dnorm(0,1e-6))
  FOR(a., LC.yr, , ? <- ?0/sigma)  # map to identified parameters
  ELSECASE
  a.00 ~ dnorm(0,1e-6)
  ENDCASE
  
  #   # Hyper-mean and SD of individual-level regression coefficients
  #   for(k in 1:numZ)
  #   {
  #     mu.b0[k] ~ dnorm(0,1e-6)
  #     sigma.b0[k] ~ dunif(0,10)
  #   }
  #   mu.b <- mu.b0/sigma         # map to identified parameters
  #   sigma.b <- sigma.b0/sigma
  
  # Process SD of individual-level regression coefficients
  IFCASE(!Z0) 
  FOR(sigma.b., LC.ind, 0, ? ~ dunif(0,10))
  FOR(sigma.b., LC.ind, , ? <- ?0/sigma)    # map to identified parameters
  ELSECASE
  #   FOR(sigma.b., LC.ind, 0, ? ~ dunif(0,10))
  FOR(sigma.b., LC.ind, , ? <- 0)  
  ENDCASE
  
  # AR(1) coef
  IFCASE(!IID)
  rho ~ dunif(-1,1)           
  ELSECASE
  rho <- 0                     # no autocorrelation if IID==1
  ENDCASE
  
  sigma.yr0 ~ dunif(0,10)      # process error SD
  sigma.yr <- sigma.yr0/sigma  # map to identified parameters
  sigma ~ dunif(0,10)          # "residual" or individual-level error SD (redundant parameterization)
  
  # Cutpoints
  for(j in 1:(M-1))
  {
    tau.vals[j] ~ dnorm(0,1e-6)
  }
  tau0 <- sort(tau.vals)
  #   tau <- (tau0 - a0.int)/sigma  # map to identified parameters
  tau <- (tau0 - a.00)/sigma    # map to identified parameters
  
  #---------------------------------------------------------------------------------------
  # COHORT-LEVEL PROCESS MODEL
  # Loop over outmigration years, evaluate likelihood of annual state (latent age indicator) 
  # under the process model. Note that initial state is a deterministic function of the
  # predictors, and must be specified separately to initialize the process model.
  #---------------------------------------------------------------------------------------
  
  #   alpha0.mu[1] <- a0.int + sum(X[1,1:numX]*a0)             # expected alpha0 for year 1
  IFCASE(!X0)
  alpha0.mu[1] <- LC(a., LC.yr, 0, 1)                      # expected alpha0 for year 1
  ELSECASE
  alpha0.mu[1] <- a.00
  ENDCASE
  
  alpha0.yr[1] <- alpha0.mu[1]   # alpha0 for year 1
  alpha0.err[1] <- 0             # alpha0 "residual" for year 1
  #   for(k in 1:numZ)
  #   {
  #     b0.yr[1,k] ~ dnorm(mu.b0[k], pow(sigma.b0[k],-2))      # individual-level regression coefs for year 1
  #   }
  
  # individual-level regression coefs for year 1
  IFCASE(!Z0)
  FOR(b., LC.ind, 0, ?[1] ~ dnorm(0,1e-6))  
  FOR(b., LC.ind, , ?[1] <- ?0[1]/sigma)    # map to identified parameters
  ELSECASE
  FOR(b., LC.ind, , ?[1] <- 0)
  ENDCASE
  
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
    #     alpha0.mu[y] <- a0.int + sum(X[y,1:numX]*a0)           # expected alpha0 for year y
    IFCASE(!X0)
    alpha0.mu[y] <- LC(a., LC.yr, 0, y)                    # expected alpha0 for year y
    ELSECASE
    alpha0.mu[y] <- a.00
    ENDCASE
    
    alpha0.yr[y] ~ dnorm(alpha0.mu[y] + rho*alpha0.err[y-1], pow(sigma.yr0,-2))   # alpha0 for year y
    alpha0.err[y] <- alpha0.yr[y] - alpha0.mu[y]           # alpha0 "residual" for year y
    #     for(k in 1:numZ)
    #     {
    #       b0.yr[y,k] ~ dnorm(mu.b0[k], pow(sigma.b0[k],-2))    # individual-level regression coefs for year y
    #     }
    
    # individual-level regression coefs for year y
    IFCASE(!Z0)
    FOR(b., LC.ind, 0, ?[y] ~ dnorm(?[y-1],sigma.?))    
    FOR(b., LC.ind, , ?[y] <- ?0[y]/sigma)                 # map to identified parameters
    ELSECASE
    FOR(b., LC.ind, , ?[y] <- 0)
    ENDCASE
    
    # partition alpha0 using cutpoints
    F.yr[y,1] <- 0
    for(j in 1:(M-1))                                   
    {
      F.yr[y,j+1] <- phi((tau0[j] - alpha0.yr[y])/sigma)
      p.yr[y,j] <- max(min(F.yr[y,j+1] - F.yr[y,j], 1), 1e-10)
    }
    p.yr[y,M] <- max(min(1 - F.yr[y,M], 1), 1e-10)
  }
  
  #   alpha.yr <- (alpha0.yr - a0.int)/sigma                   # map to identified parameters
  alpha.yr <- (alpha0.yr - a.00)/sigma                   # map to identified parameters
  
  #---------------------------------------------------------------------------------------
  # INDIVIDUAL-LEVEL PROCESS MODEL AND LIKELIHOOD
  # Loop over all individual fish. Adjust "expected" state vector for the cohort (i.e.,
  # smolt year) by individual-level covariate effects. Evaluate multinomial likelihood
  # of observed ocean age at adult return.
  #---------------------------------------------------------------------------------------

  # construct vector of ages (for calculating predicted age)
  for(j in 1:M)  
  {
    a1M[j] <- j
  }
  
  for(i in 1:length(age))
  {
    #     alpha0[i] <- alpha0.yr[year[i]] + sum(Z[i,1:numZ]*b0.yr[year[i],1:numZ])  # alpha0 for individual i
    IFCASE(!Z0)
    alpha0[i] <- alpha0.yr[year[i]] + LC(b., LC.ind, 0[year[i]], i) # alpha0 for individual i
    ELSECASE
    alpha0[i] <- alpha0.yr[year[i]]
    ENDCASE
    
    # partition alpha0 using cutpoints
    F[i,1] <- 0
    for(j in 1:(M-1))                                   
    {
      F[i,j+1] <- phi((tau0[j] - alpha0[i])/sigma)
      p[i,j] <- max(min(F[i,j+1] - F[i,j], 1), 1e-10)
    }
    p[i,M] <- max(min(1 - F[i,M], 1), 1e-10)
    
    age[i] ~ dcat(p[i,1:M])                                # likelihood for individual i
    #     cc[i] <- p[i,age[i]]==max(p[i,1:M])                    # correct classification score (0/1)
    pred[i] <- sum(a1M*p[i,1:M])                           # predicted age
    
    #     for(j in 1:i)                                          # tally concordant and discordant pairs of pred and obs age
    #     {
    #       include[(i-1)*length(age) + j] <- (i != j)*(pred[i] != pred[j])
    #       concord[(i-1)*length(age) + j] <- (pred[i] > pred[j])*(age[i] > age[j]) + (pred[i] < pred[j])*(age[i] < age[j])
    #       discord[(i-1)*length(age) + j] <- (pred[i] > pred[j])*(age[i] < age[j]) + (pred[i] < pred[j])*(age[i] > age[j])
    #     }
  }
  #   ccr <- mean(cc)                                          # overall correct classification rate
  #   somersD <- (sum(concord) + sum(discord))/sum(include)   # Somers' D coefficient of concordance
  
  # Got all that?
}

