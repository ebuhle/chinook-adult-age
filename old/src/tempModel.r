

data
{
  Y <- max(year)
  M <- max(age)
}


model
{


  a.00 ~ dnorm(0, 0.000001)
  a.CUI0 ~ dnorm(0, 0.000001)
  a.PDO0 ~ dnorm(0, 0.000001)
  a.NPGO0 ~ dnorm(0, 0.000001)
  a.ONI0 ~ dnorm(0, 0.000001)
  a.PNI0 ~ dnorm(0, 0.000001)
  a.0 <- a.00/sigma
  a.CUI <- a.CUI0/sigma
  a.PDO <- a.PDO0/sigma
  a.NPGO <- a.NPGO0/sigma
  a.ONI <- a.ONI0/sigma
  a.PNI <- a.PNI0/sigma


  sigma.b.transport0 ~ dunif(0, 10)
  sigma.b.size0 ~ dunif(0, 10)
  sigma.b.day0 ~ dunif(0, 10)
  sigma.b.flow0 ~ dunif(0, 10)
  sigma.b.temp0 ~ dunif(0, 10)
  sigma.b.transport <- sigma.b.transport0/sigma
  sigma.b.size <- sigma.b.size0/sigma
  sigma.b.day <- sigma.b.day0/sigma
  sigma.b.flow <- sigma.b.flow0/sigma
  sigma.b.temp <- sigma.b.temp0/sigma

  rho ~ dunif(-1,1)

  sigma.yr0 ~ dunif(0, 10)
  sigma.yr <- sigma.yr0/sigma
  sigma ~ dunif(0, 10)

  for(j in 1:(M-1))
  {
    tau.vals[j] ~ dnorm(0, 0.000001)
  }
  tau0 <- sort(tau.vals)
  tau <- (tau0 - a.00)/sigma


  alpha0.mu[1] <- a.00 + a.CUI0*CUI[1] + a.PDO0*PDO[1] + a.NPGO0*NPGO[1] + a.ONI0*ONI[1] + a.PNI0*PNI[1]

  alpha0.yr[1] <- alpha0.mu[1]
  alpha0.err[1] <- 0

  b.transport0[1] ~ dnorm(0, 0.000001)
  b.size0[1] ~ dnorm(0, 0.000001)
  b.day0[1] ~ dnorm(0, 0.000001)
  b.flow0[1] ~ dnorm(0, 0.000001)
  b.temp0[1] ~ dnorm(0, 0.000001)
  b.transport[1] <- b.transport0[1]/sigma
  b.size[1] <- b.size0[1]/sigma
  b.day[1] <- b.day0[1]/sigma
  b.flow[1] <- b.flow0[1]/sigma
  b.temp[1] <- b.temp0[1]/sigma

  F.yr[1,1] <- 0
  for(j in 1:(M-1))
  {
    F.yr[1,j+1] <- phi((tau0[j] - alpha0.yr[1])/sigma)
    p.yr[1,j] <- max(min(F.yr[1,j+1] - F.yr[1,j], 1), 1e-10)
  }
  p.yr[1,M] <- max(min(1 - F.yr[1,M], 1), 1e-10)

  for(y in 2:Y)
  {
    alpha0.mu[y] <- a.00 + a.CUI0*CUI[y] + a.PDO0*PDO[y] + a.NPGO0*NPGO[y] + a.ONI0*ONI[y] + a.PNI0*PNI[y]

    alpha0.yr[y] ~ dnorm(alpha0.mu[y] + rho*alpha0.err[y-1],  pow(sigma.yr0,-2))
    alpha0.err[y] <- alpha0.yr[y] - alpha0.mu[y]

    b.transport0[y] ~ dnorm(b.transport0[y-1], sigma.b.transport0)
    b.size0[y] ~ dnorm(b.size0[y-1], sigma.b.size0)
    b.day0[y] ~ dnorm(b.day0[y-1], sigma.b.day0)
    b.flow0[y] ~ dnorm(b.flow0[y-1], sigma.b.flow0)
    b.temp0[y] ~ dnorm(b.temp0[y-1], sigma.b.temp0)
    b.transport[y] <- b.transport0[y]/sigma
    b.size[y] <- b.size0[y]/sigma
    b.day[y] <- b.day0[y]/sigma
    b.flow[y] <- b.flow0[y]/sigma
    b.temp[y] <- b.temp0[y]/sigma

    F.yr[y,1] <- 0
    for(j in 1:(M-1))
    {
      F.yr[y,j+1] <- phi((tau0[j] - alpha0.yr[y])/sigma)
      p.yr[y,j] <- max(min(F.yr[y,j+1] - F.yr[y,j], 1), 1e-10)
    }
    p.yr[y,M] <- max(min(1 - F.yr[y,M], 1), 1e-10)
  }

  alpha.yr <- (alpha0.yr - a.00)/sigma


  for(j in 1:M)
  {
    a1M[j] <- j
  }

  for(i in 1:length(age))
  {
    alpha0[i] <- alpha0.yr[year[i]] + b.transport0[year[i]]*transport[i] + b.size0[year[i]]*size[i] + b.day0[year[i]]*day[i] + b.flow0[year[i]]*flow[i] + b.temp0[year[i]]*temp[i]

    F[i,1] <- 0
    for(j in 1:(M-1))
    {
      F[i,j+1] <- phi((tau0[j] - alpha0[i])/sigma)
      p[i,j] <- max(min(F[i,j+1] - F[i,j], 1), 1e-10)
    }
    p[i,M] <- max(min(1 - F[i,M], 1), 1e-10)

    age[i] ~ dcat(p[i,1:M])
    pred[i] <- sum(a1M*p[i,1:M])

  }

}

