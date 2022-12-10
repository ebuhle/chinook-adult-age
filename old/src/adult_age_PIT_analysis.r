###############################################################################################################
# FIT HIERARCHICAL MODELS OF SALMON AGE-AT-RETURN TO PIT TAG DATA:
# WILD-ORIGIN SNAKE RIVER SPRING/SUMMER CHINOOK
#
# See HB_adult_age_PIT.r for model description and JAGS specification.
###############################################################################################################

# NOTES/QUESTIONS:
# (1) For the fish tagged in fall, where are the recaps? LGR? BON? Both?
# (2) I replaced missing or zero lengths in the main adult age dataset with NAs. Is this right?
# (3) There are missing length.M and length.R obs in the M-R dataset -- WTF?
# (4) For the main adult age dataset (no M-R), I inferred tagging date from the BON.smolt.date
#     and tag.BON.days columns. Is this right?
# (5) For the first cut of the growth analysis, I did not use the marsh LGR-BON fish because 
#     (i) they are only available for three of the smolt years in the adult age data set, and
#     (ii) they appear to be growing substantially faster than the fall M-R and non-M-R fish,
#          even accounting for size.

library(lubridate)
library(R2jags)
library(rube)
library(Hmisc)
options(device=windows)


#==================================================================
# DATA
#==================================================================

# ADULT AGE AND COVARIATES
srssW.age <- read.table("srssW.age.txt", header=T, sep="\t")
srssW.age <- na.omit(srssW.age[,c("tag","age","year","LGR.smolt.date","LGR.smolt.day","route","rear.type",
                                  "length","LGR.flow","LGR.temp","CUI","CUMCUI","NPGO","ONI","PNI","PDO")])
srssW.age <- srssW.age[srssW.age$length < 165,]  # omit fish >= 165 mm
srssW.age$LGR.smolt.date <- as.Date(as.character(srssW.age$LGR.smolt.date), format="%m/%d/%Y")

# > names(srssW.age)
# [1] "tag"            "age"            "year"           "LGR.smolt.date" "LGR.smolt.day" 
# [6] "route"          "rear.type"      "length"         "LGR.adult.date" "BON.adult.date"
# [11] "LGR.flow"       "LGR.temp"       "CUI"            "CUMCUI"         "NPGO"          
# [16] "ONI"            "PNI"            "PDO"            "SST.buoy29"     "press.buoy29"  

# (Optionally) select a stratified random subset to test-drive code
# indx <- unlist(lapply(unique(srssW.age$year), function(x) sample(which(srssW.age$year==x), 100, replace=F)))
# srssW.age <- srssW.age[indx,]

# Assign global data objects for JAGS
year <- as.integer(srssW.age$year - min(srssW.age$year)) + 1
age <- as.integer(srssW.age$age)
age[age==4] <- 3  # change small number of age 4 obs to "age 3+"
X <- aggregate.data.frame(srssW.age[,c("CUI","PDO","NPGO","ONI","PNI")], by=list(srssW.age$year), mean)[,-1]
X <- scale(X)
Z <- srssW.age[,c("route","length","LGR.smolt.day","LGR.flow","LGR.temp")]
names(Z)[names(Z)=="route"] <- "transport"
names(Z)[names(Z)=="length"] <- "size"
names(Z)[names(Z)=="LGR.smolt.day"] <- "day"
names(Z)[names(Z)=="LGR.flow"] <- "flow"
names(Z)[names(Z)=="LGR.temp"] <- "temp"
Z$transport <- as.numeric(Z$transport=="T")
Z[,-1] <- scale(Z[,-1])


#==================================================================
# AGE-AT-RETURN MODELS
#==================================================================

#------------------------------------------------------------------
# ORDERED PROBIT
# No cohort-level covariates
# No individual-level covariates
# AR(1) errors
#------------------------------------------------------------------

# model-specific data objects
terms.yr <- c()
terms.ind <- c()

# initial values
inits.fn <- function() 
{
  # Note that this redundant-seeming kludge is necessary because when parallel = TRUE, 
  # apparently jags can neither take arguments nor "see" the global environment.
  terms.yr <- c()
  terms.ind <- c()
  
  inits <- list(a.00 = rnorm(1,0,1), rho0 = runif(1,-0.2,0.2), sigma.yr0 = runif(1,0,5),
                sigma = runif(1,0,1), tau.vals = runif(max(age)-1,-0.5,0.5))
  if(!is.null(terms.yr))
  {
    for(k in 1:length(terms.yr))
      inits[[paste("a.", terms.yr[k], "0", sep="")]] <- rnorm(1, 0, 1)
  }
  
  if(!is.null(terms.ind))
  {
    for(k in 1:length(terms.ind))
      inits[[paste("sigma.b.", terms.ind[k], "0", sep="")]] <- runif(1, 0, 1)
  }
  
  return(inits)
}

# Call JAGS via rube to fit the model, with appropriate optional arguments
jags.OP.X0.Z0.AR1 <- rube(model = "HB_OP_adult_age_PIT.r", 
                          data = list(X=as.data.frame(X),Z=as.data.frame(Z),year=year,age=age), 
                          inits = inits.fn, modelFile = "tempModel.r",
                          cases = c("X0","Z0","AR1"),
                          varList=list(LC.yr=NULL, 
                                       LC.ind=NULL),
                          parameters.to.save = c("rho","sigma.yr","tau","alpha.yr","p.yr","pred"),
                          modelCheck = "never", ignore="all", cullData=F, cullInits=F, cullPts=F, 
                          n.chains = 3, n.iter = 20000, n.burnin = 2000, n.thin = 18, 
                          DIC = T, parallel = T)

# Plot and print summaries of results
print(jags.OP.X0.Z0.AR1,limit=c(Inf,Inf),digits=2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:2, list(x=jags.OP.X0.Z0.AR1$sims.array), SIMPLIFY=F)), ask=T)

# Clean up
rm(list=c("terms.yr","terms.ind","inits.fn"))


#------------------------------------------------------------------
# ORDERED PROBIT
# "Full" cohort-level covariates
# No individual-level covariates
# AR(1) errors
#------------------------------------------------------------------

# model-specific data objects
terms.yr <- c("CUI","PDO","NPGO","ONI","PNI")
terms.ind <- c()

# initial values
inits.fn <- function() 
{
  # Note that this redundant-seeming kludge is necessary because when parallel = TRUE, 
  # apparently jags can neither take arguments nor "see" the global environment.
  terms.yr <- c("CUI","PDO","NPGO","ONI","PNI")
  terms.ind <- c()
  
  inits <- list(a.00 = rnorm(1,0,1), rho0 = runif(1,-0.2,0.2), sigma.yr0 = runif(1,0,5),
                sigma = runif(1,0,1), tau.vals = runif(max(age)-1,-0.5,0.5))
  if(!is.null(terms.yr))
  {
    for(k in 1:length(terms.yr))
      inits[[paste("a.", terms.yr[k], "0", sep="")]] <- rnorm(1, 0, 1)
  }
  
  if(!is.null(terms.ind))
  {
    for(k in 1:length(terms.ind))
      inits[[paste("sigma.b.", terms.ind[k], "0", sep="")]] <- runif(1, 0, 1)
  }
  
  return(inits)
}

# Call JAGS via rube to fit the model, with appropriate optional arguments
jags.OP.Xfull.Z0.AR1 <- rube(model = "HB_OP_adult_age_PIT.r", 
                          data = list(X=as.data.frame(X),Z=as.data.frame(Z),year=year,age=age), 
                          inits = inits.fn, modelFile = "tempModel.r",
                          cases = c("Xfull","Z0","AR1"),
                          varList=list(LC.yr=paste(terms.yr, collapse="+"), 
                                       LC.ind=NULL),
                          parameters.to.save = c(paste("a.", terms.yr, sep=""), 
                                                 "rho","sigma.yr","tau","alpha.yr","p.yr","pred"),
                          modelCheck = "never", ignore="all", cullData=F, cullInits=F, cullPts=F, 
                          n.chains = 3, n.iter = 20000, n.burnin = 2000, n.thin = 18, 
                          DIC = T, parallel = T)

# Plot and print summaries of results
print(jags.OP.Xfull.Z0.AR1,limit=c(Inf,Inf),digits=2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:2, list(x=jags.OP.Xfull.Z0.AR1$sims.array), SIMPLIFY=F)), ask=T)

# Clean up
rm(list=c("terms.yr","terms.ind","inits.fn"))


#------------------------------------------------------------------
# ORDERED PROBIT
# No cohort-level covariates
# "Full" individual-level covariates
# AR(1) errors
#------------------------------------------------------------------

# model-specific data objects
terms.yr <- c()
terms.ind <- c("transport","size","day","flow","temp")

# initial values
inits.fn <- function() 
{
  # Note that this redundant-seeming kludge is necessary because when parallel = TRUE, 
  # apparently jags can neither take arguments nor "see" the global environment.
  terms.yr <- c()
  terms.ind <- c("transport","size","day","flow","temp")
  
  inits <- list(a.00 = rnorm(1,0,1), rho0 = runif(1,-0.2,0.2), sigma.yr0 = runif(1,0,5),
                sigma = runif(1,0,1), tau.vals = runif(max(age)-1,-0.5,0.5))
  if(!is.null(terms.yr))
  {
    for(k in 1:length(terms.yr))
      inits[[paste("a.", terms.yr[k], "0", sep="")]] <- rnorm(1, 0, 1)
  }
  
  if(!is.null(terms.ind))
  {
    for(k in 1:length(terms.ind))
      inits[[paste("sigma.b.", terms.ind[k], "0", sep="")]] <- runif(1, 0, 1)
  }
  
  return(inits)
}

# Call JAGS via rube to fit the model, with appropriate optional arguments
jags.OP.X0.Zfull.AR1 <- rube(model = "HB_OP_adult_age_PIT.r", 
                                data = list(X=as.data.frame(X),Z=as.data.frame(Z),year=year,age=age), 
                                inits = inits.fn, modelFile = "tempModel.r",
                                cases = c("X0","Zfull","AR1"),
                                varList=list(LC.yr=NULL, 
                                             LC.ind=paste(paste(terms.ind, collapse="+"), "-1")),
                                parameters.to.save = c(paste("sigma.b.", terms.ind, sep=""), 
                                                       paste("b.", terms.ind, sep=""),
                                                       "rho","sigma.yr","tau","alpha.yr","p.yr","pred"),
                                modelCheck = "never", ignore="all", cullData=F, cullInits=F, cullPts=F, 
                                n.chains = 3, n.iter = 20000, n.burnin = 2000, n.thin = 18, 
                                DIC = T, parallel = T)

# Plot and print summaries of results
print(jags.OP.X0.Zfull.AR1,limit=c(Inf,Inf),digits=2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:2, list(x=jags.OP.X0.Zfull.AR1$sims.array), SIMPLIFY=F)), ask=T)

# Clean up
rm(list=c("terms.yr","terms.ind","inits.fn"))


#------------------------------------------------------------------
# ORDERED PROBIT
# "Full" cohort-level covariates
# "Full" individual-level covariates
# AR(1) errors
#------------------------------------------------------------------

# model-specific data objects
terms.yr <- c("CUI","PDO","NPGO","ONI","PNI")
terms.ind <- c("transport","size","day","flow","temp")

# initial values
inits.fn <- function() 
{
  # Note that this redundant-seeming kludge is necessary because when parallel = TRUE, 
  # apparently jags can neither take arguments nor "see" the global environment.
  terms.yr <- c("CUI","PDO","NPGO","ONI","PNI")
  terms.ind <- c("transport","size","day","flow","temp")
  
  inits <- list(a.00 = rnorm(1,0,1), rho0 = runif(1,-0.2,0.2), sigma.yr0 = runif(1,0,5),
                sigma = runif(1,0,1), tau.vals = runif(max(age)-1,-0.5,0.5))
  if(!is.null(terms.yr))
  {
    for(k in 1:length(terms.yr))
      inits[[paste("a.", terms.yr[k], "0", sep="")]] <- rnorm(1, 0, 1)
  }
  
  if(!is.null(terms.ind))
  {
    for(k in 1:length(terms.ind))
      inits[[paste("sigma.b.", terms.ind[k], "0", sep="")]] <- runif(1, 0, 1)
  }
  
  return(inits)
}

# Call JAGS via rube to fit the model, with appropriate optional arguments
jags.OP.Xfull.Zfull.AR1 <- rube(model = "HB_OP_adult_age_PIT.r", 
                                data = list(X=as.data.frame(X),Z=as.data.frame(Z),year=year,age=age), 
                                inits = inits.fn, modelFile = "tempModel.r",
                                cases = c("Xfull","Zfull","AR1"),
                                varList=list(LC.yr=paste(terms.yr, collapse="+"), 
                                             LC.ind=paste(paste(terms.ind, collapse="+"), "-1")),
                                parameters.to.save = c(paste("a.", terms.yr, sep=""), 
                                                       paste("sigma.b.", terms.ind, sep=""), 
                                                       paste("b.", terms.ind, sep=""),
                                                       "rho","sigma.yr","tau","alpha.yr","p.yr","pred"),
                                modelCheck = "never", ignore="all", cullData=F, cullInits=F, cullPts=F, 
                                n.chains = 3, n.iter = 20000, n.burnin = 2000, n.thin = 18, 
                                DIC = T, parallel = T)

# Plot and print summaries of results
print(jags.OP.Xfull.Zfull.AR1,limit=c(Inf,Inf),digits=2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:2, list(x=jags.OP.Xfull.Zfull.AR1$sims.array), SIMPLIFY=F)), ask=T)

# Clean up
rm(list=c("terms.yr","terms.ind","inits.fn"))


#------------------------------------------------------------------
# Model comparison
#------------------------------------------------------------------

jags.list <- list(jags.OP.X0.Z0.AR1 = jags.OP.X0.Z0.AR1, 
                  jags.OP.Xfull.Z0.AR1 = jags.OP.Xfull.Z0.AR1,
                  jags.OP.X0.Zfull.AR1 = jags.OP.X0.Zfull.AR1,
                  jags.OP.Xfull.Zfull.AR1 = jags.OP.Xfull.Zfull.AR1)

mod.tab <- data.frame(mod.no = 1:length(jags.list), model = names(jags.list), 
                      covars = c("no covariates","cohort","individual","cohort + indiv"),
                      Dbar = unlist(lapply(jags.list, function(m) m$mean$deviance)),
                      pD = unlist(lapply(jags.list, function(m) m$pD)), 
                      DIC = unlist(lapply(jags.list, function(m) m$DIC)),
                      dDIC = NA, wDIC = NA, Dxy.mean = NA, Dxy.025 = NA, Dxy.975 = NA)
row.names(mod.tab) <- NULL
mod.tab$dDIC <- mod.tab$DIC - min(mod.tab$DIC)
mod.tab$wDIC <- exp(-mod.tab$dDIC/2)/sum(exp(-mod.tab$dDIC/2))
for(j in 1:length(jags.list))
{
  Dxy <- apply(jags.list[[j]]$sims.list$pred, 1, function(x) rcorr.cens(x, age)["Dxy"])
  mod.tab$Dxy.mean[j] <- mean(Dxy)
  mod.tab$Dxy.025[j] <- quantile(Dxy, 0.025)
  mod.tab$Dxy.975[j] <- quantile(Dxy, 0.975)
  rm(Dxy)
}
mod.tab <- mod.tab[order(mod.tab$dDIC),]
rm(jags.list)



#------------------------------------------------------------------
# Plots
#------------------------------------------------------------------

###
# Posterior boxpots of cohort-level regresion coef estimates for a given model
###
mod <- jags.OP.Xfull.Zfull.AR1$sims.matrix
amat <- mod[,substring(dimnames(mod)[[2]], 1, 2) == "a."]

dev.new()
par(mar=c(5.1,6,4.1,2.1))
# png(filename="cohort_coefs.png", width=9, height=8, units="in", res=200, type="cairo-png")

bxp.dat <- boxplot(amat, plot=F)
bxp.dat$stats[3,] <- colMeans(amat)
bxp.dat$stats[c(2,4),] <- apply(amat,2,quantile,c(0.05,0.95))
bxp.dat$stats[c(1,5),] <- apply(amat,2,quantile,c(0.025,0.975))
bxp.dat$names <- dimnames(X)[[2]]
bxp(bxp.dat, xlab=expression(paste("Standardized coefficient (", gamma, ")")), ylab="", ylim=range(bxp.dat$stats), 
    show.names=TRUE, horizontal=T, whisklty=1, staplewex=0, boxwex=0.1, outpch="", 
    boxcol="black", boxfill="black", medlty="blank", medpch=16, medcex=1, medcol="white", 
    las=1, cex.axis=1, cex.lab=1.2)
mtext("Cohort-level covariates", side=3, line=1, outer=TRUE, cex=1.5)
abline(v=0, col="gray", lwd=3)
rm(list=c("bxp.dat","mod","amat"))

# dev.off()


###
# Time series of time-varying individual-level coefs
###
mod <- jags.OP.Xfull.Zfull.AR1$sims.list
blist <- mod[names(mod)[substring(names(mod),1,2)=="b."]]

dev.new(width=14,height=7)
# png(filename="indiv_coefs.png", width=14, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,3), mar=c(4.1,5.1,4.1,1.1))
for(i in 1:length(blist))
{
  plot(sort(unique(srssW.age$year)), colMeans(blist[[i]]), pch="", las=1, cex.lab=1.5, cex.axis=1.2,
       ylim=range(sapply(blist,function(x) apply(x,2,quantile,c(0.025,0.975)))),
       xlab=ifelse(i %in% 4:5, "Outmigration year", ""), 
       ylab=ifelse(i %in% c(1,4), expression(paste("Standardized coefficient (", beta[italic(t)], ")")), ""))
  ct <- col2rgb("darkgray")
  ct <- rgb(red=ct[1], green=ct[2], blue=ct[3],alpha=0.7*255, maxColorValue=255)
  polygon(c(sort(unique(srssW.age$year)), rev(sort(unique(srssW.age$year)))), 
            c(apply(blist[[i]],2,quantile,0.025), rev(apply(blist[[i]],2,quantile,0.975))),
            col=ct, border=ct)
  lines(sort(unique(srssW.age$year)), colMeans(blist[[i]]), lwd=3)
  abline(h=0, col="black")
  mtext(paste("(", LETTERS[i], ") ", substring(names(blist)[i], 3), sep=""), side=3, adj=0, line=1)
}
# dev.off()
rm(list=c("mod","blist","ct"))


###
# Paired scatterplots showing posterior correlations among coefs, for a specified year
###

yi <- 2
mod <- jags.OP.Xfull.Zfull.AR1$sims.list
alist <- mod[names(mod)[substring(names(mod),1,2) == "a."]]
blist <- mod[names(mod)[substring(names(mod),1,2) == "b."]]

dev.new(width=10,height=10)
par(oma=c(0,0,2,0))
pairs(cbind(sapply(alist, as.matrix), sapply(blist, function(x) x[,yi])))
mtext(paste("outmigration year", sort(unique(srssW.age$year))[yi]), 3, outer=T, line=1)

# all annual correlations between b.day and b.temp
sapply(1:13, function(i) cor(jags.OP.Xfull.Zfull.AR1$sims.list$b.day[,i], 
                             jags.OP.Xfull.Zfull.AR1$sims.list$b.temp[,i]))


###
# 3-panel Time series plot of cohort-level observed and estimated age frequencies,
# given migration route (I vs. T)
###

mod <- jags.OP.Xfull.Zfull.AR1$sims.list
alphaI <- mod$alpha.yr
alphaT <- alphaI + mod$b.transport
pI <- array(0, c(dim(alphaI), 3))
pI[,,1] <- sapply(1:ncol(alphaI), function(i) pnorm(mod$tau[,1], alphaI[,i], 1))
pI[,,2] <- sapply(1:ncol(alphaI), function(i) pnorm(mod$tau[,2], alphaI[,i], 1)) - pI[,,1]
pI[,,3] <- 1 - pI[,,1] - pI[,,2]
pT <- array(0, c(dim(alphaI), 3))
pT[,,1] <- sapply(1:ncol(alphaT), function(i) pnorm(mod$tau[,1], alphaT[,i], 1))
pT[,,2] <- sapply(1:ncol(alphaT), function(i) pnorm(mod$tau[,2], alphaT[,i], 1)) - pT[,,1]
pT[,,3] <- 1 - pT[,,1] - pT[,,2]

dev.new()
# png(filename="timeseriesIT_3panel.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(3,1), mar=c(4.1,4.1,0.1,2.1))

plot(unique(srssW.age$year), tapply(age==1,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.1*max(apply(pT[,,1],2,quantile,0.975))),
     xlab="", ylab=expression(italic(p[1])), las=1, pch=16, col="orange", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,1]), lwd=2, col="orange")
lines(unique(srssW.age$year), apply(pI[,,1],2,quantile,0.025), col="orange")
lines(unique(srssW.age$year), apply(pI[,,1],2,quantile,0.975), col="orange")
points(unique(srssW.age$year), tapply(age==1,list(year=year, transport=Z$transport),mean)[,"1"],
       pch=1, col="orange", cex=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,1]), lty=3, lwd=2, col="orange")
lines(unique(srssW.age$year), apply(pT[,,1],2,quantile,0.025), lty=3, col="orange")
lines(unique(srssW.age$year), apply(pT[,,1],2,quantile,0.975), lty=3, col="orange")

legend("topleft", legend=c("in-river","transported"), pch=c(16,1), lty=c(1,3), cex=1.2)

plot(unique(srssW.age$year), tapply(age==2,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.1*max(apply(pT[,,2],2,quantile,0.975))),
     xlab="", ylab=expression(italic(p[2])), las=1, pch=16, col="blue", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,2]), lwd=2, col="blue")
lines(unique(srssW.age$year), apply(pI[,,2],2,quantile,0.025), col="blue")
lines(unique(srssW.age$year), apply(pI[,,2],2,quantile,0.975), col="blue")
points(unique(srssW.age$year), tapply(age==2,list(year=year, transport=Z$transport),mean)[,"1"],
       pch=1, col="blue", cex=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,2]), lty=3, lwd=2, col="blue")
lines(unique(srssW.age$year), apply(pT[,,2],2,quantile,0.025), lty=3, col="blue")
lines(unique(srssW.age$year), apply(pT[,,2],2,quantile,0.975), lty=3, col="blue")

plot(unique(srssW.age$year), tapply(age==3,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.1*max(apply(pI[,,3],2,quantile,0.975))),
     xlab="Outmigration year", ylab=expression(italic(p[3])), las=1, pch=16, col="black", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,3]), lwd=2, col="black")
lines(unique(srssW.age$year), apply(pI[,,3],2,quantile,0.025), col="black")
lines(unique(srssW.age$year), apply(pI[,,3],2,quantile,0.975), col="black")
points(unique(srssW.age$year), tapply(age==3,list(year=year, transport=Z$transport),mean)[,"1"],
       pch=1, col="black", cex=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,3]), lty=3, lwd=2, col="black")
lines(unique(srssW.age$year), apply(pT[,,3],2,quantile,0.025), lty=3, col="black")
lines(unique(srssW.age$year), apply(pT[,,3],2,quantile,0.975), lty=3, col="black")

rm(mod);rm(alphaI);rm(pI);rm(alphaT);rm(pT)
# dev.off()


###
# 6-panel Time series plot of cohort-level observed and estimated age frequencies,
# given migration route (I vs. T)
###

mod <- jags.OP.Xfull.Zfull.AR1$sims.list
alphaI <- mod$alpha.yr
alphaT <- alphaI + mod$b.transport
pI <- array(0, c(dim(alphaI), 3))
pI[,,1] <- sapply(1:ncol(alphaI), function(i) pnorm(mod$tau[,1], alphaI[,i], 1))
pI[,,2] <- sapply(1:ncol(alphaI), function(i) pnorm(mod$tau[,2], alphaI[,i], 1)) - pI[,,1]
pI[,,3] <- 1 - pI[,,1] - pI[,,2]
pT <- array(0, c(dim(alphaI), 3))
pT[,,1] <- sapply(1:ncol(alphaT), function(i) pnorm(mod$tau[,1], alphaT[,i], 1))
pT[,,2] <- sapply(1:ncol(alphaT), function(i) pnorm(mod$tau[,2], alphaT[,i], 1)) - pT[,,1]
pT[,,3] <- 1 - pT[,,1] - pT[,,2]

dev.new()
# png(filename="timeseriesIT_6panel.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(3,2), mar=c(4.1,4.1,0.1,2.1), oma=c(0,0,3,0))

plot(unique(srssW.age$year), tapply(age==1,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.05*max(apply(pI[,,1],2,quantile,0.975))),
     xlab="", ylab=expression(italic(p[1])), las=1, pch=16, col="orange", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,1]), lwd=2, col="orange")
lines(unique(srssW.age$year), apply(pI[,,1],2,quantile,0.025), col="orange")
lines(unique(srssW.age$year), apply(pI[,,1],2,quantile,0.975), col="orange")

plot(unique(srssW.age$year), tapply(age==1,list(year=year, transport=Z$transport),mean)[,"1"], 
     ylim=c(0,1.05*max(apply(pI[,,1],2,quantile,0.975))),
     xlab="", ylab="", las=1, pch=16, col="orange", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,1]), lty=1, lwd=2, col="orange")
lines(unique(srssW.age$year), apply(pT[,,1],2,quantile,0.025), lty=1, col="orange")
lines(unique(srssW.age$year), apply(pT[,,1],2,quantile,0.975), lty=1, col="orange")

plot(unique(srssW.age$year), tapply(age==2,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.05*max(apply(pI[,,2],2,quantile,0.975))),
     xlab="", ylab=expression(italic(p[2])), las=1, pch=16, col="blue", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,2]), lwd=2, col="blue")
lines(unique(srssW.age$year), apply(pI[,,2],2,quantile,0.025), col="blue")
lines(unique(srssW.age$year), apply(pI[,,2],2,quantile,0.975), col="blue")

plot(unique(srssW.age$year), tapply(age==2,list(year=year, transport=Z$transport),mean)[,"1"], 
     ylim=c(0,1.05*max(apply(pI[,,2],2,quantile,0.975))),
     xlab="", ylab="", las=1, pch=16, col="blue", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,2]), lty=1, lwd=2, col="blue")
lines(unique(srssW.age$year), apply(pT[,,2],2,quantile,0.025), lty=1, col="blue")
lines(unique(srssW.age$year), apply(pT[,,2],2,quantile,0.975), lty=1, col="blue")

plot(unique(srssW.age$year), tapply(age==3,list(year=year, transport=Z$transport),mean)[,"0"], 
     ylim=c(0,1.05*max(apply(pI[,,3],2,quantile,0.975))),
     xlab="Outmigration year", ylab=expression(italic(p[3])), las=1, pch=16, col="black", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pI[,,3]), lwd=2, col="black")
lines(unique(srssW.age$year), apply(pI[,,3],2,quantile,0.025), col="black")
lines(unique(srssW.age$year), apply(pI[,,3],2,quantile,0.975), col="black")

plot(unique(srssW.age$year), tapply(age==3,list(year=year, transport=Z$transport),mean)[,"1"], 
     ylim=c(0,1.05*max(apply(pI[,,3],2,quantile,0.975))),
     xlab="Outmigration year", ylab="", las=1, pch=16, col="black", cex=1.2, cex.axis=1, cex.lab=1.2)
lines(unique(srssW.age$year), colMeans(pT[,,3]), lty=1, lwd=2, col="black")
lines(unique(srssW.age$year), apply(pT[,,3],2,quantile,0.025), lty=1, col="black")
lines(unique(srssW.age$year), apply(pT[,,3],2,quantile,0.975), lty=1, col="black")

mtext("in-river", side=3, at=0.27, outer=T, line=1, cex=1.2)
mtext("transported", side=3, at=0.77, outer=T, line=1, cex=1.2)

rm(mod);rm(alphaI);rm(pI);rm(alphaT);rm(pT)
# dev.off()


###
# Observed and estimated effects of length for selected years
###
mod <- jags.OP.Xfull.Zfull.AR1$sims.list
yrs <- as.numeric(sort(names(rev(sort(table(srssW.age$year))))[1:4]))
yrs.indx <- sort(rev(order(table(srssW.age$year)))[1:4])

dev.new()
# png(filename="length_curves.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
for(i in 1:4)
{  
  len.crv <- seq(min(Z$size[year==yrs.indx[i]]), max(Z$size[year==yrs.indx[i]]), length=100)
  len.crv.adj <- len.crv*sd(srssW.age$length) + mean(srssW.age$length)
  len.crv <- matrix(rep(len.crv,each=nrow(mod$alpha.yr)), nrow(mod$alpha.yr), length(len.crv))
  b.len <- mod$b.size[,yrs.indx[i]]
  alpha.crv <- matrix(mod$alpha.yr[,yrs.indx[i]], nrow=nrow(mod$alpha.yr), ncol=100) +
    matrix(b.len, nrow=nrow(mod$alpha.yr), ncol=100)*len.crv
  p.crv <- array(0, c(dim(alpha.crv), 3))
  p.crv[,,1] <- sapply(1:ncol(alpha.crv), function(i) pnorm(mod$tau[,1], alpha.crv[,i], 1))
  p.crv[,,2] <- sapply(1:ncol(alpha.crv), function(i) pnorm(mod$tau[,2], alpha.crv[,i], 1)) - p.crv[,,1]
  p.crv[,,3] <- 1 - p.crv[,,1] - p.crv[,,2]
  
  len.bin <- 10*round(srssW.age$length/10)
  p1.obs <- tapply(age[year==yrs.indx[i]]==1, len.bin[year==yrs.indx[i]], mean)
  p2.obs <- tapply(age[year==yrs.indx[i]]==2, len.bin[year==yrs.indx[i]], mean)
  p3.obs <- tapply(age[year==yrs.indx[i]]==3, len.bin[year==yrs.indx[i]], mean)
  
  plot(as.numeric(names(p1.obs)), p1.obs, xlim=range(len.bin[is.element(year, yrs.indx)]), ylim=c(0,1), 
       xlab=ifelse(i>2, "Length at LGR (mm)",""), ylab="p(age)", main=yrs[i], las=1, pch=16, col="orange", cex=1.2, cex.lab=1.2)
  lines(len.crv.adj, colMeans(p.crv[,,1]), lwd=2, col="orange")
  lines(len.crv.adj, apply(p.crv[,,1],2,quantile,0.025), col="orange")
  lines(len.crv.adj, apply(p.crv[,,1],2,quantile,0.975), col="orange")
  
  points(as.numeric(names(p2.obs)), p2.obs, pch=16, cex=1.2, col="blue")
  lines(len.crv.adj, colMeans(p.crv[,,2]), lwd=2, col="blue")
  lines(len.crv.adj, apply(p.crv[,,2],2,quantile,0.025), col="blue")
  lines(len.crv.adj, apply(p.crv[,,2],2,quantile,0.975), col="blue")
  
  points(as.numeric(names(p3.obs)), p3.obs, pch=16, cex=1.2, col="black")
  lines(len.crv.adj, colMeans(p.crv[,,3]), lwd=2, col="black")
  lines(len.crv.adj, apply(p.crv[,,3],2,quantile,0.025), col="black")
  lines(len.crv.adj, apply(p.crv[,,3],2,quantile,0.975), col="black")
  
  if(i==1)
    legend("topright", c("1","2","3"), pch=16, cex=1.2, lty=1, col=c("orange","blue","black"), title="ocean age")
  rm(list=c("len.crv","len.crv.adj","b.len","p.crv","len.bin","p1.obs","p2.obs","p3.obs"))
}
rm(list=c("mod","yrs","yrs.indx"))
# dev.off()


###
# Observed and estimated effects of migration date for selected years
###
mod <- jags.OP.Xfull.Zfull.AR1$sims.list
yrs <- as.numeric(sort(names(rev(sort(table(srssW.age$year))))[1:4]))
yrs.indx <- sort(rev(order(table(srssW.age$year)))[1:4])

dev.new()
# png(filename="timing_curves.png", width=7, height=7, units="in", res=200, type="cairo-png")
par(mfrow=c(2,2), mar=c(4.1,4.1,2.1,2.1))
for(i in 1:4)
{  
  day.crv <- seq(min(Z$day[year==yrs.indx[i]]), max(Z$day[year==yrs.indx[i]]), length=100)
  day.crv.adj <- day.crv*sd(srssW.age$LGR.smolt.day) + mean(srssW.age$LGR.smolt.day)
  day.crv <- matrix(rep(day.crv,each=nrow(mod$alpha.yr)), nrow(mod$alpha.yr), length(day.crv))
  b.day <- mod$b.day[,yrs.indx[i]]
  alpha.crv <- matrix(mod$alpha.yr[,yrs.indx[i]], nrow=nrow(mod$alpha.yr), ncol=100) +
    matrix(b.day, nrow=nrow(mod$alpha.yr), ncol=100)*day.crv
  p.crv <- array(0, c(dim(alpha.crv), 3))
  p.crv[,,1] <- sapply(1:ncol(alpha.crv), function(i) pnorm(mod$tau[,1], alpha.crv[,i], 1))
  p.crv[,,2] <- sapply(1:ncol(alpha.crv), function(i) pnorm(mod$tau[,2], alpha.crv[,i], 1)) - p.crv[,,1]
  p.crv[,,3] <- 1 - p.crv[,,1] - p.crv[,,2]
  
  day.bin <- 10*round(srssW.age$LGR.smolt.day/10)
  p1.obs <- tapply(age[year==yrs.indx[i]]==1, day.bin[year==yrs.indx[i]], mean)
  p2.obs <- tapply(age[year==yrs.indx[i]]==2, day.bin[year==yrs.indx[i]], mean)
  p3.obs <- tapply(age[year==yrs.indx[i]]==3, day.bin[year==yrs.indx[i]], mean)
  
  plot(as.numeric(names(p1.obs)), p1.obs, xlim=range(day.bin[is.element(year, yrs.indx)]), ylim=c(0,1), 
       xlab=ifelse(i>2, "Arrival day at LGR",""), ylab="p(age)", main=yrs[i], las=1, pch=16, col="orange", cex=1.2, cex.lab=1.2)
  lines(day.crv.adj, colMeans(p.crv[,,1]), lwd=2, col="orange")
  lines(day.crv.adj, apply(p.crv[,,1],2,quantile,0.025), col="orange")
  lines(day.crv.adj, apply(p.crv[,,1],2,quantile,0.975), col="orange")
  
  points(as.numeric(names(p2.obs)), p2.obs, pch=16, cex=1.2, col="blue")
  lines(day.crv.adj, colMeans(p.crv[,,2]), lwd=2, col="blue")
  lines(day.crv.adj, apply(p.crv[,,2],2,quantile,0.025), col="blue")
  lines(day.crv.adj, apply(p.crv[,,2],2,quantile,0.975), col="blue")
  
  points(as.numeric(names(p3.obs)), p3.obs, pch=16, cex=1.2, col="black")
  lines(day.crv.adj, colMeans(p.crv[,,3]), lwd=2, col="black")
  lines(day.crv.adj, apply(p.crv[,,3],2,quantile,0.025), col="black")
  lines(day.crv.adj, apply(p.crv[,,3],2,quantile,0.975), col="black")
  
  if(i==1)
    legend("topleft", c("age   ","1","2","3"), horiz=T, pch=c(-1,16,16,16), cex=0.8, lty=c(0,1,1,1), col=c("orange","blue","black"))
  rm(list=c("day.crv","day.crv.adj","b.day","p.crv","day.bin","p1.obs","p2.obs","p3.obs"))
}
rm(list=c("mod","yrs","yrs.indx"))
# dev.off()





