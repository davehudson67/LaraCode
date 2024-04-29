# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(abind)
library(lubridate)
library(bayesplot)

## load formatted data
load("tortsReady.RData")
rm(doubleSamp)
rm(effort)
rm(section_codes_in_all_occasions)
rm(tortsUK)
rm(Unknown)
rm(x)
rm(y)
rm(site)
rm(nind)
rm(nyears)
rm(nsites)
#torts <- filter(torts, Year >= 2021)

head(torts)

tortsW2018 <- torts %>%
  filter(Year == 2018)

x <- tortsW2018$DistFromTransect
B <- 60 # Strip half-width. Larger than max distance
nind <- length(x)

################################################################################

# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 500                        # Augment by 500
y <- c(rep(1, nind), rep(0, nz))     # Augmented detection indicator y
site <- c(tortsW2018$SectionCode, rep(NA, nz)) # Augmented site indicator, 
d <- c(tortsW2018$DistFromTransect, rep(NA,nz))     # Augmented distance data (with NAs)
nsites = max(site, na.rm = TRUE) 

# "model1.txt" in book code
Section8p5p2_code <- nimbleCode( {
  # Prior distributions
  beta0 ~ dunif(0,100)   # Intercept of lambda-habitat regression
  #beta1 ~ dunif(-10,10)   # Slope of log(lambda) on habitat
  alpha0 ~ dunif(-10,10)  # Intercept of log(sigma) (half-normal scale)
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dunif(-10,10)  # Slope of log(sigma) on wind
  }
  
  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind+nz)
  
  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Det. function
    mu[i] <- z[i]* p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    #site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
  }
  
  
  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 #+ beta1*habitat[s] # Linear model abundance
    #site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    # Linear model for detection 
    log(sigma[s]) <- alpha0 + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
      alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
      alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind+nz)])
  #area <- nsites*1*2*B   # Unit length == 1, half-width = B
  #D <- Ntotal/area 
}
)

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, VegTypeM = VegTypeM, nVegTypes = nVegTypes, B=B, nind=nind, nz=nz, y=y, d=d, site=site) )

# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function(){list(beta0=0, 
                         #beta1=0, 
                         alpha0 = 1, 
                         alpha = rep(0, nVegTypes),
                         z=zst,
                         ## added site inits to avoid nimble warnings
                         #site = c(rep(NA, nind),
                         #          sample(1:nsites, nz, replace = TRUE))
)}

# Parameters to save
params <- c("alpha0", "alpha", "beta0", "psi", "Ntotal")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Call nimble from R
out1 <- nimbleMCMC(code = Section8p5p2_code, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)