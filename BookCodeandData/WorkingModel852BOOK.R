# 8.5.2 Bayesian HDS using data augmentation
# ------------------------------------------------------------------------
# Recreate line transect data set
# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(mcmcplots)
library(lubridate)
library(AHMbook)

set.seed(1234)

tmp <- simHDS()                  # Line transect (default)
attach(tmp)

# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 500                        # Augment by 500
nind <- nrow(data)
y <- c(data[,2], rep(0, nz))     # Augmented detection indicator y
site <- c(data[,1], rep(NA, nz)) # Augmented site indicator, unknown (i.e., NA) for augmented inds.
d <- c(data[,5], rep(NA,nz))     # Augmented distance data (with NAs)
nsites = 100 

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, habitat=habitat, wind=wind, B=B, nind=nind, nz=nz, y=y, d=d, site=site) )
win.data$site                    # unknown site cov. for augmented inds.

# BUGS model for line transect HDS
# "model1.txt" in book code
Section8p5p2_code <- nimbleCode( {
  # Prior distributions
  beta0 ~ dunif(-10,10)   # Intercept of lambda-habitat regression
  beta1 ~ dunif(-10,10)   # Slope of log(lambda) on habitat
  alpha0 ~ dunif(-10,10)  # Intercept of log(sigma) (half-normal scale)
  alpha1 ~ dunif(-10,10)  # Slope of log(sigma) on wind
  
  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind+nz)
  
  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Det. function
    mu[i] <- z[i]* p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
  }
  
  
  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta1*habitat[s] # Linear model abundance
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    # Linear model for detection 
    log(sigma[s]) <- alpha0 + alpha1*wind[s]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind+nz)])
  area <- nsites*1*2*B   # Unit length == 1, half-width = B
  D <- Ntotal/area 
}
)

# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function(){list(beta0=0, 
                         beta1=0, 
                         alpha0=0, 
                         alpha1=0, 
                         z=zst,
                         ## added site inits to avoid nimble warnings
                         site = c(rep(NA, nind),
                                  sample(1:nsites, nz, replace = TRUE))
)}

# Parameters to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "D")

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

colnames(out1)
require(mcmcplots)
mcmcplot(out1)

Section8p5p2_compare <- compareMCMCs(
  modelInfo = list(
    code = Section8p5p2_code,
    data = win.data,
    inits = inits()
  ),
  MCMCs = c('nimble', 'jags'),
  monitors = params,
  summary = FALSE,
  burnin = nb,
  niter = 2*ni # Larger sample for more accurate ESS estimates
)

make_MCMC_comparison_pages(Section8p5p2_compare, modelNames = "Section8p5p2", dir = outputDirectory)

browseURL(file.path(outputDirectory, "Section8p5p2.html"))

# Prepare data
delta <- 0.1 # width of distance bins for approx.
midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
dclass <- d %/% delta + 1       # convert distances to cat. distances
nD <- length(midpt)             # Number of distance intervals


# Bundle and summarize data set
str( win.data <- list (y=y, dclass=dclass, site=site, midpt=midpt, delta=delta, B=B, nind=nind, nz=nz, nsites=nsites, nD=nD, habitat=habitat, wind=wind) )

# BUGS model specification for line-transect HDS (NOT point transects!)
# "model2.txt" in book code
Section8p5p2_code.line <- nimbleCode({ 
  # Prior distributions
  alpha0 ~ dunif(-10,10)
  alpha1 ~ dunif(-10,10)
  beta0 ~ dunif(-10,10)
  beta1 ~ dunif(-10,10)
  
  psi <- sum(lambda[1:nsites])/(nind+nz)     # psi is a derived parameter
  
  for(i in 1:(nind+nz)){             # Loop over individuals
    z[i] ~ dbern(psi)               # DA variables
    dclass[i] ~ dcat(pi[site[i], 1:30])  # Population distribution of dist class
    mu[i] <- z[i] * p[site[i],dclass[i]] # p depends on site AND dist class
    y[i] ~ dbern(mu[i])             # Basic Bernoulli response in DS model
    site[i] ~ dcat(site.probs[1:nsites]) # Site membership of inds
  }
  
  for(s in 1:nsites){                # Loop over sites
    # Construct cell probabilities for nG cells
    for(g in 1:nD){                    # midpt = mid point of each cell
      log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
      pi[s,g] <- delta/B              # probability of x per interval 
      f[s,g] <- p[s,g]*pi[s,g]        # pdf of observed distances
    }
    
    # not necessary   N[s]~dpois(lambda[s]) except for prediction
    N[s] ~ dpois(lambda[s])        # predict abundance at each site
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model for N
    site.probs[s] <- lambda[s]/sum(lambda[1:nsites])
    log(sigma[s]) <- alpha0 + alpha1*wind[s] # linear model for sigma
  }
  
  # Derived parameter
  Ntotal <- sum(z[1:nz])   # Also sum(N[]) which is size of a new population
  area<- nsites*1*2*B  # Unit length == 1, half-width = B
  D<- Ntotal/area 
}
)


# Inits
zst <- c(rep(1, sum(y)), rep(0, nz))
inits <- function(){list (alpha0=0, 
                          alpha1=0, 
                          beta0=0, 
                          beta1=0, 
                          z=zst,
                          ## added site inits to avoid nimble warnings
                          site = c(rep(NA, nind),
                                   sample(1:nsites, size = nz, replace = TRUE)),
                          ## added dclass inits to avoid nimble warnings
                          dclass = c(rep(NA, nind),
                                     sample(1:nD, size = nz, replace = TRUE))
) }

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run nimble from R
out2 <- nimbleMCMC(code = Section8p5p2_code.line, 
                   constants = win.data, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)

colnames(out2)
mcmcplot(out2)