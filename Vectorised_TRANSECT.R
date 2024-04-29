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

## create 3D arrays for the data... individual, year, transect
torts <- droplevels(torts)
torts$Transect <- as.numeric(as.factor(torts$TransectCode))
torts$Season <- as.numeric(as.factor(torts$season))
torts$Year <- torts$Year - min(torts$Year - 1)
torts$Site <- as.numeric(torts$SectionCode)

## extract unique combinations of Year, Season, Transect, Site
unique_combinations <- torts %>%
  distinct(Year, Season, Transect, Site)

## max number of observations at any one site in one period
torts %>%
  group_by(Site, occasion) %>%
  summarize(MaxObservations = n()) %>%
  arrange(desc(MaxObservations))

## Define the number of individuals to augment for site/occasion combination
N_aug <- 100

## Create Augmented data
AugData <- unique_combinations[rep(seq_len(nrow(unique_combinations)), each = N_aug), ]

AugData <- AugData %>%
  mutate(DistFromTransect = NA) %>%
  mutate(y = 0)

TortsA <- torts %>%
  select(Year, Transect, Season, Site, DistFromTransect) %>%
  mutate(y = 1) %>%
  full_join(AugData)

## More data needed for model
nTransects <- max(torts$Transect) # number of transects
nYears <- length(unique(torts$Year)) 
nSites <- length(unique(torts$Site))
maxSites <- max(nSites)

# Initialize the observation array
y <- c(TortsA$y)

# Initialize the distances array
x <- c(TortsA$DistFromTransect)

# Define the NIMBLE model
code <- nimbleCode({
  
  # Priors
  sigma ~ dunif(0, 100)
  
  psi <- sum(lambda[1:nTransects]) / nind
  
  # Likelihood
  for(i in 1:(nind)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    
    # Observation model
    logp[i] <- -((x[i] * x[i]) / (2 * sigma * sigma)) # Half-normal detection fct.
    p[i] <- exp(logp[i])
    
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
    
  }
  
  for(t in 1:nTransects){
    N[t] ~ dpois(lambda[t])
    log(lambda[t]) <- beta0 + beta1[t]
    transect.probs <- lambda[t] / sum(lambda[1:nTransects])
  }
  
  # Derived parameters
  Ntotal <- sum(z[1:nind])
  
})

## Constants
consts <- list(nTransects = nTransects,
               nind = length(y),
               Transect = TortsA$Transect,
               site = TortsA$Site,
               season = TortsA$Season,
               B = B)

## InitialValues
inits <- list(alpha0 = rep(0, 2),
              alpha1 = rnorm(nVegTypes, 0, 1),
              psi = runif(nTransects, 0, 1)
)

## Data
data <- list(x = x,
             y = y)

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("psi", "alpha0", "alpha1", "z"))
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 10000, 
                           nburnin = 2000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run$summary
library(mcmcplots)
mcmc_trace(run$samples, pars = c("psi", "alpha0", "alpha1"))
mcmcplot(run$samples, parms = c("psi", "alpha0", "alpha1"))