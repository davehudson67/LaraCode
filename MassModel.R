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

## create 3D arrays for the data... individual, year, transect
torts <- droplevels(torts)
torts$Transect <- as.numeric(as.factor(torts$TransectCode))
torts$Season <- as.numeric(as.factor(torts$season))
torts$Year <- torts$Year - min(torts$Year - 1)

#torts$Section <- as.numeric(torts$SectionCode)

# Create a mapping of SectionCode to Site for each TransectCode
section_to_site <- torts %>%
  group_by(TransectCode) %>%
  mutate(Site = dense_rank(SectionCode)) %>%
  select(TransectCode, SectionCode, Site) %>%
  distinct()

# Join the mapping back to the original dataframe
torts <- torts %>%
  left_join(section_to_site, by = c("TransectCode", "SectionCode"))

nTransects <- max(torts$Transect) # number of transects
nYears <- length(unique(torts$Year)) 
nSites <- aggregate(SectionCode ~ Transect, data = torts, FUN = function(x) length(unique(x)))
nSites <- nSites$SectionCode
unique_transects <- unique(torts$TransectCode)
unique_years <- unique(torts$Year)
maxSites <- max(nSites)

# Initialize visit_status array
site_visited <- array(0, dim = c(max(torts$Transect), maxSites, nYears, nseason))

# Populate the site_visited array
for (i in 1:nrow(torts)) {
  site_visited[torts$Transect[i], torts$Site[i], torts$Year[i], torts$Season[i]] <- 1
}

# Max observations by site
torts %>%
  group_by(Year, Season, TransectCode, SectionCode) %>%
  summarize(NumObservations = n()) %>%
  ungroup() %>%
  summarize(MaxObservations = max(NumObservations))

Nind <- 80  # Number of augmented and real individuals per site/season/year
nSitesTot <- nrow(VegTypeM)

# Initialize the observation array
y <- array(0, dim = c(nTransects, maxSites, nYears, nseason, Nind))

# Populate the observation array
for (i in 1:nrow(torts)) {
  transect_index <- which(unique_transects == torts$TransectCode[i])
  site_index <- as.numeric(torts$SectionCode[i])  # Adjust as necessary
  year_index <- which(unique_years == torts$Year[i])
  season_index <- ifelse(torts$Season[i] == "Wet", 1, 2)
  
  if (transect_index <= nTransects && site_index <= maxSites && year_index <= nYears && season_index <= 2) {
    id_index <- as.numeric(torts$ID[i])
    if (id_index <= Nind) {
      y[transect_index, site_index, year_index, season_index, id_index] <- 1
    }
  }
}

# Initialize the distances array
x <- array(0, dim = c(nTransects, maxSites, nYears, nseason, Nind))

# Populate the distances array
for (i in 1:nrow(torts)) {
  transect_index <- which(unique_transects == torts$TransectCode[i])
  site_index <- as.numeric(torts$SectionCode[i])
  year_index <- which(unique_years == torts$Year[i])
  season_index <- ifelse(torts$Season[i] == "Wet", 1, 2)
  id_index <- as.numeric(torts$ID[i])
  
  if (transect_index <= nTransects && site_index <= maxSites && year_index <= nYears && season_index <= 2 && id_index <= Nind) {
    x[transect_index, site_index, year_index, season_index, id_index] <- torts$DistFromTransect[i]
  }
}

x[x == 0] <- NA

# Define the NIMBLE model
model <- nimbleCode({
  
  # Priors
  alphaV0 ~ dunif(0, 10)
  for(v in 1:nVegTypes){
    alphaV1[v] ~ dnorm(0, 1)
  }
  for(h in 1:nSeasons){
    betaS1[h] ~ dnorm(0, 1)
  }
  for(j in 1:nTransects){
    psi[j] ~ dunif(0, 1)
  }
  
  # Observation model - site level vegetation effect and season effect
  for (s in 1:nSitesTot) {
    for (se in 1:nSeasons) {
      log(sigma[s, se]) <- alphaV0 + sum(alphaV1[1:nVegTypes] * VegTypeM[s, 1:nVegTypes]) + betaS1[se]
    }
  }

  # Likelihood
  for (t in 1:nTransects) {
    for (s in 1:nSites[t]) {
      for (yr in 1:nYears) {
        for (se in 1:nSeasons) {
          for (k in site_visited[t, s, yr, se]) {
            
            # Data augmentation at the site level
            for (i in 1:(nind[t, s, yr, se])) {
              z[t, s, yr, se, i] ~ dbern(psi[t])            # DA variable for each augmented individual
              x[t, s, yr, se, i] ~ dunif(0, B)              # Distribution of distances
              
              logp[t, s, yr, se, i] <- -((x[t, s, yr, se, i] * x[t, s, yr, se, i]) / (2 * sigma[s, se] * sigma[s, se])) # Detection function
              p[t, s, yr, se, i] <- exp(logp[t, s, yr, se, i])
              
              mu[t, s, yr, se, i] <- z[t, s, yr, se, i] * p[t, s, yr, se, i]
              y[t, s, yr, se, i] ~ dbern(mu[t, s, yr, se, i])             # Observation model
            }
          }
        }
      }
    }
  }
  
  # Derived pop estimates
  for (t in 1:nTransects) {
    for (yr in 1:nYears) {
      for (se in 1:nSeasons) {
        PopEstimates[t, yr, se] <- sum(z[t, 1:nSitesTot, yr, se, 1:Nind])
          }
        }
    }
})

## Constants
consts <- list(nTransects = nTransects,
               nSites = nSites,
               nSitesTot = nrow(VegTypeM),
               nYears = nYears,
               nVegTypes = nVegTypes,
               nSeasons = 2,
               nind = array(Nind, dim = c(nTransects, maxSites, nYears, nseason)),
               B = B,
               site_visited = site_visited,
               VegTypeM = VegTypeM,
               Nind = Nind,
               nSitesTot = nSitesTot)

## InitialValues
inits <- list(alphaV0 = 0,
              alphaV1 = rnorm(nVegTypes, 0, 1),
              betaS1 = rnorm(nseason, 0, 1),
              psi = runif(nTransects, 0, 1)
)

## Data
data <- list(x = x,
             y = y)

# Build model
model <- nimbleModel(model, constants = consts, data = data, inits = inits)

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("PopEstimates", "psi", "alphaV0", "alphaV1", "betaS1"))
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


