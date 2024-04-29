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

torts2019 <- torts %>%
  filter(Year == 2019)

# Split data by season
torts2019_se <- split(torts2019, torts2019$season)

data_aug <- torts2019 %>%
  group_by(season) %>%
  do({
    season_data <- .
    sites_visited <- levels(season_data$SectionCode)[table(season_data$SectionCode) > 0]
    nz <- length(sites_visited) * 20  # Augment by 20 at each site
    aug_site <- rep(sites_visited, each = 20)
    aug_site_codes <- factor(aug_site, levels = levels(season_data$SectionCode))
    aug_dist <- rep(NA, nz)
    
    # Ensure season column is included for augmented data, matching the current group's season
    augmented_season <- rep(unique(season_data$season), times = nz)
    
    # Combine the original and augmented data for the season, now including season column
    data.frame(
      season = c(season_data$season, augmented_season),
      SectionCode = factor(c(season_data$SectionCode, aug_site_codes), levels = levels(torts2019$SectionCode)),
      DistFromTransect = c(season_data$DistFromTransect, aug_dist),
      y = c(rep(1, nrow(season_data)), rep(0, nz))  # Distinguish augmented data
    )
  }) %>%
  ungroup()

## Extract data
d <- data_aug$DistFromTransect
site <- as.numeric(data_aug$SectionCode)
nsites = max(site, na.rm = TRUE) 
y <- data_aug$y
nind <- sum(y)
nz <- length(y) - nind
season <- as.numeric(as.factor(data_aug$season))

################################################################################

code <- nimbleCode({
  # Prior distributions
  beta0 ~ dunif(0,100)   # Intercept of lambda
  
  for(k in 1:nseason){
  alpha0[k] ~ dunif(-10,10)  # Intercept of log(sigma) for each season
  }
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dunif(-10,10)  # Slope of log(sigma) on vegetation type
  }
  
  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind+nz)
  
  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){# i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i] * d[i] / (2 * sigma[site[i], season[i]] * sigma[site[i], season[i]])) # Det. function
    mu[i] <- z[i] * p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
  }
  
  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0

    # Linear model for detection
    for(seas in 1:nseason){
    log(sigma[s, seas]) <- alpha0[seas] + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
      alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
      alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind+nz)])
}
  })

# Bundle and summarize data set
data <- list(y=y, d=d)

consts <- list(nsites=nsites, VegTypeM = VegTypeM, nVegTypes = nVegTypes, B=B, nind=nind, nz=nz, site=site, season= season, nseason=nseason)

# Inits
zst <- y
inits <- function(){list(beta0=0, 
                         #beta1=0, 
                         alpha0 = rep(1, 2), 
                         alpha = rep(0, nVegTypes),
                         sigma = matrix(0, nsites, nseason),
                         z=zst)}

# Parameters to save
params <- c("alpha0", "alpha", "beta0", "psi", "Ntotal", "N")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
                   data = data,
                   constants = consts, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   samplesAsCodaMCMC = TRUE)


mcmcplot(out1)
