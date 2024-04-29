# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(abind)
library(lubridate)
library(bayesplot)
library(mcmcplots)

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

## 2019
torts2019 <- torts %>%
  filter(Year == 2019) %>%
  droplevels()

b <- 60

# Split data by season
torts2019_dry <- torts2019 %>%
  filter(season == "Dry") %>%
  droplevels() %>%
  mutate(SectionCode = as.numeric(SectionCode)) %>%
  mutate(TransectCode = as.factor(TransectCode)) %>%
  droplevels() %>%
  mutate(TransectCode = as.numeric(TransectCode))
max(torts2019_dry$SectionCode)

torts2019_wet <- torts2019 %>%
  filter(season == "Wet") %>%
  droplevels() %>%
  mutate(SectionCode = as.numeric(SectionCode)) %>%
  mutate(TransectCode = as.factor(TransectCode)) %>%
  droplevels() %>%
  mutate(TransectCode = as.numeric(TransectCode))
max(torts2019_wet$SectionCode)

sites_season <- c(length(unique(torts2019_dry$SectionCode)), length(unique(torts2019_wet$SectionCode)))

data_aug <- bind_rows(torts2019_dry, torts2019_wet) %>%
  group_by(season) %>%
  do({
    season_data <- .
    # Directly work with numeric values to find unique visited sites
    sites_visited <- unique(season_data$SectionCode)
    nz <- length(sites_visited) * 20  # Augment by 20 at each site
    aug_site <- rep(sites_visited, each = 20)
    aug_dist <- rep(NA, nz)
    
    # Ensure season column is included for augmented data, matching the current group's season
    augmented_season <- rep(unique(season_data$season), times = nz)
    
    # Combine the original and augmented data for the season
    data.frame(
      season = c(season_data$season, augmented_season),
      SectionCode = c(season_data$SectionCode, aug_site), # Keep as numeric
      DistFromTransect = c(season_data$DistFromTransect, aug_dist),
      y = c(rep(1, nrow(season_data)), rep(0, nz))  # Distinguish augmented data
    )
  }) %>%
  ungroup()

## Extract data
dd <- data_aug$DistFromTransect[data_aug$season == "Dry"]
dw <- data_aug$DistFromTransect[data_aug$season == "Wet"]
add_dd <- rep(0, times = length(dw) - length(dd))
dd <- c(dd, add_dd)
d  <- array(c(dd, dw), dim = c(2, length(dd)))

sited <- as.numeric(data_aug$SectionCode[data_aug$season == "Dry"])
sitew <- as.numeric(data_aug$SectionCode[data_aug$season == "Wet"])
add_sd <- rep(1, times = length(sitew) - length(sited))
sited <- c(sited, add_sd)
site <- array(c(sited, sitew), dim = c(2, length(sited)))
nsitesMax <- max(site)

yd <- data_aug$y[data_aug$season == "Dry"]
yw <- data_aug$y[data_aug$season == "Wet"]
add_yd <- rep(0, times = length(yw) - length(yd))
yd <- c(yd, add_yd)
y <- array(c(yd, yw), dim = c(2, length(yd)))

nind <- c(sum(data_aug$y[data_aug$season == "Dry"]), sum(data_aug$y[data_aug$season == "Wet"]))
nz <- c(6517 - nind[1], 6517 - nind[2])

## Sort Veg Type Arrays
VegTypeM <- torts2019_dry %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM_dry <- model.matrix(~ -1 + VegType, data = VegTypeM)

VegTypeM <- torts2019_wet %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM_wet <- model.matrix(~ -1 + VegType, data = VegTypeM)

rowstoadd <- matrix(NA, nrow = 4, ncol = 8)
VegTypeM_wet <- rbind(VegTypeM_wet, rowstoadd)

VegTypeA <- array(c(VegTypeM_dry, VegTypeM_wet), dim = c(263, 8, 2))

## 2020
torts2020 <- torts %>%
  filter(Year == 2020) %>%
  droplevels()

# Split data by season
torts2020_dry <- torts2020 %>%
  filter(season == "Dry") %>%
  droplevels() %>%
  mutate(SectionCode = as.numeric(SectionCode)) %>%
  mutate(TransectCode = as.factor(TransectCode)) %>%
  droplevels() %>%
  mutate(TransectCode = as.numeric(TransectCode))
max(torts2020_dry$SectionCode)

torts2020_wet <- torts2020 %>%
  filter(season == "Wet") %>%
  droplevels() %>%
  mutate(SectionCode = as.numeric(SectionCode)) %>%
  mutate(TransectCode = as.factor(TransectCode)) %>%
  droplevels() %>%
  mutate(TransectCode = as.numeric(TransectCode))
max(torts2020_wet$SectionCode)

sites_season <- array(c(sites_season, length(unique(torts2020_dry$SectionCode)), length(unique(torts2020_wet$SectionCode))), dim = c(2,2))

data_aug <- bind_rows(torts2020_dry, torts2020_wet) %>%
  group_by(season) %>%
  do({
    season_data <- .
    # Directly work with numeric values to find unique visited sites
    sites_visited <- unique(season_data$SectionCode)
    nz <- length(sites_visited) * 20  # Augment by 20 at each site
    aug_site <- rep(sites_visited, each = 20)
    aug_dist <- rep(NA, nz)
    
    # Ensure season column is included for augmented data, matching the current group's season
    augmented_season <- rep(unique(season_data$season), times = nz)
    
    # Combine the original and augmented data for the season
    data.frame(
      season = c(season_data$season, augmented_season),
      SectionCode = c(season_data$SectionCode, aug_site), # Keep as numeric
      DistFromTransect = c(season_data$DistFromTransect, aug_dist),
      y = c(rep(1, nrow(season_data)), rep(0, nz))  # Distinguish augmented data
    )
  }) %>%
  ungroup()

## Extract data
dd <- data_aug$DistFromTransect[data_aug$season == "Dry"]
dw <- data_aug$DistFromTransect[data_aug$season == "Wet"]
add_dd <- rep(0, times = length(dw) - length(dd))
dd <- c(dd, add_dd)
d  <- array(c(dd, dw), dim = c(2, length(dd)))

sited <- as.numeric(data_aug$SectionCode[data_aug$season == "Dry"])
sitew <- as.numeric(data_aug$SectionCode[data_aug$season == "Wet"])
add_sd <- rep(1, times = length(sitew) - length(sited))
sited <- c(sited, add_sd)
site <- array(c(sited, sitew), dim = c(2, length(sited)))
nsitesMax <- max(site)

yd <- data_aug$y[data_aug$season == "Dry"]
yw <- data_aug$y[data_aug$season == "Wet"]
add_yd <- rep(0, times = length(yw) - length(yd))
yd <- c(yd, add_yd)
y <- array(c(yd, yw), dim = c(2, length(yd)))

nind <- c(sum(data_aug$y[data_aug$season == "Dry"]), sum(data_aug$y[data_aug$season == "Wet"]))
nz <- c(6517 - nind[1], 6517 - nind[2])

## Sort Veg Type Arrays
VegTypeM <- torts2020_dry %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM_dry <- model.matrix(~ -1 + VegType, data = VegTypeM)

VegTypeM <- torts2020_wet %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM_wet <- model.matrix(~ -1 + VegType, data = VegTypeM)

rowstoadd <- matrix(NA, nrow = 4, ncol = 8)
VegTypeM_wet <- rbind(VegTypeM_wet, rowstoadd)

VegTypeA <- array(c(VegTypeM_dry, VegTypeM_wet), dim = c(263, 8, 2))






















code <- nimbleCode({
  
  # Prior distributions
  # Intercept of lambda
  
  #for (t in 1:nSites){
  #  beta[t] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  #}
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  }
  
  for (year in 1:nyear){
    for (f in 1:nseason){
      
      # Intercept of log(sigma) for each season
      alpha0[f] ~ dunif(-10, 10) 
      
      # psi is a derived parameter under DA for stratified populations
      psi[year, f] <- sum(lambda[year, f, 1:nsites[f]]) / (nind[year, f] + nz[year, f])
      
      # 'Likelihood' (sort of...)
      for (i in 1:(nind[year, f] + nz[year, f])){               # i is index for individuals
        z[year, f, i] ~ dbern(psi[year, f])                    # Data augmentation variables
        d[year, f, i] ~ dunif(0, B)                      # distance uniformly distributed
        p[year, f, i] <- exp(-d[year, f, i] * d[year, f, i] / (2 * sigma[year, f, site[year, f, i]] * sigma[year, f, site[year, f, i]])) # Det. function
        mu[year, f, i] <- z[year, f, i] * p[year, f, i]              # 'straw man' for WinBUGS
        y[year, f, i] ~ dbern(mu[year, f, i])                  # basic Bernoulli random variable
      }
      
      for (s in 1:nsites[year, f]){                    # s is index for sites
        # next line not necessary, but allows to make predictions
        beta0[year, f, s] ~ dunif(0, 100)
        N[year, f, s] ~ dpois(lambda[year, f, s])              # Realized abundance at site s
        log(lambda[year, f, s]) <- beta0[year, s] #+ beta[1] * TransectA[f, 1, s] + beta[2] * TransectA[f, 2, s] + beta[3] * TransectA[f, 3, s] + beta[4] * TransectA[f, 4, s] +
        #beta[5] * TransectA[f, 5, s] + beta[6] * TransectA[f, 6, s] + beta[7] * TransectA[f, 7, s] + beta[8] * TransectA[f, 8, s] +
        #beta[9] * TransectA[f, 9, s] + beta[10] * TransectA[f, 10, s] + beta[11] * TransectA[f, 11, s] + beta[12] * TransectA[f, 12, s] +
        #beta[13] * TransectA[f, 13, s] + beta[14] * TransectA[f, 14, s] + beta[7] * TransectA[f, 15, s] + beta[8] * TransectA[f, 15, s] +
        #beta[16] * TransectA[f, 16, s] + beta[17] * TransectA[f, 17, s]
        
        # Linear model for detection
        log(sigma[year, f, s]) <- alpha0[f] + alpha[1] * VegTypeA[s, 1, f] + alpha[2] * VegTypeA[s, 2, f] + alpha[3] * VegTypeA[s, 3, f] +
          alpha[4] * VegTypeA[s, 4, f] + alpha[5] * VegTypeA[s, 5, f] + alpha[6] * VegTypeA[s, 6, f] +
          alpha[7] * VegTypeA[s, 7, f] + alpha[8] * VegTypeA[s, 8, f]
      }
    }
  }
  
  # Derived parameter: total population size across all site
  Ntotal2019_dry <- sum(z[1, 1, 1:6517])
  Ntotal2019_wet <- sum(z[1, 2, 1:6517])
  
})