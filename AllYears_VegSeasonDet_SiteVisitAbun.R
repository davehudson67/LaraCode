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
#torts <- filter(torts, Year >= 2021)

head(torts)
levels(torts$occasion)
torts <- filter(torts, Year != 2018, .preserve = TRUE)

data_aug <- torts %>%
  group_by(Year, season) %>%
  do({
    group_data <- .
    # Determine sites visited in the current year and season
    sites_visited <- unique(group_data$SectionCode)
    
    # Number of augmentations per site visited
    nz <- length(sites_visited) * 20
    
    # Prepare augmented site data
    aug_site <- rep(sites_visited, each = 20)
    aug_site_codes <- factor(aug_site, levels = unique(group_data$SectionCode))
    aug_dist <- rep(NA, nz)
    
    augmented_year <- rep(unique(group_data$Year), times = nz)
    augmented_season <- rep(unique(group_data$season), times = nz)
    
    # Combine original and augmented data
    data.frame(
      Year = c(group_data$Year, augmented_year),
      season = c(group_data$season, augmented_season),
      SectionCode = factor(c(group_data$SectionCode, aug_site_codes), levels = unique(torts$SectionCode)),
      DistFromTransect = c(group_data$DistFromTransect, aug_dist),
      y = c(rep(1, nrow(group_data)), rep(0, nz))  # Distinguish augmented data
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

## Site visited indicator 
# Create a complete dataset with all combinations of site, season, and year
full_combinations <- expand.grid(
  SectionCode = unique(torts$SectionCode),
  season = unique(torts$season),
  Year = unique(torts$Year)
)

# Merge with the original dataset to identify visited combinations
visited_combinations <- torts %>%
  select(SectionCode, season, Year) %>%
  distinct() %>%
  mutate(visited = 1)  # Mark these as visited

# Join to ensure all combinations are present, filling in 'visited' as 0 where not matched
all_combinations <- merge(full_combinations, visited_combinations, by = c("SectionCode", "season", "Year"), all.x = TRUE) %>%
  replace_na(list(visited = 0))

# Convert the data frame to a 3D array
all_combinations <- all_combinations %>%
  arrange(SectionCode, season, Year)

# Convert 'visited' column to a matrix, then to an array
visited_matrix <- matrix(all_combinations$visited, nrow = length(unique(all_combinations$SectionCode)), byrow = TRUE)

# Convert matrix to a 3D array with dimensions for site, season, year
dim(visited_matrix) <- c(
  length(unique(all_combinations$SectionCode)), # Site dimension
  length(unique(all_combinations$season)),      # Season dimension
  length(unique(all_combinations$Year))         # Year dimension
)

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
      p[i] <- exp(-d[i] * d[i] / (2 * sigma[site[i], season[i], year[i]] * sigma[site[i], season[i], year[i]])) # Det. function
      mu[i] <- z[i] * p[i]                  # 'straw man' for WinBUGS
      y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    }
  
  # Linear models for abundance and for detection
  for (year in 1:nyears){
    for(s in 1:nsites){                    # s is index for sites
      # Model for abundance
      # next line not necessary, but allows to make predictions
      N[s, year, season] ~ dpois(lambda[s, year])              # Realized abundance at site s
      log(lambda[s, year]) <- beta0 * sitesvisited[s, year]
      
      # Linear model for detection
      for(seas in 1:nseason){
        log(sigma[s, seas, year]) <- alpha0[seas] + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
          alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
          alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
      }
    }
  }
  # Derived parameter: total population size across all sites
  Ntotal2019 <- sum(z[1:12766])
  Ntotal2020 <- sum(z[12767:24118])
  Ntotal2021 <- sum(z[24119:35848])
  Ntotal2022 <- sum(z[35849:39953])
})

#data_aug$Year <- as.factor(data_aug$Year)
#summary(data_aug$Year)
Ndata <- visited_matrix
Ndata[Ndata == 1] <- NA

# Bundle and summarize data set
data <- list(y=y, d=d, N=Ndata)

consts <- list(nsites=nsites, VegTypeM = VegTypeM, nVegTypes = nVegTypes, B=B, nind=nind, 
               nz=nz, site=site, season= season, nseason=nseason, sitesvisited=sitesvisited)

# Inits
zst <- y
inits <- function(){list(beta0=0, 
                         #beta1=0, 
                         alpha0 = rep(1, 2), 
                         alpha = rep(0, nVegTypes),
                         sigma = matrix(0, nsites, nseason),
                         z=zst)}

# Parameters to save
params <- c("alpha0", "alpha", "beta0", "psi", "N", "Ntotal")

# MCMC settings
ni <- 80000   ;   nb <- 12000   ;   nt <- 2   ;   nc <- 2

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
                   data = data,
                   constants = consts, 
                   inits = inits, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   nchains = nc,
                   samplesAsCodaMCMC = TRUE)

saveRDS(out1, "C2019_VSD_SVA_80Kout.rds")
mcmcplot(out1)