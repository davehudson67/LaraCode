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

saveRDS(torts,  "TortsReady.rds")

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

## Extract data
nYears <- length(unique(data_aug$Year))
nSeasons <- length(unique(data_aug$season))

# Assuming 'individual' is implicitly defined by the number of rows per Year and Season
data_aug <- data_aug %>% 
  group_by(Year, season) %>% 
  mutate(individual = row_number())
nIndividuals <- max(data_aug$individual)

# Initialize the 3D arrays
distance <- array(NA, dim = c(nYears, nSeasons, nIndividuals))
y <- array(NA, dim = c(nYears, nSeasons, nIndividuals))
site <- array(NA, dim = c(nYears, nSeasons, nIndividuals))

# Populate the 3D arrays
for (i in 1:nrow(data_aug)) {
  yearIndex <- as.integer(factor(data_aug$Year[i], levels = unique(data_aug$Year)))
  seasonIndex <- as.integer(factor(data_aug$season[i], levels = unique(data_aug$season)))
  individualIndex <- data_aug$individual[i]
  
  distance[yearIndex, seasonIndex, individualIndex] <- data_aug$DistFromTransect[i]
  y[yearIndex, seasonIndex, individualIndex] <- data_aug$y[i]
  site[yearIndex, seasonIndex, individualIndex] <- data_aug$SectionCode[i]
}

# Calculate the number of observations per year/season
obs_counts <- data_aug %>%
  group_by(Year, season) %>%
  summarise(n_individuals = n(), .groups = 'drop')

# Create a matrix to hold these counts
years <- sort(unique(data_aug$Year))
seasons <- sort(unique(data_aug$season))

# Initialize the matrix with NA or 0
individuals_matrix <- matrix(NA, nrow = nYears, ncol = nSeasons,
                             dimnames = list(years, seasons))

# Populate the matrix with the counts
for (row in 1:nrow(obs_counts)) {
  year <- obs_counts$Year[row]
  season <- obs_counts$season[row]
  count <- obs_counts$n_individuals[row]
  
  individuals_matrix[as.character(year), as.character(season)] <- count
} 

## Create transect site matrix
TransectM <- torts %>% 
  select(TransectCode, SectionCode) %>%          # Select only the relevant columns
  distinct() %>%                     # Remove duplicates if any
  mutate(value = 1) %>%              # Add a dummy variable column filled with 1s
  spread(key = TransectCode, value = value, fill = 0)

TransectM <- as.matrix(TransectM)
TransectM <- TransectM[,-1]

################################################################################

code <- nimbleCode({
  
  # Prior distributions
  for (v in 1:nVegTypes){
    alpha[v] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  }
  
  for (seas in 1:nseason){
    alpha0[seas] ~ dunif(-10, 10)  # Intercept of log(sigma) for each season
    beta0[seas] ~ dunif(0, 100)   # Intercept of lambda
  }
  
  for (year in 1:nyears){
    for (f in 1:nseason){
    
      # psi is a derived parameter under DA for stratified populations
      psi[year, f] <- sum(lambda[year, f, 1:nsites]) / (nind[year, f])
      
      # 'Likelihood' (sort of...)
      for (i in 1:(nind[year, f])){               # i is index for individuals
        z[year, f, i] ~ dbern(psi[year, f])                    # Data augmentation variables
        d[year, f, i] ~ dunif(0, B)                      # distance uniformly distributed
        p[year, f, i] <- exp(-d[year, f, i] * d[year, f, i] / (2 * sigma[f, site[year, f, i]] * sigma[f, site[year, f, i]])) # Det. function
        mu[year, f, i] <- z[year, f, i] * p[year, f, i]              # 'straw man' for WinBUGS
        y[year, f, i] ~ dbern(mu[year, f, i])                  # basic Bernoulli random variable
      }
      
      for (s in 1:nsites){                    # s is index for sites
        # next line not necessary, but allows to make predictions
        N[year, f, s] ~ dpois(lambda[year, f, s])              # Realized abundance at site s
        log(lambda[year, f, s]) <- beta0[f] + beta1[1] * transect[s, 1] + beta1[2] * transect[s, 2] + 
          beta1[3] * transect[s, 3] + beta1[4] * transect[s, 4] + beta1[5] * transect[s, 5] + beta1[6] * transect[s, 6] + 
          beta1[7] * transect[s, 7] + beta1[8] * transect[s, 8] + beta1[9] * transect[s, 9] + beta1[10] * transect[s, 10] + 
          beta1[11] * transect[s, 11] + beta1[12] * transect[s, 12] + beta1[13] * transect[s, 13] + beta1[14] * transect[s, 14] + 
          beta1[15] * transect[s, 15] + beta1[16] * transect[s, 16] + beta1[17] * transect[s, 17]
      }
    }
  }
  
  for (dseason in 1:nseason){
    for (dsite in 1:nsites){
      # Linear model for detection
      log(sigma[dseason, dsite]) <- alpha0[dseason] + alpha[1] * VegTypeA[dsite, 1] + alpha[2] * VegTypeA[dsite, 2] + alpha[3] * VegTypeA[dsite, 3] +
        alpha[4] * VegTypeA[dsite, 4] + alpha[5] * VegTypeA[dsite, 5] + alpha[6] * VegTypeA[dsite, 6] +
        alpha[7] * VegTypeA[dsite, 7] + alpha[8] * VegTypeA[dsite, 8]
    }
  }
  
})

# Bundle and summarize data set
data <- list(y=y, d=distance)

consts <- list(nsites=nrow(TransectM), VegTypeA = VegTypeM, nVegTypes = nVegTypes, B=B, nind=individuals_matrix, 
               nseason=nseason, nyears = nYears, transect = TransectM, site = site)
# Inits
zst <- y
d_init <- d
d_init[is.na(d_init)] <- runif(length(d_init[is.na(d_init)]), 0, B)
N_init <- array(ceiling(runif(nseason * nsitesMax, 1, 30)), dim = c(nseason, nsitesMax))
sigma_init <- array(1, dim = c(nseason, nsitesMax))
lambda_init <- array(exp(1), dim = c(nseason, nsitesMax))

inits <- list(beta0 = 1, 
              alpha0 = rep(1, 2), 
              alpha = rep(0, nVegTypes),
              N = N_init,
              z = zst,
              d = d_init,
              sigma = sigma_init)
#              lamda = lambda_init)

# Parameters to save
params <- c("alpha0", "alpha", "psi", "Ntotal")

# Build model
model <- nimbleModel(code, constants = consts, data = data)#, inits = inits)
model$initializeInfo()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
config$removeSamplers("alpha")
config$addSampler(target = c("alpha"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("psi", "alpha", "alpha0", "beta0", "Ntotal_dry", "Ntotal_wet"))
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 80000, 
                           nburnin = 12000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, "C2019_VSD_SVA_80Kout.rds")

parms <- c("psi", "z", "alpha", "alpha0")
mcmcplot(run)#, parms = parms)

# Create the summary table
summary_table <- torts %>%
  group_by(Year, season, VegType) %>%
  summarise(Observations = n(), .groups = 'drop') %>%
  arrange(Year, season, VegType)

summary_table <- torts %>%
  unite("YearSeason", Year, season, sep = "") %>% # Combine Year and Season into a single column
  group_by(YearSeason, VegType) %>%
  summarise(Observations = n(), .groups = 'drop') %>%
  pivot_wider(names_from = VegType, values_from = Observations, values_fill = list(Observations = 0)) # Reshape the data


