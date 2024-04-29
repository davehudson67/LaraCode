# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(abind)
library(lubridate)
library(bayesplot)

# Read in distance sampling data and remove entries with no transect ID and no distance data
torts <- read.csv("distancetorts.v2.csv", header = TRUE) %>%
  drop_na(TransectID) %>%
  drop_na(DistFromTransect)

# deal with unknown VegTypes - recording errors?? Could add this as a latent variable if unknown.
# for now just changed both to CHP for simplicity.
torts$VegType[torts$VegType == "?"] <- "CHP"
torts$VegType[torts$VegType == "#N/A"] <- "CHP"
torts$SectionCode <- as.factor(torts$SectionCode)
torts$VegType <- as.factor(torts$VegType)

# Remove entries left over from before 2018 and with distance >50m
torts <- torts %>%
  mutate(DateOfTransect = dmy(torts$DateOfTransect)) %>%
  filter(year(DateOfTransect) > 2018) %>%
  filter(DistFromTransect <= 50) %>% 
  filter(month(DateOfTransect) != 6) %>%
  filter(month(DateOfTransect) != 7) %>%
  mutate(season = ifelse(month(DateOfTransect) < 6, "1", "2")) %>%
  mutate(year = year(DateOfTransect)) %>%
  unite("occasion", year, season, sep = "", remove = FALSE) %>%
  arrange(DateOfTransect) %>%
  mutate(occasion = as.numeric(as.factor(occasion))) %>%
  mutate(y = 1) %>%
  select(VegType, DistFromTransect, DateOfTransect, SectionCode, y, season, occasion) %>%
  droplevels()

# Create a function to augment data for each season
data_aug <- function(torts, year, n_aug) {
  # Filter the data for the specified year and season
  torts_year_w <- torts %>%
    filter(year(DateOfTransect) == year) %>%
    filter(season == "1") %>%
    droplevels()
  
  torts_year_d <- torts %>%
    filter(year(DateOfTransect) == year) %>%
    filter(season == "2") %>%
    droplevels()
  
  # Get the number of observations in the year
  n_ind_w <- nrow(torts_year_w)
  n_aug_w <- n_aug - n_ind_w
  
  n_ind_d <- nrow(torts_year_d)
  n_aug_d <- n_aug - n_ind_d
  
  # Create a vector of observed distances
  x_w <- c(torts_year_w$DistFromTransect, rep(NA, n_aug_w))
  x_d <- c(torts_year_d$DistFromTransect, rep(NA, n_aug_d))
  x <- cbind(x_w, x_d)
  
  # Create a vector of binary indicators for real or augmented data
  y_w <- c(rep(1, n_ind_w), rep(0, n_aug_w))
  y_d <- c(rep(1, n_ind_d), rep(0, n_aug_d))
  y <- cbind(y_w, y_d)
  
  # Create a vector of section codes
  site_w <- c(torts_year_w$SectionCode, rep(NA, n_aug_w))
  site_d <- c(torts_year_d$SectionCode, rep(NA, n_aug_d))
  site <- cbind(site_w, site_d)
  
  # Return the data in a list
  return(list(x = x, y = y, site = site))
}


# Create a function to cycle over years and create 3-dimensional arrays for x, y, and site
create_array <- function(torts, n_aug) {
  
  # Create a list to store the 3-dimensional arrays
  arrays <- list()
  
  # Cycle over years
  for (year in 2019:2022) {
    # Split the data for the specified year
    data <- data_aug(torts, year, n_aug)
    
    # Create 3-dimensional arrays for x, y, and site
    array_x <- abind(data$x, along = 2)
    array_y <- abind(data$y, along = 2)
    array_site <- abind(data$site, along = 2)
    
    # Store the 3-dimensional arrays in the list
    arrays[[year]] <- list(x = array_x, y = array_y, site = array_site)
  }
  
  # Return the list of 3-dimensional arrays
  return(arrays[2019:2022])
}

# Create the 3-dimensional arrays
arrays <- create_array(torts, 5000)

# Function to convert a list of matrices to a 3D array
convert_to_3D_array <- function(matrices_list) {
  # Bind matrices along the third dimension to create a 3D array
  return(abind(matrices_list, along = 3))
}

# Extract and convert the 'x', 'y', and 'site' variables for each year and season
x <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$x)))
y <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$y)))
site <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$site)))
available <- y
available[available == 0] <- NA

# Additional required data
B <- 50
nsites <- length(levels(torts$SectionCode))
VegTypeM <- torts %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts$VegType))
nseason <- length(levels(as.factor(torts$season)))
nind <- 5000 # augmented and real for each sampling occasion
nyears <- 4

#tidy
rm(arrays)

code <- nimbleCode( {
  
  beta0 ~ dunif(0, 5)
  for (v in 1:nVegTypes){
    beta[v] ~ dnorm(0, 1)
  }
  
  for (r in 1:nseasons){
    sigma[r] ~ dnorm(0, 1)
    available_prob[r] ~ dunif(0, 1)
  }
  
  # 'Likelihood' (sort of...)
  for(t in 1:nyears){
    
    for(k in 1:nseasons){
      
      # psi is a derived parameter under DA for stratified populations
      psi[k, t] <- sum(lambda[1:nsites, t]) / nind       # proportion of real individuals in the population
      #psiS[k, t] ~ dunif(0, 1)  
      
      for(i in 1:nind){                                  # i is index for individuals
        z[i, k, t] ~ dbern(psi[k, t])                    # Data augmentation variables (are you real?!)
        x[i, k, t] ~ dunif(0, B)                         # distance uniformly distributed
        p[i, k, t] <- exp(-x[i, k, t] * x[i, k, t] / (2 * sigma[k] * sigma[k])) # Det. function (if available)
        mu[i, k, t] <- (z[i, k, t] * (available[i, k, t])) * p[i, k, t]  # expected prob of observing an individual
        y[i, k, t] ~ dbern(mu[i, k, t])                  # actual observations - basic Bernoulli random variable
        #site[i, k, t] ~ dcat(site.probs[1:nsites, t])   # Population distribution among sites
        available[i, k, t] ~ dbern(available_prob[k])    # is individual available for detection
      } # individual
    } # season
    
    for(s in 1:nsites){                          # s is index for sites
      N[s, t] ~ dpois(lambda[s, t])              # Realized abundance at site s
      log(lambda[s, t]) <- beta0 + beta[1] * VegTypeM[s, 1] + beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] +
        beta[4] * VegTypeM[s, 4] + beta[5] * VegTypeM[s, 5] + beta[6] * VegTypeM[s, 6] +
        beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8] + beta[9] * VegTypeM[s, 9] +
        beta[10] * VegTypeM[s, 10] + beta[11] * VegTypeM[s, 11]
      
      #site.probs[s, t] <- lambda[s, t] / sum(lambda[1:nsites, t])
    } # sites
  } # years
  
  # Derived parameter: total population size across all sites
  Ntotal2019_W <- sum(z[1:nind, 1, 1])
  Ntotal2019_D <- sum(z[1:nind, 2, 1])
  Ntotal2020_W <- sum(z[1:nind, 1, 2])
  Ntotal2020_D <- sum(z[1:nind, 2, 2])
  Ntotal2021_W <- sum(z[1:nind, 1, 3])
  Ntotal2021_D <- sum(z[1:nind, 2, 3])
  Ntotal2022_W <- sum(z[1:nind, 1, 4])
  Ntotal2022_D <- sum(z[1:nind, 2, 4])
}
)

# Bundle data set
data <- list(x = x,
             y = y,
             available = available, 
             site = site)

consts <- list(nind = nind,
               B = B,
               nsites = nsites,
               nseasons = nseason,
               VegTypeM = VegTypeM,
               nyears = nyears)

# Inits
zst <- y

x_init <- x
xNA <- is.na(x_init)
x_init[xNA] <- runif(sum(xNA), 0.1, B)
x_init[!xNA] <- NA

site_init <- site
siteNA <- is.na(site)
site_init[siteNA] <- sample(1:nsites, sum(siteNA), replace = TRUE)
site_init[!siteNA] <- NA

available_init <- available 
availableNA <- is.na(available_init)
available_init[availableNA] <- rbinom(sum(availableNA), 1, 0.5)
available_init[!availableNA] <- NA

inits <- list(beta = runif(11, 0, 0.2),
              beta0 = 0.5,
              available = available_init,
              available_prob = c(1, 1), 
              psiS = matrix(runif(8, 0, 1), nrow = 2),
              z = zst,
              x = x_init, 
              N = matrix(rpois(nsites * 4, exp(5)), ncol = 4),
              #site = site_init,
              sigma = rnorm(2, 0, 1))

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
#model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
config$removeSamplers("beta", "beta0", "available_prob")
config$addSampler(target = c("beta", "beta0"), type = 'AF_slice')
config$addSampler(target = c("available_prob"), type = 'AF_slice')
#config$addSampler(target = c("psiS"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("sigma", "psi", "psiS", "Ntotal2019_W", "Ntotal2019_D",
                     "Ntotal2020_W", "Ntotal2020_D",
                     "Ntotal2021_W", "Ntotal2021_D",
                     "Ntotal2022_W", "Ntotal2022_D", "beta", "beta0",
                     "available_prob"))
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 50000, 
                           nburnin = 12000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, "Outputs/VegSeasonYear_TE_50K_AFslice.rds")
run <- readRDS("Outputs/VegSeasonYear_TE_50K_AFslice.rds")
run$summary
mcmc_trace(run$samples)
library(MCMCvis)

