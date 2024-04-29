# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(abind)
library(lubridate)
library(bayesplot)

## load formatted data
load("tortsReady.RData")

## try one season
torts19W <- torts %>%
  filter(occasion == "2019_Wet")

x <- c(torts19W$DistFromTransect, rep(NA, times = 5000 - nrow(torts19W)))
y <- c(rep(1, nrow(torts19W)), rep(0, times = 5000 - nrow(torts19W)))

# nimble code
code <- nimbleCode( {
  
  # Priors
  sigma ~ dunif(0, 100)  # Half-normal scale
  psi ~ dunif(0, 1)       # DA parameter 
  
  # Likelihood
  for(i in 1:(nind)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    
    # Observation model
    #logp[i] <- -((x[i] / sigma)) # negative exponential fct
    logp[i] <- -((x[i] * x[i]) / (2 * sigma * sigma)) # Half-normal detection fct. - BOOK
    p[i] <- exp(logp[i])
    
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
  }
  # Derived quantities
  N <- sum(z[1:5000])  # Population size
}
)

# Bundle data set
data <- list(x = x,
             y = y)

consts <- list(nind = nind,
               B = B)
# Inits
z_init <- y
z_init[z_init == 0] <- rbinom(5000 - nrow(torts19W), 1, 0.5)

x_init <- x
x_init[is.na(x_init)] <- runif(5000 - nrow(torts19W), 0.1, B)

inits <- list(sigma = runif(1, 0, 100),
              psi = runif(1, 0, 1),
              x = x_init,
              z = z_init)

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
#model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("sigma", "psi", "N"))
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
plot(run$samples)

##########################################################################################################
##########################################################################################################

## Model 2
## try one season
torts19W <- torts %>%
  filter(occasion == "2019_Wet") %>%
  droplevels()

nsites <- length(levels(torts19W$SectionCode))
VegTypeM <- torts19W %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts19W$VegType))
nind <- 5000 # augmented and real for each sampling occasion

x <- c(torts19W$DistFromTransect, rep(NA, times = 5000 - nrow(torts19W)))
y <- c(rep(1, nrow(torts19W)), rep(0, times = 5000 - nrow(torts19W)))

## Still going for one season, but now site-specific vegetation effects
site <- c(torts19W$SectionCode, rep(NA, times = 5000 - nrow(torts19W)))
nsites <- length(levels(torts19W$SectionCode)) 

# nimble code
code <- nimbleCode({
  # Priors
  sigma ~ dunif(0, 100)  # Half-normal scale for detection distance
  alpha0 ~ dunif(0, 10)  # Intercept for log-scale site-level abundance
  
  for (v in 1:nVegTypes) {
    beta[v] ~ dnorm(0, 1)  # Effect of vegetation type on abundance
  }
  
  # Transect-level abundance model
  for (t in 1:nTransects) {
    lambdaTransect[t] ~ dpois(transectLambdaRate) # Total abundance for each transect
  }
  
  # Vegetation effect on site-level abundance within each transect
  for (s in 1:nsites) {
    log(lambda[s]) <- alpha0 + sum(beta[1:nVegTypes] * VegTypeM[s, 1:nVegTypes])
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
  }
  
  # Likelihood for individual observations
  for (i in 1:nind) {
    transect[i] ~ dcat(transect.probs[1:nTransects]) # Assigning transect to individual
    site[i] ~ dcat(site.probs[1:nsites]) # Assigning site within transect to individual
    
    # Process model
    z[i] ~ dbern(lambdaTransect[transect[i]] / nind) # DA variables based on transect abundance
    
    x[i] ~ dunif(0, B) # Distribution of distances
    
    # Observation model influenced by season and vegetation
    logp[i] <- -((x[i] * x[i]) / (2 * sigma * sigma)) # Half-normal detection function
    p[i] <- exp(logp[i])
    
    mu[i] <- z[i] * p[i] * seasonEffect[season[i]] * vegetationEffect[site[i]]
    y[i] ~ dbern(mu[i]) # Bernoulli measurement error process
  }
  
  # Derived quantities
  Ntot <- sum(z[1:nind]) # Total estimated population size
})


# Bundle data set
data <- list(x = x,
             y = y,
             site = site)

consts <- list(nind = nind,
               B = B, 
               nsites = nsites,
               VegTypeM = VegTypeM,
               nVegTypes = nVegTypes)
# Inits
z_init <- y
z_init[z_init == 0] <- rbinom(5000 - nrow(torts19W), 1, 0.5)

x_init <- x
x_init[is.na(x_init)] <- runif(5000 - nrow(torts19W), 0.1, B)

site_init <- site
site_init[is.na(site_init)] <- sample(1:nsites, sum(is.na(site)), replace = TRUE)

inits <- list(sigma = runif(1, 0, 100),
              alpha0 = 0,
              beta = rnorm(nVegTypes, 0, 1),
              x = x_init,
              z = z_init,
              site = site_init,
              N = rep(10, nsites))

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
#model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("sigma", "psi", "alpha0", "beta", "N"))
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
mcmc_trace(run$samples, pars = c("alpha0", "beta[1]", "sigma"))

samplesS <- as.matrix(run$summary$all.chains)
sum(samplesS[1:259, 1])
##########################################################################################################
##########################################################################################################

## Model 3
## try one season
torts19W <- torts %>%
  filter(occasion == "2019_Wet") %>%
  droplevels()

nsites <- length(levels(torts19W$SectionCode))
VegTypeM <- torts19W %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts19W$VegType))
nind <- 5000 # augmented and real for each sampling occasion

x <- c(torts19W$DistFromTransect, rep(NA, times = 5000 - nrow(torts19W)))
y <- c(rep(1, nrow(torts19W)), rep(0, times = 5000 - nrow(torts19W)))

## Still going for one season, but now site-specific vegetation effects
site <- c(torts19W$SectionCode, rep(NA, times = 5000 - nrow(torts19W)))
nsites <- length(levels(torts19W$SectionCode)) 

# nimble code
code <- nimbleCode( {
  
  # Priors
  alpha0 ~ dunif(0, 10)
  beta0 ~ dunif(0, 10)
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dnorm(0, 1)
    beta[v] ~ dnorm(0, 1)
    
  }
  
  psi <- sum(lambda[1:nsites] / nind)
  
  # Likelihood
  for(i in 1:(nind)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    
    # Observation model
    #logp[i] <- -((x[i] / sigma)) # negative exponential fct
    logp[i] <- -((x[i] * x[i]) / (2 * sigma[site[i]] * sigma[site[i]])) # Half-normal detection fct. - BOOK
    p[i] <- exp(logp[i])
    
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
    site[i] ~ dcat(site.probs[1:nsites])
  }
  
  for(s in 1:nsites){
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta[1] * VegTypeM[s, 1] + beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] +
      beta[4] * VegTypeM[s, 4] + beta[5] * VegTypeM[s, 5] + beta[6] * VegTypeM[s, 6] +
      beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8]
    
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    log(sigma[s]) <- alpha0 + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
      alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
      alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
  }

}
)

# Bundle data set
data <- list(x = x,
             y = y,
             site = site)

consts <- list(nind = nind,
               B = B, 
               nsites = nsites,
               VegTypeM = VegTypeM,
               nVegTypes = nVegTypes)
# Inits
z_init <- y
z_init[z_init == 0] <- rbinom(5000 - nrow(torts19W), 1, 0.5)

x_init <- x
x_init[is.na(x_init)] <- runif(5000 - nrow(torts19W), 0.1, B)

site_init <- site
site_init[is.na(site_init)] <- sample(1:nsites, sum(is.na(site)), replace = TRUE)

inits <- list(alpha0 = 1,
              beta0 = 0,
              alpha = rnorm(nVegTypes, 0, 1),
              beta = rnorm(nVegTypes, 0, 1),
              x = x_init,
              z = z_init,
              site = site_init,
              N = rep(10, nsites))

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
#model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("alpha0", "alpha", "beta0", "beta", "N", "sigma"))
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
mcmc_trace(run$samples, pars = c("alpha0", "beta0", "beta[1]", "alpha[1]"))

samplesS <- as.matrix(run$summary$all.chains)
sum(samplesS[1:259, 1]) #3540
##########################################################################################################
##########################################################################################################

##########################################################################################################
##########################################################################################################

## Model 4
## try one year... so 2 seasons
torts19 <- torts %>%
  filter(occasion == "2019_Wet" | occasion == "2019_Dry") %>%
  droplevels()

nsites <- length(levels(torts19W$SectionCode))
VegTypeM <- torts19 %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts19$VegType))
nind <- 5000 # augmented and real for each sampling occasion

data_aug <- function(torts, year, n_aug) {

  # Filter the data for the specified year and season
  torts_year_w <- torts %>%
    filter(Year == year, season == "Wet") %>%
    droplevels()
  
  torts_year_d <- torts %>%
    filter(Year == year, season == "Dry") %>%
    droplevels()
  
  # Calculate the number of augmented observations needed
  n_ind_w <- nrow(torts_year_w)
  n_aug_w <- n_aug - n_ind_w
  
  n_ind_d <- nrow(torts_year_d)
  n_aug_d <- n_aug - n_ind_d
  
  # Create vectors of observed and augmented distances
  x_w <- c(torts_year_w$DistFromTransect, rep(NA, n_aug_w))
  x_d <- c(torts_year_d$DistFromTransect, rep(NA, n_aug_d))
  
  # Create vectors of binary indicators for real or augmented data
  y_w <- c(rep(1, n_ind_w), rep(0, n_aug_w))
  y_d <- c(rep(1, n_ind_d), rep(0, n_aug_d))
  
  # Create vectors of section codes
  site_w <- c(torts_year_w$SectionCode, rep(NA, n_aug_w))
  site_d <- c(torts_year_d$SectionCode, rep(NA, n_aug_d))
  
  # Return the data in a list
  return(list(x = cbind(x_w, x_d), y = cbind(y_w, y_d), site = cbind(site_w, site_d)))
}

create_array <- function(torts, year, n_aug) {
  # Split the data for the specified year
  data <- data_aug(torts, year, n_aug)
  
  # Return the data as a list of arrays
  return(list(x = data$x, y = data$y, site = data$site))
}

arrays <- create_array(torts19, 2019, 5000)

# Extracting the 2D arrays
x <- arrays$x
y <- arrays$y
site <- arrays$site
nsites <- length(levels(torts19$SectionCode)) 

# nimble code
code <- nimbleCode({
  
  # Priors for season-specific parameters
  for(season in 1:nseason) {
    alpha0[season] ~ dunif(0, 10)
    
    for (v in 1:nVegTypes) {
      alpha[v, season] ~ dnorm(0, 1)
      beta[v] ~ dnorm(0, 1)
    }
  }
  
  beta0 ~ dunif(0, 10)
  
  # Calculate total detection probability
  psi <- sum(lambda[1:nsites] / nind)
  
  # Likelihood for individual observations
  for(i in 1:nind) {
    for(season in 1:nseason) {
      # Process model
      z[i, season] ~ dbern(psi)
      x[i, season] ~ dunif(0, B)
      site[i, season] ~ dcat(site.probs[1:nsites])
      
      # Observation model
      logp[i, season] <- -((x[i, season] * x[i, season]) / (2 * sigma[site[i, season], season] * sigma[site[i, season], season]))
      p[i, season] <- exp(logp[i, season])
      mu[i, season] <- z[i, season] * p[i, season]
      y[i, season] ~ dbern(mu[i, season])
    }
  }
  
  # Site-specific models
  for(s in 1:nsites) {
    N[s] ~ dpois(lambda[s]) # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta[1] * VegTypeM[s, 1] +  beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] + beta[4] * VegTypeM[s, 4] +
                              beta[5] * VegTypeM[s, 5] +  beta[6] * VegTypeM[s, 6] + beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8]
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    for(season in 1:nseason) {   
      log(sigma[s, season]) <- alpha0[season] + alpha[1] * VegTypeM[s, 1] +  alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] + alpha[4] * VegTypeM[s, 4] +
                                                alpha[5] * VegTypeM[s, 5] +  alpha[6] * VegTypeM[s, 6] + alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
    }
  }
})


# Bundle data set
data <- list(x = x,
             y = y,
             site = site)

consts <- list(nind = nind,
               B = B, 
               nsites = nsites,
               VegTypeM = VegTypeM,
               nVegTypes = nVegTypes,
               nseason = nseason)
# Inits
z_init <- y
z_init[z_init == 0] <- rbinom(length(which(z_init == 0)), 1, 0.5)

x_init <- x
x_init[is.na(x_init)] <- runif(length(which(y == 0)), 0.1, B)

site_init <- site
site_init[is.na(site_init)] <- sample(1:nsites, sum(is.na(site)), replace = TRUE)

inits <- list(alpha0 = rep(1, 2),
              beta0 = rep(0, 2),
              alpha = matrix(rnorm(nVegTypes * nseason, 0, 1), nrow = 2),
              beta = rnorm(nVegTypes, 0, 1),
              x = x_init,
              z = z_init,
              site = site_init,
              N = rep(10, nsites))

# Build model
model <- nimbleModel(code, constants = consts, data = data)
#model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("alpha0", "alpha", "beta0", "beta", "N", "sigma"))
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
mcmc_trace(run$samples, pars = c("alpha0", "beta0", "beta[1]", "alpha[1]"))

samplesS <- as.matrix(run$summary$all.chains)
sum(samplesS[1:259, 1]) #3540






code <- nimbleCode( {
  
  alpha0 ~ dunif(0, 10)
  beta0 ~ dunif(0, 10)
  
  for (v in 1:nVegTypes){
    beta[v] ~ dnorm(0, 10)
    alpha[v] ~ dnorm(0, 10)
  }
  
  # 'Likelihood'
  for(t in 1:nyears){
    for(k in 1:nseasons){
      
      #sigma[k, t] ~ dunif(0, 100)
      # psi is a derived parameter under DA for stratified populations
      psi[k, t] <- sum(lambda[1:nsites]) / nind
      psiS[k, t] ~ dunif(0, 1)  
      
      for(i in 1:nind){                       # i is index for individuals
        z[i, k, t] ~ dbern(psi[k, t])                    # Data augmentation variables (are you real?!)
        x[i, k, t] ~ dunif(0, B)                         # distance uniformly distributed
        p[i, k, t] <- exp(-x[i, k, t] * x[i, k, t] / (2 * sigma[site[i, k, t]] * sigma[site[i, k, t]])) # Det. function
        mu[i, k, t] <- z[i, k, t] * p[i, k, t]           # 'straw man' for WinBUGS
        y[i, k, t] ~ dbern(mu[i, k, t])                  # basic Bernoulli random variable
        site[i, k, t] ~ dcat(site.probs[1:nsites])       # Population distribution among sites
      }
    }
  }
  
  ## Models for abundance and for detection
  # Abundance dependent on VegType
  for(s in 1:nsites){                    # s is index for sites
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta[1] * VegTypeM[s, 1] + beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] +
      beta[4] * VegTypeM[s, 4] + beta[5] * VegTypeM[s, 5] + beta[6] * VegTypeM[s, 6] +
      beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8]
    
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    log(sigma[s]) <- alpha0 + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
      alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
      alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8]
  }
  
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

inits <- list(beta = rep(0, 8),
              beta0 = 1,
              alpha = rep(0, 8),
              alpha0 = 1,
              psiS = matrix(runif(8, 0, 1), nrow = 2),
              z = zst,
              x = x_init, 
              N = rpois(nsites, exp(5)),
              site = site_init)
#sigma = matrix(runif(8, 0, 100), nrow = 2))

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("sigma", "psi", "psiS", "Ntotal2019_W", "Ntotal2019_D",
                     "Ntotal2020_W", "Ntotal2020_D",
                     "Ntotal2021_W", "Ntotal2021_D",
                     "Ntotal2022_W", "Ntotal2022_D", "beta", "beta0"))
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

