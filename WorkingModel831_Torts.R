# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
#library(mcmcplots)
library(lubridate)
library(bayesplot)

# Read in distance sampling data and remove entries with no transect ID and no distance data
torts <- read.csv("distancetorts.v2.csv", header = TRUE) %>%
  drop_na(TransectID) %>%
  drop_na(DistFromTransect)

# Remove entries left over from before 2018 and with distance >50m
torts <- torts %>%
  mutate(DateOfTransect = dmy(torts$DateOfTransect)) %>%
  filter(year(DateOfTransect) > 2018) %>%
  filter(DistFromTransect <= 50)

# plot data
hist(torts$DistFromTransect, main = "total tortoises", xlab = "Distance from transect (m)")

# select one year and one season
torts1 <- torts %>%
  filter(year(DateOfTransect) == 2019) %>%
  filter(month(DateOfTransect) < 6)

# plot data
hist(torts1$DistFromTransect, main = "total tortoises", xlab = "Distance from transect (m)")

# Get data and do data-augmentation
# Observed distances (meters)
x <- torts1$DistFromTransect

B <- 50 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
nz <- 5000 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented


# nimble code
code <- nimbleCode( {
  
  # Priors
  sigma ~ dunif(0, 1000)  # Half-normal scale
  psi ~ dunif(0, 1)       # DA parameter 
  
  # Likelihood
  for(i in 1:(nind+nz)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    
    # Observation model
    logp[i] <- -((x[i]*x[i])/(2*sigma*sigma)) # Half-normal detection fct. - BOOK
    #logp[i] <- -((x[i] * sigma)) # negative exponential fct  
    #logp[i] <- -(x[i] / sigma)
    p[i] <- exp(logp[i])
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
  }
  # Derived quantities
  N <- sum(z[1:(nind + nz)])  # Population size
  #D <- N / 2.515             # Density, with A = 60 km^2 when B = 500
}
)

# Bundle data set
data <- list(x = x,
             y = y)

consts <- list(nind = nind,
               nz = nz,
               B = B)

# Inits
zst <- y
inits <- list(psi = runif(1),
              z = zst,
              sigma = runif(1, 1, 100))

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
model$initializeInfo()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a", "b", "c1"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("sigma", "psi", "N", "x", "z"))
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

saveRDS(run, "2019W_Basic")

mcmc_trace(run$samples, pars = c("sigma", "psi", "N"))

4735/3.45

1372*115

