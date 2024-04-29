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

torts$VegType[torts$VegType == "?"] <- "CHP"
torts$SectionCode <- as.factor(torts$SectionCode)
torts$VegType <- as.factor(torts$VegType)

# Remove entries left over from before 2018 and with distance >50m
torts <- torts %>%
  mutate(DateOfTransect = dmy(DateOfTransect)) %>%
  filter(year(DateOfTransect) > 2018) %>%
  filter(DistFromTransect <= 50) %>%
  select(VegType, DistFromTransect, DateOfTransect, SectionCode)

# plot data
hist(torts$DistFromTransect, main = "total tortoises", xlab = "Distance from transect (m)")

# select one year and one season
torts1 <- torts %>%
  filter(year(DateOfTransect) == 2019) %>%
  filter(month(DateOfTransect) < 6) %>%
  droplevels()

# plot data
hist(torts1$DistFromTransect, main = "total tortoises", xlab = "Distance from transect (m)")

# Get data and do data-augmentation
# Observed distances (meters)
x <- torts1$DistFromTransect

# Analysis of continuous data using data augmentation (DA)
nz <- 5000 # Augment observed data with nz = 200 zeroes
nind <- nrow(torts1)
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented
B <- 50 # Strip half-width. Larger than max distance
site <- c(torts1$SectionCode, rep(NA, nz))
nsites <- length(levels(torts1$SectionCode))
VegTypeM <- torts1 %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts1$VegType))

# nimble code
code <- nimbleCode( {
  
  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind+nz)
  
  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    x[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-x[i] * x[i] / (2 * sigma[site[i]] * sigma[site[i]])) # Det. function
    mu[i] <- z[i] * p[i]                 # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
  }
  
  for (v in 1:nVegTypes){
  beta[v] ~ dnorm(0, 10)   # Intercept of lambda-habitat regression
  alpha[v] ~ dnorm(0, 10) # Intercept of log(sigma) (half-normal scale)
  }
  
  beta0 ~ dunif(0, 10)
  alpha0 ~ dunif(0, 10)
  
  # Linear models for abundance and for detection
  for(s in 1:nsites){                    # s is index for sites
    
    # Model for abundance
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta[1] * VegTypeM[s, 1] + beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] +
                      beta[4] * VegTypeM[s, 4] + beta[5] * VegTypeM[s, 5] + beta[6] * VegTypeM[s, 6] +
                      beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8] + beta[9] * VegTypeM[s, 9] +
                      beta[10] * VegTypeM[s, 10] + beta[11] * VegTypeM[s, 11]
    
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    
    # model for detection 
    log(sigma[s]) <- alpha0 + alpha[1] * VegTypeM[s, 1] + alpha[2] * VegTypeM[s, 2] + alpha[3] * VegTypeM[s, 3] +
                     alpha[4] * VegTypeM[s, 4] + alpha[5] * VegTypeM[s, 5] + alpha[6] * VegTypeM[s, 6] +
                     alpha[7] * VegTypeM[s, 7] + alpha[8] * VegTypeM[s, 8] + alpha[9] * VegTypeM[s, 9] +
                     alpha[10] * VegTypeM[s, 10] + alpha[11] * VegTypeM[s, 11]
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind+nz)])
  #area <- nsites*1*2*B   # Unit length == 1, half-width = B
  #D <- Ntotal/area 
  }
)

# Bundle data set
data <- list(x = x,
             y = y,
             site = site)

consts <- list(nind = nind,
               nz = nz,
               B = B,
               nsites = nsites,
               VegTypeM = VegTypeM)

# Inits
zst <- y
inits <- list(beta = rep(1, 11),
              alpha = rep(1, 11),
              alpha0 = 0,
              beta0 = 0,
              z = zst,
              x = c(rep(NA, nind), runif(nz, 0, B)),
              N = rpois(nsites, exp(5)),
              site = c(rep(NA, nind),sample(1:nsites, nz, replace = TRUE)))

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
config$addMonitors(c("sigma", "psi", "Ntotal", "alpha", "beta", "alpha0", "beta0"))
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

saveRDS(run, "Outputs/2019W_VegType.rds")

mcmc_trace(run$samples, pars = c("sigma", "psi", "N"))

4735/3.45

1372*115

