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
  select(VegType, DistFromTransect, DateOfTransect, SectionCode, y, season, occasion)

# plot data
ggplot(torts, aes(x = DistFromTransect, colour = season, fill = season)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ year(DateOfTransect))

# select one year
torts1 <- torts %>%
  filter(year(DateOfTransect) == 2019) %>%
  droplevels()

# Get data and do data-augmentation
nz <- 5000 # Augment observed data with nz = 200 zeroes
nind <- nrow(torts1)
nindW <- length(torts1$DistFromTransect[torts1$season == 1])
nindD <- length(torts1$DistFromTransect[torts1$season == 2])
nind <- nindD + nindW

# Observed distances (meters)
x <- c(torts1$DistFromTransect, rep(NA, nz))
xW <- c(torts1$DistFromTransect[torts1$season == 1], rep(NA, nz - nindW))
xD <- c(torts1$DistFromTransect[torts1$season == 2], rep(NA, nz - nindD))
x <- cbind(xW, xD)

# Real or Augmented
y <- c(rep(1, nind), rep(0, nz))
yW <- c(rep(1, nindW ), rep(0, nz - nindW)) # Augmented inds. have y=0 by definition
yD <- c(rep(1, nindD ), rep(0, nz - nindD)) # Augmented inds. have y=0 by definition
y <- cbind(yW, yD)

# Analysis of continuous data using data augmentation (DA)
B <- 50 # Strip half-width. Larger than max distance
site <- c(torts1$SectionCode, rep(NA, nz))
siteW <- c(torts1$SectionCode[torts1$season == 1], rep(NA, nz - nindW))
siteD <- c(torts1$SectionCode[torts1$season == 2], rep(NA, nz - nindD))
site <- cbind(siteW, siteD)

nsites <- length(levels(torts1$SectionCode))

VegTypeM <- torts1 %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts1$VegType))

nseason <- 2
seasonW <- c(as.numeric(torts1$season)[torts1$season == 1], rep(NA, nz - nindW))
seasonD <- c(as.numeric(torts1$season)[torts1$season == 2], rep(NA, nz - nindD))
season <- cbind(seasonW, seasonD)

code <- nimbleCode( {
  
  beta0 ~ dunif(0, 10)
  for (v in 1:nVegTypes){
    beta[v] ~ dnorm(0, 10)   # Intercept of lambda-habitat regression
  }
  
  # 'Likelihood' (sort of...)
  for(k in 1:nseasons){
  
    sigma[k] ~ dunif(0, 100)
      
    # psi is a derived parameter under DA for stratified populations
    psi[k] <- sum(lambda[1:nsites]) / 5000
    psiS[k] ~ dunif(0, 1)  
    
        for(i in 1:5000){                       # i is index for individuals
          z[i, k] ~ dbern(psi[k])                    # Data augmentation variables (are you real?!)
          x[i, k] ~ dunif(0, B)                   # distance uniformly distributed
          p[i, k] <- exp(-x[i, k] * x[i, k] / (2 * sigma[k] * sigma[k])) # Det. function
          mu[i, k] <- z[i, k] * p[i, k]                 # 'straw man' for WinBUGS
          y[i, k] ~ dbern(mu[i, k])                  # basic Bernoulli random variable
          site[i, k] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
      }
  }
  
  ## Models for abundance and for detection
  # Abundance dependent on VegType
  for(s in 1:nsites){                    # s is index for sites
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta[1] * VegTypeM[s, 1] + beta[2] * VegTypeM[s, 2] + beta[3] * VegTypeM[s, 3] +
      beta[4] * VegTypeM[s, 4] + beta[5] * VegTypeM[s, 5] + beta[6] * VegTypeM[s, 6] +
      beta[7] * VegTypeM[s, 7] + beta[8] * VegTypeM[s, 8] + beta[9] * VegTypeM[s, 9] +
      beta[10] * VegTypeM[s, 10] + beta[11] * VegTypeM[s, 11]
    
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
  }
  
  # Derived parameter: total population size across all sites
  Ntotal2019_W <- sum(z[1:5000, 1])
  Ntotal2019_D <- sum(z[1:5000, 2])
}
)

# Bundle data set
data <- list(x = x,
             y = y,
             site = site)
             #season = season)

consts <- list(nind = nind,
               nz = nz,
               B = B,
               nsites = nsites,
               nseasons = nseason,
               VegTypeM = VegTypeM)

# Inits
zst <- y
xinitW <- c(rep(NA, nindW), runif(nz - nindW, 0.1, 50))
xinitD <- c(rep(NA, nindD), runif(nz - nindD, 0.1, 50))
xinit <- cbind(xinitW, xinitD)
siteinitW <- c(rep(NA, nindW), sample(1:nsites, nz - nindW, replace = TRUE))
siteinitD <- c(rep(NA, nindD), sample(1:nsites, nz - nindD, replace = TRUE))
siteinit <- cbind(siteinitW, siteinitD)

inits <- list(beta = rep(0, 11),
              beta0 = 1,
              #alpha0 = c(100, 100),
              #alpha1 = dnorm(2, 0, 1),
              psiS = runif(2, 0, 1),
              z = zst,
              x = xinit, 
              N = rpois(nsites, exp(5)),
              site = siteinit,
              sigma = runif(2, 0, 100))

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
config$addMonitors(c("sigma", "psi", "psiS", "Ntotal2019_W", "Ntotal2019_D", "beta", "beta0"))
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

saveRDS(run, "Outputs/2019WD_VegSeason.rds")
run$summary
mcmc_trace(run$samples)

