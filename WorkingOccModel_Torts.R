# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
#library(mcmcplots)
library(lubridate)
library(bayesplot)

set.seed(15)

# Read in distance sampling data and remove entries with no transect ID and no distance data
torts <- read.csv("distancetorts.v2.csv", header = TRUE) %>%
  drop_na(TransectID) %>%
  drop_na(DistFromTransect)

# Remove entries left over from before 2018 and with distance >50m
torts <- torts %>%
  mutate(DateOfTransect = dmy(torts$DateOfTransect)) %>%
  filter(year(DateOfTransect) > 2018) %>%
  filter(DistFromTransect <= 50) %>% 
  filter(month(DateOfTransect) != 6) %>%
  filter(month(DateOfTransect) != 7) %>%
  select(DateOfTransect, TransectCode, VegType, DistFromTransect) %>%
  mutate(season = ifelse(month(DateOfTransect) < 6, "1", "2")) %>%
  mutate(year = year(DateOfTransect)) %>%
  unite("occasion", year, season, sep = "", remove = FALSE) %>%
  arrange(DateOfTransect) %>%
  mutate(occasion = as.numeric(as.factor(occasion))) %>%
  mutate(y = 1)

# plot data
ggplot(torts, aes(x = DistFromTransect, colour = season, fill = season)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~year)

# data augmentation
M <- 6000 # meta-population 

tortS <- torts %>%
  group_by(occasion) %>%
  count()
  
Obs <- tortS$n
Aug <- M - Obs

## get rid of unwanted info
torts <- select(torts, VegType, DistFromTransect, occasion, y)
names <- colnames(torts)

## occasion 1
torts1 <- filter(torts, occasion == 1)
torts1DA <- as.data.frame(matrix(nrow = Aug[1], ncol = ncol(torts)))
colnames(torts1DA) <- names
torts1DA$occasion <- 1
torts1DA$y <- 0
torts1 <- rbind(torts1, torts1DA)
rm(torts1DA)

## occasion 2
torts2 <- filter(torts, occasion == 2)
torts2DA <- as.data.frame(matrix(nrow = Aug[2], ncol = ncol(torts)))
colnames(torts2DA) <- names
torts2DA$occasion <- 2
torts2DA$y <- 0
torts2 <- rbind(torts2, torts2DA)
rm(torts2DA)

## occasion 3
torts3 <- filter(torts, occasion == 3)
torts3DA <- as.data.frame(matrix(nrow = Aug[3], ncol = ncol(torts)))
colnames(torts3DA) <- names
torts3DA$occasion <- 3
torts3DA$y <- 0
torts3 <- rbind(torts3, torts3DA)
rm(torts3DA)

## create matrices of info required
x <- data.frame(torts1$DistFromTransect, torts2$DistFromTransect, torts3$DistFromTransect)
colnames(x) <- c("2019W", "2019D", "2020W")

y <- data.frame(torts1$y, torts2$y, torts3$y)
colnames(y) <- c("2019W", "2019D", "2020W")

B <- 50 # Strip half-width. Larger than max distance
nseason <- 3

# nimble code
code <- nimbleCode({
  
  for(seas in 1:nseason){
    
    # Priors
    sigma[seas] ~ dunif(0, 1000)  # Half-normal scale
    psi[seas] ~ dunif(0, 1)       # DA parameter 
    
  for(i in 1:M){
    
      # Likelihood

        # Process model
        z[i, seas] ~ dbern(psi[seas])   # DA variables
        x[i, seas] ~ dunif(0, B)  # Distribution of distances
        
        # Observation model
        logp[i, seas] <- -((x[i, seas]*x[i, seas])/(2 * sigma[seas] * sigma[seas])) # Half-normal detection fct. - BOOK
        #logp[i] <- -(x[i] / sigma) # negative exponential
        p[i, seas] <- exp(logp[i, seas])
        mu[i, seas] <- z[i, seas] * p[i, seas]
        y[i, seas] ~ dbern(mu[i, seas]) # Simple Bernoulli measurement error process
      
    }
  }
  
  # Derived quantities
  N_2019W <- sum(z[1:M, 1])  # Population size
  N_2019D <- sum(z[1:M, 2])  # Population size
  N_2020W <- sum(z[1:M, 3])  # Population size

}
)

# Bundle data set
data <- list(x = x,
             y = y)

consts <- list(M = M,
               B = B,
               nseason = nseason)

# Inits
zst <- y
inits <- list(psi = runif(3),
              #z = zst,
              sigma = runif(3, 1, 100))

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
config$addMonitors(c("sigma", "psi", "N_2019W", "N_2019D", "N_2020W"))
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
saveRDS(run, "3Occ_Basic.rds")
run <- readRDS("Outputs/3Occ_Basic.rds")

mcmc_trace(run$samples, pars = c("sigma[1]", "psi[1]", "N_2019W", "N_2019D", "N_2020W"))
