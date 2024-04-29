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

tortsW2018 <- torts %>%
  filter(Year == 2018)
# 8.3.1 Bayesian analysis of line transect data
# ------------------------------------------------------------------------
# Get data and do data-augmentation
# Observed distances (meters)

x <- tortsW2018$DistFromTransect
B <- 60 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
nz <- 200 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented

# Analysis of binned data using data augmentation
delta <- 5                # Width of distance bins
xg <- seq(0, B, delta)     # Make the interval cut points
dclass <- x %/% delta + 1  # Convert distances to distance category
nD <- length(xg) -1        # N intervals = length(xg) if max(x) = B

# Bundle data
# Note data changed to include dclass, nG, bin-width delta and midpt
midpt <- xg[-1] - delta/2  # Interval mid-points
str( win.data <- list (nind=nind, nz=nz, dclass=dclass, y=y, B=B,
                       delta=delta, nD=nD, midpt=midpt) )   # Bundle and summarize

# BUGS model specification
# This code is called "model2.txt" in the book code.
Section8p3p1_code_model2 <- nimbleCode({ 
  # Priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 1000)
  
  # Likelihood
  # construct conditional detection probability and Pr(x) for each bin
  for(g in 1:nD){        # midpt = mid point of each cell
    log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
    pi[g] <- delta / B  # probability of x in each interval
  }
  
  for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)             # model for individual covariates
    dclass[i] ~ dcat(pi[1:nD])    # population distribution of distance class
    mu[i] <- z[i] * p[dclass[i]]  # p depends on distance class
    y[i] ~ dbern(mu[i])
  }
  # Derived quantities: Population size and density
  N <- sum(z[1:(nind+nz)])
  D <- N / 60
}
)

# Inits function
zst <- y # DA variables start at observed value of y
inits <- function(){ list (psi=runif(1), 
                           z=zst, 
                           sigma=runif(1,40,200),
                           dclass = c(
                             rep(as.numeric(NA), nind),
                             sample(1:win.data$nD, win.data$nz, replace = TRUE))
) }

# Parameters to save
params <- c("N", "sigma", "D")

# Unleash WinBUGS and summarize posteriors
out.1 <- nimbleMCMC(code = Section8p3p1_code_model2,
                    nchains = 2,
                    constants = win.data, 
                    inits = inits,
                    monitors = params,
                    nburnin = 1000,
                    niter = 11000,
                    samplesAsCodaMCMC = TRUE)

mcmcplot(out.1)
