#source("Chapter_8_setup.R")

# 8.3.1 Bayesian analysis of line transect data
# ------------------------------------------------------------------------
# Get data and do data-augmentation
# Observed distances (meters)

x <- torts$DistFromTransect

B <- 50 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
#nz <- 200 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind)) # Augmented inds. have y=0 by definition
#x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented

# Bundle and summarize data set
str( win.data <- list(nind=nind, y=y, B=B) )
code <- nimbleCode( {
  
  # Priors
  sigma ~ dunif(0,1000)  # Half-normal scale
  psi ~ dunif(0,1)       # DA parameter
  
  # Likelihood
  for(i in 1:(nind)){
    
    # Process model
    #z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    
    # Observation model
    logp[i] <- -((x[i] * x[i]) / (2 * sigma * sigma)) # Half-normal detection fct.
    p[i] <- exp(logp[i])
    
    mu[i] <- p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
  }
  # Derived quantities
  N <- sum(z[1:(nind)]) # Population size
  D <- N / 60                # Density, with A = 60 km^2 when B = 500
}
)

# Inits
zst <- y
inits <- function(){ list (psi=runif(1), 
                           z=zst, 
                           sigma=runif(1,40,200)) }

# Params to save
params <- c("N", "sigma", "D")

# Experience the raw power of BUGS and summarize marginal posteriors
out <- nimbleMCMC(code = code, 
                  constants = win.data, 
                  inits = inits,
                  monitors = params,
                  samplesAsCodaMCMC = TRUE,
                  nburnin = 1000,
                  niter = 11000)
