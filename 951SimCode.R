library(nimble)
library(bayesplot)

code <- nimbleCode({ # point count with DA
  # Prior distributions
  beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
  mean.lam <- exp(beta0)
  beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
  phi ~ dunif(0,1)        # Probability of availability
  sigma ~ dunif(0.01,5)   # Distance function parameter
  
  # Detection probs for each distance interval and related things
  for(b in 1:nD){ # loops over each distance interval
    
    # detection prob for each dist. interval using the half-normal detection fct.
    log(g[b]) <- -midpt[b] * midpt[b] / (2 * sigma * sigma)     
    
    # radial density function - how the density changes with distance conditional cell prob. for each distance interval
    f[b] <- (2 * midpt[b] * delta) / (B * B)                   
    
    # This line calculates the cell probability for each distance interval
    cellprobs[b] <- g[b] * f[b] 

    # Conditional cell probability - the prob. of detecting an animal at that distance, given that it was detected 
    #   somewhere within the range of distance intervals.
    cellprobs.cond[b] <- cellprobs[b] / sum(cellprobs[1:nD])
  }
  
    # This line calculates the cell probability for the distance interval beyond the last interval.
    #   (i.e., outside the range of nD intervals).
    cellprobs[nD + 1] <- 1 - sum(cellprobs[1:nD])                
                                                               
  for (s in 1:nsites) { # loops over sites
    
    for (k in 1:K) { # loops over replicate survey occasions
      
      # prob. of detection 
      pdet[s, k] <- sum(cellprobs[1:nD])   # Distance class probabilities
      
      # marginal prob.
      pmarg[s, k] <- pdet[s, k] * phi      # Marginal probability
      
      # Model part 4: distance class frequencies
      # multivariate observation model
      y3d[s, 1:nD, k] ~ dmulti(cellprobs.cond[1:nD], nobs[s, k])
      
      # Model part 3: total number of detections:
      # counts of detections follow binomial distribution
      nobs[s, k] ~ dbin(pmarg[s, k], M[s])
      # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
     
      # Model part 2: Availability. Not used in this model but simulated.
      Navail[s,k] ~ dbin(phi, M[s]) 
    }  # end k loop

    # Model part 1: Abundance model
    M[s] ~ dpois(lambda[s])    
    log(lambda[s]) <- beta0 + beta1 * habitat[s]
  }  # End s loop
  
  
  # Derived quantities
  Mtot <- sum(M[1:nsites])
  for(k in 1:K){ 
    Ntot[k] <- sum(Navail[1:nsites, k])
  }
} # End model
)

y3d <- y4d[, , , 1]

## Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion
data <- list(y3d = y3d, habitat = Sim$habitat, B = B, nobs = nobs)
consts <- list(nsites = nsites, K = K, nD = nD, midpt = midpt, delta = delta)

# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max) + 2
inits <- list(M=Mst, sigma = 1.0, phi=0.9, beta0=log(2), beta1=0.5)


params <- c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 5   ;   nc <- 3

# Build model
model <- nimbleModel(code, constants = consts, data = data, inits = inits)
model$initializeInfo()
#warnings()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
config$removeSamplers(c("beta0", "beta1"))
config$addSampler(target = c("beta0", "beta1"), type = 'AF_slice')
#config$addSampler(target = c("M"), type = 'AF_slice')


config$resetMonitors()
config$addMonitors(c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot"))
#config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 1000000, 
                           nburnin = 50000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run$summary


mcmc_trace(run$samples)


