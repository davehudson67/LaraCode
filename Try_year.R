code <- nimbleCode({
  
  # Prior distributions
  # Intercept of lambda
  
  #for (t in 1:nSites){
  #  beta[t] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  #}
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  }
  
  for (year in 1:nyear){
    for (f in 1:nseason){
      
      # Intercept of log(sigma) for each season
      alpha0[f] ~ dunif(-10, 10) 
      
      # psi is a derived parameter under DA for stratified populations
      psi[year, f] <- sum(lambda[year, f, 1:nsites[f]]) / (nind[year, f] + nz[year, f])
      
      # 'Likelihood' (sort of...)
      for (i in 1:(nind[year, f] + nz[year, f])){               # i is index for individuals
        z[year, f, i] ~ dbern(psi[year, f])                    # Data augmentation variables
        d[year, f, i] ~ dunif(0, B)                      # distance uniformly distributed
        p[year, f, i] <- exp(-d[year, f, i] * d[year, f, i] / (2 * sigma[year, f, site[year, f, i]] * sigma[year, f, site[year, f, i]])) # Det. function
        mu[year, f, i] <- z[year, f, i] * p[year, f, i]              # 'straw man' for WinBUGS
        y[year, f, i] ~ dbern(mu[year, f, i])                  # basic Bernoulli random variable
      }
      
      for (s in 1:nsites[year, f]){                    # s is index for sites
        # next line not necessary, but allows to make predictions
        beta0[year, f, s] ~ dunif(0, 100)
        N[year, f, s] ~ dpois(lambda[year, f, s])              # Realized abundance at site s
        log(lambda[year, f, s]) <- beta0[year, s] #+ beta[1] * TransectA[f, 1, s] + beta[2] * TransectA[f, 2, s] + beta[3] * TransectA[f, 3, s] + beta[4] * TransectA[f, 4, s] +
        #beta[5] * TransectA[f, 5, s] + beta[6] * TransectA[f, 6, s] + beta[7] * TransectA[f, 7, s] + beta[8] * TransectA[f, 8, s] +
        #beta[9] * TransectA[f, 9, s] + beta[10] * TransectA[f, 10, s] + beta[11] * TransectA[f, 11, s] + beta[12] * TransectA[f, 12, s] +
        #beta[13] * TransectA[f, 13, s] + beta[14] * TransectA[f, 14, s] + beta[7] * TransectA[f, 15, s] + beta[8] * TransectA[f, 15, s] +
        #beta[16] * TransectA[f, 16, s] + beta[17] * TransectA[f, 17, s]
        
        # Linear model for detection
        log(sigma[year, f, s]) <- alpha0[f] + alpha[1] * VegTypeA[s, 1, f] + alpha[2] * VegTypeA[s, 2, f] + alpha[3] * VegTypeA[s, 3, f] +
          alpha[4] * VegTypeA[s, 4, f] + alpha[5] * VegTypeA[s, 5, f] + alpha[6] * VegTypeA[s, 6, f] +
          alpha[7] * VegTypeA[s, 7, f] + alpha[8] * VegTypeA[s, 8, f]
      }
    }
  }
  
  # Derived parameter: total population size across all site
  Ntotal2019_dry <- sum(z[1, 1, 1:6517])
  Ntotal2019_wet <- sum(z[1, 2, 1:6517])
  
})