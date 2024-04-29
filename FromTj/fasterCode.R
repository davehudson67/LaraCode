## load libraries
library(tidyverse)
library(nimble)

###########################################################
#######              FASTER VERSION                 #######
###########################################################

## load data
load("data.RData")

## Extract data from DA version and amend to non-DA version
site <- site[!is.na(d)]
season <- season[!is.na(d)]
d <- d[!is.na(d)]

## check all sites present
stopifnot(all((sort(unique(site)) - 1:nsites) == 0))
y <- rep(1, nind)
B <- 50 #Max observation distance
nseason <- 2
stopifnot(all((sort(unique(season)) - 1:nseason) == 0))

## amend dummy variables (this should remove extraneous column
## and shift dummy variables)
VegTypeM <- VegTypeM[, -1]
nVegTypes <- nVegTypes - 1

## augment at each site / season
aug <- 200
M <- tapply(d, paste0(site, "_", season), length) + aug
M <- data.frame(nam = names(M), M = M) %>%
  separate(nam, c("site", "season"), sep = "_") %>%
  mutate(across(c(site, season), as.numeric)) %>%
  arrange(site, season)

## generate index matrix for missing site/season combinations 
site_inds <- matrix(NA, nsites, nseason)
for(i in 1:nrow(M)) {
  site_inds[M$site[i], M$season[i]] <- i
}

## model code
code <- nimbleCode({

    ## positive constant so zeros trick works
    C <- 10000

    ## loop over observed individuals
    for(i in 1:nind) {
        ## -(log-likelihood) + constant (C must be chosen so that rate > 0 always)
        rate[i] <- -((-d[i] * d[i] / (2 * sigma[site_inds[site[i], season[i]]] * sigma[site_inds[site[i], season[i]]])) - log(pbar[site_inds[site[i], season[i]]]) - log(B)) + C
        y[i] ~ dpois(rate[i])
    }

    # loop over sites and seasons
    for(s in 1:nsite_ind) {
        ## likelihood terms
        n[s] ~ dbin(pbar[s], N[s])
        N[s] ~ dbin(psi[s], M[s])
        
        ## integrated detection probability
        pbar[s] <- sqrt_2_pi * sigma[s] * (pnorm(B, mean = 0, sd = sigma[s]) - 0.5) / B
        
        ## priors
        lsigma[s] <- alpha0[seasona_ind[s]] + alpha1[1] * VegTypeM[sitea_ind[s], 1] + 
          alpha1[2] * VegTypeM[sitea_ind[s], 2] + alpha1[3] * VegTypeM[sitea_ind[s], 3] +
          alpha1[4] * VegTypeM[sitea_ind[s], 4] + alpha1[5] * VegTypeM[sitea_ind[s], 5] +
          alpha1[6] * VegTypeM[sitea_ind[s], 6] + alpha1[7] * VegTypeM[sitea_ind[s], 7]
        sigma[s] <- exp(lsigma[s])
        psi[s] ~ dunif(0, 1)
    }
    for(s in 1:nVegType) {
        alpha1[s] ~ dnorm(0, sd = 100)
    }
    for(s in 1:nseason) {
        alpha0[s] ~ dnorm(0, sd = 100)
    }
    
    # Derived parameter: total population size across all sites
    Ntotal <- sum(N[1:nsite_ind])
})

# Bundle and summarise data set
data <- list(y = rep(0, length(d)), d = d, n = M$M - aug)

consts <- list(
    B = B, 
    nind = nind, 
    nsite_ind = length(M$M), 
    site = site, 
    season = season,  
    nseason = nseason, 
    M = M$M,
    seasona_ind = M$season,
    sitea_ind = M$site,
    site_inds = site_inds,
    sqrt_2_pi = sqrt(2 * pi),
    VegTypeM = VegTypeM,
    nVegType = nVegTypes)
    
## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
inits <- function(nVegType, nseason, M, aug, model) {
    valid <- 0
    while(valid == 0) {
        inits1 <- list(
            alpha0 = rnorm(nseason, 0, sd = 1),
            alpha1 = rnorm(nVegType, 0, sd = 1),
            psi = runif(length(M), 0, 1),
            N = round(runif(length(M), M - aug, M))
        )
        ## check validity of initial values
        model$setInits(inits1)
        valid <- ifelse(!is.finite(model$calculate()), 0, 1)
    }
    inits1
}
inits1 <- inits(nVegTypes, nseason, M$M, aug, model)
    
# Parameters to save
params <- c("sigma", "psi", "N", "Ntotal", "alpha0", "alpha1")

# MCMC settings
ni <- 20000
nb <- 10000
nt <- 2
nc <- 3

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
    data = data,
    constants = consts, 
    inits = inits1, 
    monitors = params,
    nburnin = nb, 
    niter = ni,
    samplesAsCodaMCMC = TRUE)

## trace plots for different subsets
plot(out1[, "Ntotal"])
plot(out1[, grep("N\\[", colnames(out1))])
plot(out1[, grep("sigma\\[", colnames(out1))])
plot(out1[, grep("psi", colnames(out1))])
plot(out1[, grep("alpha", colnames(out1))])
    
## check whether number of augmented variables is high enough
## (after convergence these should all be TRUE, else the required
## number of alive animals bounces up against the maximum set. Note
## that this would be impractical if using the full DA version I 
## suspect. If you don't do this then it underestimates the
## population size).
table(apply(rbind(M$M, out1[, grep("N\\[", colnames(out1))]), 2, function(x) all(x[-1] < x[1])))

## correlation between the alphas
cor(out1[, grep("alpha", colnames(out1))])

###########################################################
#######    FASTER VERSION WITH SINGLE INTERCEPT     #######
###########################################################

## load data
load("data.RData")

## Extract data from DA version and amend to non-DA version
site <- site[!is.na(d)]
season <- season[!is.na(d)]
d <- d[!is.na(d)]

## check all sites present
stopifnot(all((sort(unique(site)) - 1:nsites) == 0))
y <- rep(1, nind)
B <- 50 #Max observation distance
nseason <- 2
stopifnot(all((sort(unique(season)) - 1:nseason) == 0))

## amend dummy variables (this should remove extraneous column
## and shift dummy variables)
VegTypeM <- VegTypeM[, -1]
nVegTypes <- nVegTypes - 1

## augment at each site / season
aug <- 200
M <- tapply(d, paste0(site, "_", season), length) + aug
M <- data.frame(nam = names(M), M = M) %>%
  separate(nam, c("site", "season"), sep = "_") %>%
  mutate(across(c(site, season), as.numeric)) %>%
  arrange(site, season)

## generate index matrix for missing site/season combinations 
site_inds <- matrix(NA, nsites, nseason)
for(i in 1:nrow(M)) {
  site_inds[M$site[i], M$season[i]] <- i
}


## model code
code <- nimbleCode({
  
  ## positive constant so zeros trick works
  C <- 10000
  
  ## loop over observed individuals
  for(i in 1:nind) {
    ## -(log-likelihood) + constant (C must be chosen so that rate > 0 always)
    rate[i] <- -((-d[i] * d[i] / (2 * sigma[site_inds[site[i], season[i]]] * sigma[site_inds[site[i], season[i]]])) - log(pbar[site_inds[site[i], season[i]]]) - log(B)) + C
    y[i] ~ dpois(rate[i])
  }
  
  # loop over sites and seasons
  for(s in 1:nsite_ind) {
    ## likelihood terms
    n[s] ~ dbin(pbar[s], N[s])
    N[s] ~ dbin(psi[s], M[s])
    
    ## integrated detection probability
    pbar[s] <- sqrt_2_pi * sigma[s] * (pnorm(B, mean = 0, sd = sigma[s]) - 0.5) / B
    
    ## priors
    lsigma[s] <- alpha0 + alpha1[1] * VegTypeM[sitea_ind[s], 1] + 
      alpha1[2] * VegTypeM[sitea_ind[s], 2] + alpha1[3] * VegTypeM[sitea_ind[s], 3] +
      alpha1[4] * VegTypeM[sitea_ind[s], 4] + alpha1[5] * VegTypeM[sitea_ind[s], 5] +
      alpha1[6] * VegTypeM[sitea_ind[s], 6] + alpha1[7] * VegTypeM[sitea_ind[s], 7]
    sigma[s] <- exp(lsigma[s])
    psi[s] ~ dunif(0, 1)
  }
  alpha0 ~ dnorm(0, sd = 100)
  for(s in 1:nVegType) {
    alpha1[s] ~ dnorm(0, sd = 100)
  }
  
  # Derived parameter: total population size across all sites
  Ntotal <- sum(N[1:nsite_ind])
})

# Bundle and summarise data set
data <- list(y = rep(0, length(d)), d = d, n = M$M - aug)

consts <- list(
  B = B, 
  nind = nind, 
  nsite_ind = length(M$M), 
  site = site, 
  season = season,
  M = M$M,
  sitea_ind = M$site,
  site_inds = site_inds,
  sqrt_2_pi = sqrt(2 * pi),
  VegTypeM = VegTypeM,
  nVegType = nVegTypes)

## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
inits <- function(nVegType, M, aug, model) {
  valid <- 0
  while(valid == 0) {
    inits1 <- list(
      alpha0 = rnorm(1, 0, sd = 1),
      alpha1 = rnorm(nVegType, 0, sd = 1),
      psi = runif(length(M), 0, 1),
      N = round(runif(length(M), M - aug, M))
    )
    ## check validity of initial values
    model$setInits(inits1)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  inits1
}
inits1 <- inits(nVegTypes, M$M, aug, model)

# Parameters to save
params <- c("sigma", "psi", "N", "Ntotal", "alpha0", "alpha1")

# MCMC settings
ni <- 100000
nb <- 10000
nt <- 2
nc <- 3

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
  data = data,
  constants = consts, 
  inits = inits1, 
  monitors = params,
  nburnin = nb, 
  niter = ni,
  samplesAsCodaMCMC = TRUE)

## trace plots for different subsets
plot(out1[, "Ntotal"])
plot(out1[, grep("N\\[", colnames(out1))])
plot(out1[, grep("sigma\\[", colnames(out1))])
plot(out1[, grep("psi", colnames(out1))])
plot(out1[, grep("alpha", colnames(out1))])

## save for later comparison
save_out1 <- out1[, grep("alpha", colnames(out1))]

## check whether number of augmented variables is high enough
## (after convergence these should all be TRUE, else the required
## number of alive animals bounces up against the maximum set. Note
## that this would be impractical if using the full DA version I 
## suspect. If you don't do this then it underestimates the
## population size).
table(apply(rbind(M$M, out1[, grep("N\\[", colnames(out1))]), 2, function(x) all(x[-1] < x[1])))

## correlation between the alphas
cor(out1[, grep("alpha", colnames(out1))])

############################################################
#######      FASTER VERSION WITH NO INTERCEPT       #######
############################################################

## load data
load("data.RData")

## Extract data from DA version and amend to non-DA version
site <- site[!is.na(d)]
season <- season[!is.na(d)]
d <- d[!is.na(d)]

## check all sites present
stopifnot(all((sort(unique(site)) - 1:nsites) == 0))
y <- rep(1, nind)
B <- 50 #Max observation distance
nseason <- 2
stopifnot(all((sort(unique(season)) - 1:nseason) == 0))

## amend VegTypeM
VegTypeM <- apply(VegTypeM, 1, function(x) which(x == 1))

## augment at each site / season
aug <- 200
M <- tapply(d, paste0(site, "_", season), length) + aug
M <- data.frame(nam = names(M), M = M) %>%
  separate(nam, c("site", "season"), sep = "_") %>%
  mutate(across(c(site, season), as.numeric)) %>%
  arrange(site, season)

## generate index matrix for missing site/season combinations 
site_inds <- matrix(NA, nsites, nseason)
for(i in 1:nrow(M)) {
  site_inds[M$site[i], M$season[i]] <- i
}

## model code
code <- nimbleCode({
  
  ## positive constant so zeros trick works
  C <- 10000
  
  ## loop over observed individuals
  for(i in 1:nind) {
    ## -(log-likelihood) + constant (C must be chosen so that rate > 0 always)
    rate[i] <- -((-d[i] * d[i] / (2 * sigma[site_inds[site[i], season[i]]] * sigma[site_inds[site[i], season[i]]])) - log(pbar[site_inds[site[i], season[i]]]) - log(B)) + C
    y[i] ~ dpois(rate[i])
  }
  
  # loop over sites and seasons
  for(s in 1:nsite_ind) {
    ## likelihood terms
    n[s] ~ dbin(pbar[s], N[s])
    N[s] ~ dbin(psi[s], M[s])
    
    ## integrated detection probability
    pbar[s] <- sqrt_2_pi * sigma[s] * (pnorm(B, mean = 0, sd = sigma[s]) - 0.5) / B
    
    ## priors
    lsigma[s] <- alpha[VegTypeM[sitea_ind[s]]]
    sigma[s] <- exp(lsigma[s])
    psi[s] ~ dunif(0, 1)
  }
  for(s in 1:nVegType) {
    alpha[s] ~ dnorm(0, sd = 100)
  }
  
  # Derived parameter: total population size across all sites
  Ntotal <- sum(N[1:nsite_ind])
})

# Bundle and summarise data set
data <- list(y = rep(0, length(d)), d = d, n = M$M - aug)

consts <- list(
  B = B, 
  nind = nind, 
  nsite_ind = length(M$M), 
  site = site, 
  season = season,
  M = M$M,
  sitea_ind = M$site,
  site_inds = site_inds,
  sqrt_2_pi = sqrt(2 * pi),
  VegTypeM = VegTypeM,
  nVegType = nVegTypes)

## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
inits <- function(nVegType, M, aug, model) {
  valid <- 0
  while(valid == 0) {
    inits1 <- list(
      alpha = rnorm(nVegType, 0, sd = 1),
      psi = runif(length(M), 0, 1),
      N = round(runif(length(M), M - aug, M))
    )
    ## check validity of initial values
    model$setInits(inits1)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  inits1
}
inits1 <- inits(nVegTypes, M$M, aug, model)

# Parameters to save
params <- c("sigma", "psi", "N", "Ntotal", "alpha", "pbar")

# MCMC settings
ni <- 20000
nb <- 10000
nt <- 2
nc <- 3

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
  data = data,
  constants = consts, 
  inits = inits1, 
  monitors = params,
  nburnin = nb, 
  niter = ni,
  samplesAsCodaMCMC = TRUE)

## trace plots for different subsets
plot(out1[, "Ntotal"])
plot(out1[, grep("N\\[", colnames(out1))])
plot(out1[, grep("sigma\\[", colnames(out1))])
plot(out1[, grep("psi", colnames(out1))])
plot(out1[, grep("alpha", colnames(out1))])

## compare to earlier model
for(j in 1:7) save_out1[, paste0("alpha1[", j, "]")] <- save_out1[, "alpha0"] + save_out1[, paste0("alpha1[", j, "]")]
summary(save_out1)
summary(out1[, grep("alpha", colnames(out1))])
colnames(save_out1) <- paste0("alpha[", 1:8, "]")
p <- as_tibble(save_out1) %>%
  mutate(model = 1) %>%
  rbind(mutate(as_tibble(out1[, grep("alpha", colnames(out1))]), model = 2)) %>%
  group_by(model) %>%
  summarise(across(everything(), list(
    mn = ~mean(.), 
    LCI = ~quantile(., probs = 0.025), 
    UCI = ~quantile(., probs = 0.975)
  ))) %>%
  pivot_longer(!model) %>%
  separate(name, c("var", "type"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(model = factor(model)) %>%
  ggplot() +
    geom_point(aes(x = var, y = mn, colour = model)) +
    geom_errorbar(aes(x = var, ymin = LCI, ymax = UCI, colour = model))
p

## check whether number of augmented variables is high enough
## (after convergence these should all be TRUE, else the required
## number of alive animals bounces up against the maximum set. Note
## that this would be impractical if using the full DA version I 
## suspect. If you don't do this then it underestimates the
## population size).
table(apply(rbind(M$M, out1[, grep("N\\[", colnames(out1))]), 2, function(x) all(x[-1] < x[1])))

## correlation between the alphas
cor(out1[, grep("alpha", colnames(out1))])

## how many observations in each site
tapply(tapply(d, site, length), VegTypeM, sum)

## check model fits using posterior predictive distributions
out2 <- as.matrix(out1) %>%
  as.data.frame() %>%
  select(starts_with("pbar") | starts_with("N")) %>%
  as.matrix()
out3 <- list()
for(i in 1:nrow(M)) {
  temp <- out2[, paste0("pbar[", i, "]")]
  tempN <- out2[, paste0("N[", i, "]")]
  temp <- rbinom(length(temp), size = tempN, prob = temp)
  out3[[i]] <- data.frame(
    mean = mean(temp), 
    LCI = quantile(temp, prob = 0.025), 
    UCI = quantile(temp, prob = 0.975)
  )
}
out3 <- bind_rows(out3, .id = "site") %>%
  mutate(site = as.numeric(site)) %>%
  mutate(season = M$season[site]) %>%
  mutate(site = M$site[site])
rownames(out3) <- NULL

## data
Ninit <- tapply(d, paste0(site, "_", season), length)
n <- matrix(NA, nsites, nseason)
Ninit_nam <- strsplit(names(Ninit), "_")
for(i in 1:length(Ninit_nam)) {
  temp <- as.numeric(Ninit_nam[[i]])
  n[temp[1], temp[2]] <- Ninit[i]
}
p <- as.data.frame(n) %>%
  set_names(c("s1", "s2")) %>%
  mutate(site = 1:n()) %>%
  pivot_longer(!site, names_to = "season") %>%
  mutate(season = as.numeric(gsub("s", "", season))) %>%
  filter(!is.na(value)) %>%
  inner_join(out3, by = c("site", "season")) %>%
  mutate(inside = value > LCI & value < UCI) %>%
  ggplot() +
    geom_point(aes(x = site, y = mean)) +
    geom_errorbar(aes(x = site, ymin = LCI, ymax = UCI)) +
    geom_point(aes(x = site, y = value, colour = inside)) +
    facet_wrap(~season, ncol = 1)
p
