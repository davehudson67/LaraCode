## load libraries
library(tidyverse)
library(nimble)

###########################################################
#######             ORIGINAL DA VERSION             #######
###########################################################

## load data
load("data.RData")

## amend vegetation type
VegTypeM <- apply(VegTypeM, 1, function(x) which(x == 1))

## extract maximum size at each site
M <- tapply(d, paste0(site, "_", season), length)
M <- tibble(nam = names(M), M = M) %>%
  separate(nam, c("site", "season"), sep = "_") %>%
  mutate(across(c(site, season), as.numeric)) %>%
  arrange(site, season)

## generate index matrix for observed site/season combinations 
site_inds <- matrix(NA, nsites, nseason)
for(i in 1:nrow(M)) {
  site_inds[M$site[i], M$season[i]] <- i
}

code <- nimbleCode({
  
  # 'Likelihood' (sort of...)
  for(i in 1:(nind+nz)){# i is index for individuals
    z[i] ~ dbern(psi[site_inds[site[i], season[i]]])                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i] * d[i] / (2 * sigma[site_inds[site[i], season[i]]] * sigma[site_inds[site[i], season[i]]])) # Det. function
    mu[i] <- z[i] * p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
  }
  
  # priors
  for(s in 1:nsite_ind) {
    lsigma[s] <- alpha[VegTypeM[sitea_ind[s]]]
    sigma[s] <- exp(lsigma[s])
    psi[s] ~ dunif(0, 1)
  }
  for(s in 1:nVegType) {
    alpha[s] ~ dnorm(0, sd = 100)
  }
  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind+nz)])
})

# Bundle and summarize data set
z <- numeric(0)
z[!is.na(d)] <- 1
z[is.na(d)] <- NA
data <- list(y = y, d = d, z = z)

consts <- list(
  B = B, 
  nind = nind, 
  nz = nz, 
  nsite_ind = length(M$M), 
  site = site, 
  season = season,
  sitea_ind = M$site,
  site_inds = site_inds,
  VegTypeM = VegTypeM,
  nVegType = nVegTypes)

## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
inits <- function(nVegTypes, nsite_ind, nind, nz, B, d, model) {
  valid <- 0
  while(valid == 0) {
    z <- rbinom(nind + nz, size = 1, prob = 0.5)
    d1 <- runif(nind + nz, 0, B)
    z[!is.na(d)] <- NA
    d1[!is.na(d)] <- NA
    inits1 <- list(
      alpha = runif(nVegTypes, 5, 10),
      psi = runif(nsite_ind, 0.1, 1),
      z = z,
      d = d1
    )
    ## check validity of initial values
    model$setInits(inits1)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  inits1
}
inits1 <- inits(nVegTypes, length(M$M), nind, nz, B, d, model)

# Parameters to save
params <- c("alpha", "psi", "sigma", "Ntotal", "z", "mu")

# MCMC settings
ni <- 20000
nb <- 10000
nt <- 2
nc <- 3

# Call nimble from R
## took 1738.192 secs
system.time(
  out1 <- nimbleMCMC(code = code, 
    data = data,
    constants = consts, 
    inits = inits1, 
    monitors = params,
    nburnin = nb, 
    niter = ni,
    samplesAsCodaMCMC = TRUE)
)

## plot traces
plot(out1[, "Ntotal"])
plot(out1[, grep("sigma", colnames(out1))])
plot(out1[, grep("psi", colnames(out1))])
plot(out1[, grep("alpha", colnames(out1))])

## extract posterior samples for number of individuals in each 
## site/season combination
N <- apply(out1[, grep("z", colnames(out1))], 1, function(z, g) {
  tapply(z, g, sum)
}, g = paste0(site, "_", season))

## reorder rows to match M
inds <- match(paste0(M$site, "_", M$season), rownames(N))
stopifnot(all((sort(inds) - 1:nrow(M)) == 0))
N <- N[inds, ]
stopifnot(identical(paste0(M$site, "_", M$season), rownames(N)))

## check whether number of augmented variables is high enough
## (after convergence these should all be TRUE, else the required
## number of alive animals bounces up against the maximum set. Note
## that this would be impractical if using the full DA version I 
## suspect. If you don't do this then it underestimates the
## population size).
table(apply(cbind(M$M, N), 1, function(x) all(x[-1] < x[1])))

## summarise population sizes
N <- apply(N, 1, function(x) {
  c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
})
N <- t(N)
colnames(N) <- c("mean", "LCI", "UCI")
N <- as_tibble(N) %>%
  mutate(site = M$site) %>%
  mutate(season = M$season)

## generate posterior predictive summaries for number detected
Ny <- list()
for(i in 1:(nind + nz)) {
  Ny[[i]] <- rbinom(
    nrow(out1), 
    size = out1[, paste0("z[", i, "]")], 
    prob = out1[, paste0("mu[", i, "]")]
  )
}
Ny <- do.call("cbind", Ny)
Ny <- apply(Ny, 1, function(z, g) {
  tapply(z, g, sum)
}, g = paste0(site, "_", season))

## reorder rows to match M
inds <- match(paste0(M$site, "_", M$season), rownames(Ny))
stopifnot(all((sort(inds) - 1:nrow(M)) == 0))
Ny <- Ny[inds, ]
stopifnot(identical(paste0(M$site, "_", M$season), rownames(Ny)))

## extract posterior summaries for number of detected animals in each
## site / season combination
Ny <- apply(Ny, 1, function(x){
  c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
})
Ny <- t(Ny)
colnames(Ny) <- c("mean", "LCI", "UCI")
Ny <- as_tibble(Ny) %>%
  mutate(site = M$site) %>%
  mutate(season = M$season)

## observed data
Ninit <- cbind(d, paste0(site, "_", season))
Ninit <- Ninit[!is.na(Ninit[, 1]), ]
Ninit <- tapply(Ninit[, 1], Ninit[, 2], length)
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
  inner_join(Ny, by = c("site", "season")) %>%
  mutate(inside = value > LCI & value < UCI) %>%
  ggplot() +
  geom_point(aes(x = site, y = mean)) +
  geom_errorbar(aes(x = site, ymin = LCI, ymax = UCI)) +
  geom_point(aes(x = site, y = value, colour = inside)) +
  facet_wrap(~season, ncol = 1)
p

############################################################
#######              FASTER VERSION                  #######
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

## amend dummy variables (this should remove extraneous column
## and shift dummy variables)
VegTypeM <- apply(VegTypeM, 1, function(x) which(x == 1))

## augment at each site / season
aug <- 20
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
## takes 65 seconds, so speed-up
## of 1738 / 65 = 27 times faster
system.time(
  out2 <- nimbleMCMC(code = code, 
    data = data,
    constants = consts, 
    inits = inits1, 
    monitors = params,
    nburnin = nb, 
    niter = ni,
    samplesAsCodaMCMC = TRUE)
)

## trace plots for different subsets
plot(out2[, "Ntotal"])
plot(out2[, grep("N\\[", colnames(out2))])
plot(out2[, grep("sigma\\[", colnames(out2))])
plot(out2[, grep("psi", colnames(out2))])
plot(out2[, grep("alpha", colnames(out2))])

## check whether number of augmented variables is high enough
## (after convergence these should all be TRUE, else the required
## number of alive animals bounces up against the maximum set. Note
## that this would be impractical if using the full DA version I 
## suspect. If you don't do this then it underestimates the
## population size).
table(apply(rbind(M$M, out2[, grep("N\\[", colnames(out2))]), 2, function(x) all(x[-1] < x[1])))

## generate posterior predictive summaries for number detected
Ny <- list()
for(i in 1:nrow(M)) {
  Ny[[i]] <- rbinom(
    nrow(out2), 
    size = out2[, paste0("N[", i, "]")], 
    prob = out2[, paste0("pbar[", i, "]")]
  )
}
Ny <- do.call("cbind", Ny)

## extract posterior summaries for number of detected animals in each
## site / season combination
Ny <- apply(Ny, 2, function(x){
  c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
})
Ny <- t(Ny)
colnames(Ny) <- c("mean", "LCI", "UCI")
Ny <- as_tibble(Ny) %>%
  mutate(Ny, site = M$site) %>%
  mutate(season = M$season)

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
  inner_join(Ny, by = c("site", "season")) %>%
  mutate(inside = value > LCI & value < UCI) %>%
  ggplot() +
    geom_point(aes(x = site, y = mean)) +
    geom_errorbar(aes(x = site, ymin = LCI, ymax = UCI)) +
    geom_point(aes(x = site, y = value, colour = inside)) +
    facet_wrap(~season, ncol = 1)
p

##################################################
#######            COMPARISONS             #######
##################################################

## compare alphas
summary(out1[, grep("alpha", colnames(out1))])
summary(out2[, grep("alpha", colnames(out2))])
p <- as_tibble(out1[, grep("alpha", colnames(out1))]) %>%
  mutate(model = "DA") %>%
  rbind(mutate(as_tibble(out2[, grep("alpha", colnames(out2))]), model = "Faster")) %>%
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

## compare population sizes in each site/season combination
p <- apply(out2[, grep("N\\[", colnames(out2))], 2, function(x) {
      c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }) %>%
  t() %>%
  as.data.frame() %>%
  set_names(c("mean", "LCI", "UCI")) %>%
  mutate(site = M$site) %>%
  mutate(season = M$season) %>%
  mutate(model = "Faster") %>%
  rbind(mutate(N, model = "DA")) %>%
  ggplot(aes(x = site, colour = model)) +
    geom_point(aes(y = mean)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI)) +
    facet_wrap(~season)
p

## compare population sizes in each site/season combination
p <- apply(out2[, grep("N\\[", colnames(out2))], 2, function(x) {
    c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }) %>%
  t() %>%
  as.data.frame() %>%
  set_names(c("mean", "LCI", "UCI")) %>%
  mutate(site = M$site) %>%
  mutate(season = M$season) %>%
  mutate(model = "Faster") %>%
  rbind(mutate(N, model = "DA")) %>%
  select(mean, site, season, model) %>%
  pivot_wider(names_from = model, values_from = mean) %>%
  ggplot(aes(x = DA, y = Faster)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
p

## compare overall population sizes
p <- rbind(
    tibble(N = out1[, "Ntotal"], model = "DA"),
    tibble(N = out2[, "Ntotal"], model = "Faster")
  ) %>%
  ggplot() +
    geom_density(aes(x = N, fill = model), alpha = 0.3)
  

