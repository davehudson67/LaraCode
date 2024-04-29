############################################################
#######      FASTER VERSION WITH INTERCEPT           #######
############################################################

# Packages and data
library(tidyverse)
library(nimble)
library(mcmcplots)

## load formatted data
torts <- readRDS("Tortsi2018.rds")
head(torts)

torts <- torts %>%
  filter(Year != 2018) %>%
  droplevels() %>%
  mutate(SectionCode = as.factor(as.numeric(SectionCode))) %>%
  mutate(site = as.numeric(SectionCode)) %>%
  mutate(season = as.numeric(as.factor(season))) %>%
  mutate(year = as.numeric(as.factor(Year))) %>%
  mutate(Occ = as.numeric(as.factor((occasion))))

## Adjust veg type levels
levels(torts$VegType) <- sapply(levels(torts$VegType), function(x) {
  if (x == "CRB") "EXP" else if (x == "CSB") "SND" else if (x == "MSC") "SMX" else x
})

## Extract data from DA version and amend to non-DA version
site <- as.numeric(torts$SectionCode)
season <- as.numeric(as.factor(torts$season))
d <- torts$DistFromTransect
nsites <- length(unique(site))
nind <- nrow(torts)
VegType_vector <- as.numeric(torts$VegType)
nVegType <- length(unique(torts$VegType))
year <- as.numeric(as.factor(torts$Year))
nyears <- length(unique(year))
noccasions <- length(unique(torts$Occ))
occasion <- torts$Occ
transect <- torts %>%
  mutate(TransectCode = as.numeric(as.factor(TransectCode))) %>%
  distinct(SectionCode, .keep_all = TRUE) %>%
  select(SectionCode, TransectCode) %>%
  arrange(SectionCode)
nTransects <- length(unique(torts$TransectCode))
VegType <- torts %>%
  mutate(VegType = as.numeric(VegType)) %>%
  distinct(SectionCode, .keep_all = TRUE) %>%
  select(SectionCode, VegType) %>%
  arrange(SectionCode)

VegType$VegType <- as.factor(VegType$VegType)
VegTypeM <- model.matrix(~ VegType - 1, data = VegType)
VegTypeM <- VegTypeM[, -1]
nVegType <- nVegType - 1

## check all sites present
stopifnot(all((sort(unique(site)) - 1:nsites) == 0))
y <- rep(1, nind)
B <- 50 #Max observation distance
nseasons <- 2
stopifnot(all((sort(unique(occasion)) - 1:noccasions) == 0))

## augment at each site / season
aug <- 200
M <- tapply(d, paste0(site, "_", season, "_", year, "_", occasion), length) + aug
M <- data.frame(nam = names(M), M = M) %>%
  separate(nam, c("site", "season", "year", "occasion"), sep = "_") %>%
  mutate(across(c(site, season, year, occasion), as.numeric)) %>%
  arrange(site, occasion)

## generate index matrix for missing site/season combinations 
site_inds <- matrix(NA, nsites, noccasions)
for(i in 1:nrow(M)) {
  site_inds[M$site[i], M$occasion[i]] <- i
}

## generate index matrix for missing site/year combinations 
M$SiteYear <- with(M, paste(site, year, sep ="_"))
SiteYear <- M %>%
  distinct(site, year) %>%
  arrange(site, year) %>%
  mutate(global_seq = row_number()) %>%
  select(site, year, global_seq) %>%
  pivot_wider(names_from = year, values_from = global_seq)

SiteYear[is.na(SiteYear)] <- NA
SiteYear <- as.matrix(SiteYear[,-1])  # Exclude the site column

## map n to N - allow for different N for each year
ntoN <- M %>%
  distinct(site, year) %>%
  arrange(site, year) %>%
  mutate(seq <- row_number()) %>%
  rename(seq = 'seq <- row_number()')

ntoN <- M %>%
  left_join(ntoN %>%
              select(site, year, seq), by = c("site", "year"))

## Convert the relevant data (site, season, year, global_seq) into a matrix format
ntoN <- c(with(ntoN, matrix(seq, nrow = nrow(ntoN), ncol = 1)))
M$ntoN <- ntoN

##
M <- M %>%
  group_by(site, year) %>%
  mutate(maxMsite = max(M, na.rm = TRUE)) %>%
  ungroup()

maxMsite <- M %>%
  group_by(site) %>%
  summarise(max = max(M, na.rm = TRUE))
maxMsite <- maxMsite$max

## model code
code <- nimbleCode({
  
  ## positive constant so zeros trick works
  C <- 10000
  
  ## loop over observed individuals
  for(i in 1:nind) {
    ## -(log-likelihood) + constant (C must be chosen so that rate > 0 always)
    rate[i] <- -((-d[i] * d[i] / (2 * sigma[site_inds[site[i], occasion[i]]] * sigma[site_inds[site[i], occasion[i]]])) - log(pbar[site_inds[site[i], occasion[i]]]) - log(B)) + C
    y[i] ~ dpois(rate[i])
  }
  
  for(k in 1:nsites_yr){
    N[k] ~ dbin(psi[k], M[k])
    psi[k] ~ dunif(0, 1)
  }
  
  # loop over sites and seasons
  for(s in 1:nSamplingOccasions) {
    ## likelihood terms
    n[s] ~ dbin(pbar[s], N[ntoN[s]])
    
    ## integrated detection probability
    pbar[s] <- sqrt_2_pi * sigma[s] * (pnorm(B, mean = 0, sd = sigma[s]) - 0.5) / B
    
    ## priors
    lsigma[s] <- alpha0 + (alphaS * season_ind[s]) + alpha1[1] * VegTypeM[sitea_ind[s], 1] + 
      alpha1[2] * VegTypeM[sitea_ind[s], 2] + alpha1[3] * VegTypeM[sitea_ind[s], 3] +
      alpha1[4] * VegTypeM[sitea_ind[s], 4] + alpha1[5] * VegTypeM[sitea_ind[s], 5] +
      alpha1[6] * VegTypeM[sitea_ind[s], 6] + alpha1[7] * VegTypeM[sitea_ind[s], 7] +
      beta1[1] * season_ind[s] * VegTypeM[sitea_ind[s], 1] +
      beta1[2] * season_ind[s] * VegTypeM[sitea_ind[s], 2] +
      beta1[3] * season_ind[s] * VegTypeM[sitea_ind[s], 3] +
      beta1[4] * season_ind[s] * VegTypeM[sitea_ind[s], 4] +
      beta1[5] * season_ind[s] * VegTypeM[sitea_ind[s], 5] +
      beta1[6] * season_ind[s] * VegTypeM[sitea_ind[s], 6] +
      beta1[7] * season_ind[s] * VegTypeM[sitea_ind[s], 7]
    sigma[s] <- exp(lsigma[s])
  }
  
  alpha0 ~ dnorm(0, sd = 100)
  
  #  for(seas in 1:nSeasons) {
  alphaS ~ dnorm(0, sd = 100)
  # }
  
  for(v in 1:nVegType) {
    alpha1[v] ~ dnorm(0, sd = 100)
  }
  
  for(j in 1:nVegType) {
    beta1[j] ~ dnorm(0, sd = 100)
  }
  
  
  # Derived parameter: total population size across all sites
  #Ntotal <- sum(N[1:nsites])
})

# Bundle and summarise data set
data <- list(y = rep(0, length(d)), d = d, n = M$M - aug)

consts <- list(
  nind = nind,
  B = B, 
  site = site,
  nsites_yr = length(unique(ntoN)),
  nSamplingOccasions = nrow(M),
  ntoN = ntoN,
  occasion = occasion,
  M = M$maxMsite,
  site_inds = site_inds,
  sitea_ind = M$site,
  season_ind = M$season - 1,
  sqrt_2_pi = sqrt(2 * pi),
  VegTypeM = VegTypeM,
  nVegType = 7)


## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
inits <- function(nVegType, M, aug, model, nsites) {
  valid <- 0
  while(valid == 0) {
    inits1 <- list(
      alpha0 = rnorm(1, 0, 1),
      alpha1 = rnorm(7, 0, 1),
      alphaS = rnorm(1, 0, 1),
      beta1 = rnorm(7, 0, 1)
#      psi = runif(max(ntoN), 0, 1)
     # N = round(runif(max(ntoN), maxMsite - aug, maxMsite - 20))
    )
    ## check validity of initial values
    model$setInits(inits1)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  inits1
}
inits1 <- inits(nVegType, M$M, aug, model, nsites)

# Parameters to save
params <- c("sigma", "psi", "N", "alpha1", "alphaS", "pbar", "alpha0", "beta1")

# MCMC settings
ni <- 20000
nb <- 10000
nt <- 2
nc <- 2

# Call nimble from R
out1 <- nimbleMCMC(code = code, 
                   data = data,
                   constants = consts, 
                   inits = inits1, 
                   monitors = params,
                   nburnin = nb, 
                   niter = ni,
                   nchains = nc,
                   #summary = TRUE,
                   samplesAsCodaMCMC = TRUE)

#out1$summary
mcmcplot(out1, parms = c("alpha1", "alphaS", "beta1", "alpha0"))

out2 <- as.matrix(out1) %>%
  as.data.frame() %>%
  select(starts_with("pbar") | starts_with("N")) %>%
  as.matrix()
out3 <- list()

for (i in 1:nrow(M)) {
  # Get the appropriate index for N using the mapping
  N_idx <- M$site[i]
  
  temp <- out2[, paste0("pbar[", i, "]")]
  tempN <- out2[, paste0("N[", N_idx, "]")]
  
  # Simulate data based on the pbar and corresponding N
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
p + xlim(0, 80)
#which(p$inside == FALSE)
