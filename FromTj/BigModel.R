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
  #filter(Year == 2019) %>%
  droplevels() %>%
  mutate(SectionCode = as.factor(as.numeric(SectionCode))) %>%
  mutate(site = as.numeric(SectionCode)) %>%
  mutate(season = as.numeric(as.factor(season))) %>%
  mutate(year = as.numeric(as.factor(Year))) %>%
  mutate(Occ = as.numeric(as.factor((occasion))))

levels(torts$occasion)

## load data
#load("data.RData")

## Extract data from DA version and amend to non-DA version
site <- as.numeric(torts$SectionCode)
season <- as.numeric(as.factor(torts$season))
d <- torts$DistFromTransect
nsites <- length(unique(site))
nind <- nrow(torts)
VegType <- as.numeric(torts$VegType)
nVegType <- length(unique(torts$VegType))
year <- as.numeric(as.factor(torts$Year))
nyear <- length(unique(year))
noccasions <- length(unique(torts$Occ))
occasion <- torts$Occ
transect <- torts %>%
  mutate(TransectCode = as.numeric(as.factor(TransectCode))) %>%
  distinct(SectionCode, .keep_all = TRUE) %>%
  select(SectionCode, TransectCode)
nTransects <- length(unique(torts$TransectCode))
VegType <- torts %>%
  mutate(VegType = as.numeric(VegType)) %>%
  distinct(SectionCode, .keep_all = TRUE) %>%
  select(SectionCode, VegType)

## check all sites present
stopifnot(all((sort(unique(site)) - 1:nsites) == 0))
y <- rep(1, nind)
B <- 50 #Max observation distance
nseason <- 2
stopifnot(all((sort(unique(season)) - 1:nseason) == 0))

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
    lsigma[s] <- (alpha1 * season_ind[s]) + alpha[VegType[sitea_ind[s]]] + beta[Transect[sitea_ind[s]]] * season_ind[s]
    sigma[s] <- exp(lsigma[s])
  }
  for(s in 1:nVegType) {
    alpha[s] ~ dnorm(0, sd = 100)
  }
  for(r in 1:nTransects){
    beta[r] ~ dnorm(0, sd = tau)
  }
  alpha1 ~ dnorm(0, sd = 100)
  tau ~ dunif(0, 100)
  
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
  Transect = transect$TransectCode,
  nTransects = nTransects,
  M = M$maxMsite,
  sitea_ind = M$site,
  season_ind = M$season - 1,
  site_inds = site_inds,
  sqrt_2_pi = sqrt(2 * pi),
  VegType = VegType$VegType,
  nVegType = nVegType)

## set up model with no initial values
model <- nimbleModel(code = code, constants = consts, data = data)

## check one can generate valid initial values
nsites_yr = length(unique(ntoN))
Minit <- M %>%
  select(ntoN, maxMsite) %>%
  distinct(ntoN, .keep_all = T)
maxMsite <- Minit$maxMsite

inits <- function(nVegType, aug, model, nsites_yr, maxMsite) {
  valid <- 0
  while(valid == 0) {
    inits1 <- list(
      alpha = rnorm(nVegType, 0, sd = 1),
      alpha1 = rnorm(1, 0, 1),
      beta = rnorm(nTransects, 0, 1),
      tau = runif(1, 0, 100),
      psi = runif(nsites_yr, 0, 1),
      N = floor(runif(nsites_yr, maxMsite - aug, maxMsite - 20))
    )
    ## check validity of initial values
    model$setInits(inits1)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  inits1
}
inits1 <- inits(nVegType, aug, model, nsites_yr, maxMsite)

# Parameters to save
params <- c("sigma", "psi", "N", "alpha", "alpha1", "pbar")

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
warnings()
saveRDS(out1, "BigModelSamples.rds")
mcmcplot(out1, parms = c("N", "alpha1", "psi", "alpha", "pbar"))

out1<- readRDS("BigModelSamples.rds")

## Collect yearly estimates
DevP <- M %>%
  group_by(year) %>%
  distinct(ntoN)

T2019 <- DevP$ntoN[DevP$year == 1]
T2020 <- DevP$ntoN[DevP$year == 2]
T2021 <- DevP$ntoN[DevP$year == 3]
T2022 <- DevP$ntoN[DevP$year == 4]

## Area sampled (km2) - taken from Lara's email, max sampled sites in either season
A2018 <- 445 * 50 / 10000
A2019 <- 440 * 50 / 10000
A2020 <- 395 * 50 / 10000
A2021 <- 417 * 50 / 10000
A2022 <- 445 * 50 / 10000
TotIsland <- 155.5

## Combine MCMC output
outM <- as.data.frame(as.matrix(out1))

## Select only N
N <-  outM %>%
  select(grep("N\\[", names(outM)))

## Sum N for each year
N_2019 <- rowSums(N[, T2019])
N_2020 <- rowSums(N[, T2020])
N_2021 <- rowSums(N[, T2021])
N_2022 <- rowSums(N[, T2022])

D2019 <- N_2019 / A2019 * TotIsland
D2020 <- N_2020 / A2020 * TotIsland
D2021 <- N_2021 / A2021 * TotIsland
D2022 <- N_2022 / A2022 * TotIsland

hist(D2019)
hist(D2020)
hist(D2021)
hist(D2022)

## Posterior predictive checks
## check model fits using posterior predictive distributions
out2 <- as.matrix(out1) %>%
  as.data.frame() %>%
  select(starts_with("pbar") | starts_with("N")) %>%
  as.matrix()

out3 <- list()

for (i in 1:nrow(M)) {
  # Get the appropriate index for N using the mapping
  N_idx <- ntoN[i]
  
  temp <- out2[, paste0("pbar[", i, "]")]
  tempN <- out2[, paste0("N[", N_idx, "]")]  # Use mapped index here
  
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
  mutate(occasion = M$occasion) %>%
  mutate(site = M$site[site])
rownames(out3) <- NULL

## data
Ninit <- tapply(d, paste0(site, "_", occasion), length)
n <- matrix(NA, nsites, noccasions)
Ninit_nam <- strsplit(names(Ninit), "_")
for(i in 1:length(Ninit_nam)) {
  temp <- as.numeric(Ninit_nam[[i]])
  n[temp[1], temp[2]] <- Ninit[i]
}
p <- as.data.frame(n) %>%
  set_names(c("o1", "o2", "o3", "o4", "o5", "o6", "o7", "o8")) %>%
  mutate(site = 1:n()) %>%
  pivot_longer(!site, names_to = "occasion") %>%
  mutate(occasion = as.numeric(gsub("o", "", occasion))) %>%
  filter(!is.na(value)) %>%
  inner_join(out3, by = c("site", "occasion")) %>%
  mutate(inside = value > LCI & value < UCI) %>%
  ggplot() +
  geom_point(aes(x = site, y = mean)) +
  geom_errorbar(aes(x = site, ymin = LCI, ymax = UCI)) +
  geom_point(aes(x = site, y = value, colour = inside)) +
  facet_wrap(~occasion, ncol = 1)
p + xlim(0, 100)

summary(p$inside)
1517/139
139/1517

which(p$inside == FALSE & (p$occasion == 1 | p$occasion == 2))
which(p$inside == FALSE & (p$occasion == 3 | p$occasion == 4))
which(p$inside == FALSE & (p$occasion == 5 | p$occasion == 6))
which(p$inside == FALSE & (p$occasion == 7 | p$occasion == 8))
