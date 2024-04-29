# Load necessary libraries
library(dplyr)
library(ggplot2)
library(nimble)

# Set seed for reproducibility
set.seed(123)

# Define parameters
M <- 2500              # Super population size
Tr <- 10               # Number of transects
S <- 10                # Number of sections per transect
max_distance <- 50     # Maximum observation distance in meters
years <- 5             # Number of years to simulate
veg_types <- 1:8       # Vegetation types
seasons <- 1:2         # Season types (1 for dry, 2 for wet)

# New alpha parameters for the log-linear model of sigma
# Season effect (alpha0 for dry and wet)
alpha0 <- c(log(5), log(2))  # Using log to match the model's log-linear nature

# Vegetation type effects (alpha[1] to alpha[8])
# Assuming these are to be defined or calculated elsewhere
alpha <- runif(8, -0.2, 0.2)  # Example: Random effects for vegetation types on a log scale

N_sites <- S * Tr

# Initialize VegTypeA with zeros
VegTypeA <- matrix(0, nrow=N_sites, ncol=length(veg_types))

# Assign exactly one vegetation type per site
for(i in 1:N_sites) {
  # Randomly select one vegetation type for each site
  veg_type_index <- sample(length(veg_types), 1)
  VegTypeA[i, veg_type_index] <- 1
}

# Updated detection function to use the new sigma model
half_normal_detection <- function(distance, sector, season) {
  # Calculate log(sigma) using the new model
  log_sigma <- alpha0[season] + sum(alpha * VegTypeA[sector, ])
  
  # Convert log_sigma to sigma
  sigma_adjusted <- exp(log_sigma)
  
  # Calculate detection probability with the adjusted sigma
  prob <- exp(-(distance^2) / (2 * sigma_adjusted^2))
  return(prob)
}

# Simulate data function, adjusted for the revised detection model
simulate_data <- function(year) {
  # Adjust population for trend
  M_year <- round(M * 1.05^(year - 1))
  
  # Generate observations
  observations <- data.frame()
  for (t in 1:Tr) {
    for (s in 1:S) {
      # Calculate unique section number across transects
      unique_section_number <- ((t - 1) * S) + s
      
      section_pop <- round(M_year / (Tr * S))
      for (i in 1:section_pop) {
        distance <- runif(1, 0, max_distance)
        veg_type <- sample(veg_types, 1)
        season <- sample(seasons, 1)
        detection_prob <- half_normal_detection(distance, veg_type, season)
        detected <- runif(1) < detection_prob
        if (detected) {
          observations <- rbind(observations, data.frame(Year = year, Transect = t, Section = unique_section_number, Distance = distance, Vegetation = veg_type, Season = season))
        }
      }
    }
  }
  
  return(observations)
}

# Generate data for each year
all_data <- data.frame()
for (year in 1:years) {
  yearly_data <- simulate_data(year)
  all_data <- rbind(all_data, yearly_data)
}

# View the structure of the generated data
str(all_data)

# You can then analyze the simulated data, for example, by visualizing the detection probability by distance
ggplot(all_data, aes(x = Distance)) + geom_histogram(bins = 30) + facet_wrap(~Year) + ggtitle("Detected Individuals by Distance and Year")

#all_data$Transect <- as.factor(all_data$Transect)
#all_data$Section <- as.factor(all_data$Section)

### Try my model to see if it can find starting parameters!?
year1 <- all_data %>%
  filter(Year == 1) %>%
  mutate(y = 1)

## Augment some data to the sites that were visited
data_aug <- bind_rows(year1) %>%
  group_by(Season) %>%
  do({
    season_data <- .
    # Directly work with numeric values to find unique visited sites
    sites_visited <- unique(season_data$Section)
    nz <- length(sites_visited) * 20  # Augment by 20 at each site
    aug_site <- rep(sites_visited, each = 20)
    aug_dist <- rep(NA, nz)
    
    # Ensure season column is included for augmented data, matching the current group's season
    augmented_season <- rep(unique(season_data$Season), times = nz)
    
    # Combine the original and augmented data for the season
    data.frame(
      Season = c(season_data$Season, augmented_season),
      Section = c(season_data$Section, aug_site), # Keep as numeric
      Distance = c(season_data$Distance, aug_dist),
      y = c(season_data$y, rep(0, nz))  # Distinguish augmented data
    )
  }) %>%
  ungroup()

## Extract data
dd <- data_aug$Distance[data_aug$Season == 1]
dw <- data_aug$Distance[data_aug$Season == 2]
add_dw <- rep(0, times = length(dd) - length(dw))
dw <- c(dw, add_dw)
d  <- array(c(dd, dw), dim = c(2, length(dd)))

sited <- as.numeric(data_aug$Section[data_aug$Season == 1])
sitew <- as.numeric(data_aug$Section[data_aug$Season == 2])
add_sw <- rep(1, times = length(sited) - length(sitew))
sitew <- c(sitew, add_sw)
site <- array(c(sited, sitew), dim = c(2, length(sited)))
nsitesMax <- max(site)

yd <- data_aug$y[data_aug$Season == 1]
yw <- data_aug$y[data_aug$Season == 2]
add_yw <- rep(0, times = length(yd) - length(yw))
yw <- c(yw, add_yw)
y <- array(c(yd, yw), dim = c(2, length(yd)))

nind <- c(sum(data_aug$y[data_aug$Season == 1]), sum(data_aug$y[data_aug$Season == 2]))
nz <- c(1826 - nind[1], 1826 - nind[2])

code <- nimbleCode({
  
  # Prior distributions
  beta0 ~ dunif(0, 100)   # Intercept of lambda
  
  for (v in 1:nVegTypes){
    alpha[v] ~ dnorm(0, 0.1)  # Slope of log(sigma) on vegetation type
  }
  
  for (f in 1:nseason){
    
    # Intercept of log(sigma) for each season
    alpha0[f] ~ dunif(-10, 10) 
    
    # psi is a derived parameter under DA for stratified populations
    psi[f] <- sum(lambda[f, 1:nsites[f]]) / (nind[f] + nz[f])
    
    # 'Likelihood' (sort of...)
    for (i in 1:(nind[f] + nz[f])){               # i is index for individuals
      z[f, i] ~ dbern(psi[f])                    # Data augmentation variables
      d[f, i] ~ dunif(0, B)                      # distance uniformly distributed
      p[f, i] <- exp(-d[f, i] * d[f, i] / (2 * sigma[f, site[f, i]] * sigma[f, site[f, i]])) # Det. function
      mu[f, i] <- z[f, i] * p[f, i]              # 'straw man' for WinBUGS
      y[f, i] ~ dbern(mu[f, i])                  # basic Bernoulli random variable
    }
    
    for (s in 1:nsites[f]){                    # s is index for sites
      # next line not necessary, but allows to make predictions
      N[f, s] ~ dpois(lambda[f, s])              # Realized abundance at site s
      log(lambda[f, s]) <- beta0
      
      
      # Linear model for detection
      log(sigma[f, s]) <- alpha0[f] + alpha[1] * VegTypeA[s, 1] + alpha[2] * VegTypeA[s, 2] + alpha[3] * VegTypeA[s, 3] +
        alpha[4] * VegTypeA[s, 4] + alpha[5] * VegTypeA[s, 5] + alpha[6] * VegTypeA[s, 6] +
        alpha[7] * VegTypeA[s, 7] + alpha[8] * VegTypeA[s, 8]
    }
  }
  
  # Derived parameter: total population size across all site
  Ntotal_dry <- sum(z[1, 1:2133])
  Ntotal_wet <- sum(z[2, 1:2133])
  
})

# Bundle and summarize data set
data <- list(y=y, d=d)

length(unique(data_aug$Section[data_aug$Season == 1]))

consts <- list(nsites=c(100, 100), VegTypeA = VegTypeA, nVegTypes = 8, B=50, nind=nind, 
               nz=nz, site=site, nseason=2)
# Inits
zst <- y
d_init <- d
d_init[is.na(d_init)] <- runif(length(d_init[is.na(d_init)]), 0, B)
N_init <- array(ceiling(runif(nseason * nsitesMax, 1, 30)), dim = c(nseason, nsitesMax))
sigma_init <- array(1, dim = c(nseason, nsitesMax))
lambda_init <- array(exp(1), dim = c(nseason, nsitesMax))

inits <- list(beta0 = 1, 
              alpha0 = rep(1, 2), 
              alpha = rep(0, nVegTypes),
              N = N_init,
              z = zst,
              d = d_init,
              sigma = sigma_init)
#              lamda = lambda_init)

# Parameters to save
params <- c("alpha0", "alpha", "psi", "Ntotal")

# Build model
model <- nimbleModel(code, constants = consts, data = data) #, inits = inits)
model$initializeInfo()

## compile the model
cIndicatorModel <- compileNimble(model)

## configure MCMC
config <- configureMCMC(model)
config$removeSamplers("alpha")
config$addSampler(target = c("alpha"), type = 'AF_slice')

config$resetMonitors()
config$addMonitors(c("psi", "alpha", "alpha0", "beta0", "Ntotal_dry", "Ntotal_wet"))
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 80000, 
                           nburnin = 12000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, "C2019_VSD_SVA_80Kout.rds")

parms <- c("psi", "z", "alpha", "alpha0")
mcmcplot(run)#, parms = parms)

