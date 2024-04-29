##### PACKAGES AND INTRO #####
library(unmarked)
library(R2WinBUGS)
library(jagsUI)
library(tidyverse)
library(nimble)
library(mcmcplots)
library(lubridate)

# Read in distance sampling data and remove entries with no transect ID and no distance data
torts <- read.csv("distancetorts.csv", header = TRUE) %>%
  drop_na(TransectID) %>%
  drop_na(DistFromTransect)

# Subset the data to remove entries left over from before 2018
#torts <- subset(torts, format.Date(DateOfTransect, "%Y") != c("1998", "2017"))

# Subset the data to remove entries with distance >50m
#torts <- subset(torts, DistFromTransect <= 50)

## or using dplyr and lubridate
torts <- torts %>%
  mutate(DateOfTransect = dmy(torts$DateOfTransect)) %>%
  filter(year(DateOfTransect) > 2018) %>%
  filter(DistFromTransect <= 50)

##### MAXIMUM LIKELIHOOD #####

# Define negative exponential detection function
g <- function(x, sig) exp(-x/sig)

# Full likelihood
Lfull <- function(parm) {
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  N <- length(torts$DistFromTransect)+ n0
  pbar <- integrate(g, 0, 100, sig=sigma)$value/100
  -1*( lgamma(N+1) - lgamma(n0+1) + sum(
    log(g(torts$DistFromTransect,sig=sigma)/100)) + n0*log(1-pbar) )
}

# Call optim to maximize full likelihood
optim(c(log(30), log(4)), Lfull, hessian=TRUE)

##### CONVERT ESTIMATES OF SIMPLE DENSITY - NO EXTRA PREDICTORS #####

# To convert estimates of density - divide the estimates of N by the area of the transect
n <- length(torts$DistFromTransect)

length(unique(torts$SectionCode))

# Transect width = 50m either side = 0.1km total
# 440 transect sections (50m intervals along transect): total distance = 22km, total area = 22*0.1 = 2.2km^2
n0hat <- exp(11.32)
Nhat.full <- n + n0hat
Dhat.full <- Nhat.full/(22*0.1)
