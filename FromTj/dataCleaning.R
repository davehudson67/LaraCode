# Packages and data
library(tidyverse)

## load formatted data
#load("tortsReady.RData")
torts <- readRDS("TortsReady.rds")
head(torts)

torts2019 <- torts %>%
  filter(Year == 2019) %>%
  droplevels() %>%
  mutate(SectionCode = as.factor(as.numeric(SectionCode)))

data_aug <- torts2019 %>%
  group_by(season) %>%
  do({
    season_data <- .
    sites_visited <- levels(season_data$SectionCode)[table(season_data$SectionCode) > 0]
    nz <- length(sites_visited) * 20  # Augment by 20 at each site
    aug_site <- rep(sites_visited, each = 20)
    aug_site_codes <- factor(aug_site, levels = levels(season_data$SectionCode))
    aug_dist <- rep(NA, nz)
    
    # Ensure season column is included for augmented data, matching the current group's season
    augmented_season <- rep(unique(season_data$season), times = nz)
    
    # Combine the original and augmented data for the season, now including season column
    data.frame(
      season = c(season_data$season, augmented_season),
      SectionCode = factor(c(season_data$SectionCode, aug_site_codes), levels = levels(torts2019$SectionCode)),
      DistFromTransect = c(season_data$DistFromTransect, aug_dist),
      y = c(rep(1, nrow(season_data)), rep(0, nz))  # Distinguish augmented data
    )
  }) %>%
  ungroup()

## Extract data
d <- data_aug$DistFromTransect
site <- as.numeric(data_aug$SectionCode)
nsites = max(site, na.rm = TRUE) 
y <- data_aug$y
nind <- sum(y)
nz <- length(y) - nind
season <- as.numeric(as.factor(data_aug$season))
B <- 50 #Max observation distance
nseason <- 2

### Create VegType Matrix
#VegTypeM <- torts %>%
#  select(SectionCode, VegType) %>%
#  distinct() %>%
#  as.data.frame()
#VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
VegTypeM <- readRDS("VegTypeM.rds")
nVegTypes <- length(levels(torts$VegType))

## save data out
save(d, site, nsites, y, nind, nz, season, B, nseason, VegTypeM, nVegTypes, file = "data.RData")
