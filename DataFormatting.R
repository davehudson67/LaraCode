# Packages and data
library(unmarked)
library(tidyverse)
library(nimble)
library(abind)
library(lubridate)
library(bayesplot)

# Read in distance sampling data and remove entries with no transect ID and no distance data
tortsR <- read.csv("distancetorts.v2.csv", header = TRUE)
torts <- tortsR %>%
  drop_na(TransectID) %>%
  drop_na(DistFromTransect) %>%
  filter(TransectCode != "YMA") #remove Polymnie Island

# deal with unknown VegTypes - recording errors?? Could add this as a latent variable if unknown.
# for now just changed both to CHP for simplicity.
Unknown <- which(torts$VegType == "?")
Unknown <- c(Unknown, which(torts$VegType == "#N/A"))
tortsUK <- torts[c(Unknown), ]

## remove unknown vegtypes
torts <- anti_join(torts, tortsUK)

## select required columns
torts <- select(torts, ID, TransectID, TransectCode, DateOfTransect, SectionCode, VegType, DistFromTransect)

## adjust columns
torts$SectionCode <- as.factor(torts$SectionCode)
torts$VegType <- as.factor(torts$VegType)
torts <- droplevels(torts)

# Remove entries left over from before 2018 and with distance >50m
torts <- torts %>%
  mutate(DateOfTransect = dmy(DateOfTransect)) %>%
  filter(year(DateOfTransect) >= 2018) %>%
  filter(DistFromTransect <= 50) %>%
  mutate(season = ifelse(month(DateOfTransect) >= 5 & month(DateOfTransect) <= 10, "Dry", "Wet")) %>%
  mutate(Year = ifelse(season == "Wet" & month(DateOfTransect) %in% 1:4, year(DateOfTransect) - 1, year(DateOfTransect))) %>%
  unite("occasion", Year, season, sep = "_", remove = FALSE) %>%
  arrange(SectionCode, occasion, DateOfTransect)

# Creating a helper DataFrame for counting unique visit dates
fdates <- torts %>%
  group_by(SectionCode, occasion) %>%
  reframe(UniqueVisitDates = unique(DateOfTransect)) %>%
  ungroup()

fdates <- fdates %>%
  group_by(SectionCode, occasion) %>%
  add_tally() %>%
  arrange(UniqueVisitDates) %>%
  mutate(Count = row_number()) %>%
  ungroup() %>%
  arrange(SectionCode)

# Joining the helper DataFrame back to the original DataFrame
joined_df <- left_join(torts, fdates, by = c("SectionCode", "occasion", "DateOfTransect" = "UniqueVisitDates"))
joined_df <- arrange(joined_df, SectionCode, occasion, Count)

doubleSamp <- filter(joined_df, n == 2)

## for simplicity remove the second sample point for now.
doubleSamp <- filter(doubleSamp, Count == 2)
torts <- anti_join(torts, doubleSamp)

## check distribution of observations
torts$occasion <- as.factor(torts$occasion)
summary(torts$occasion)

## check number of sections covered during each season
head(torts)
effort <- torts %>%
  group_by(occasion) %>%
  summarise(UniqueSectionCodes = n_distinct(SectionCode))

# Create a function to augment data for each season
data_aug <- function(torts, year, n_aug) {
  
  # Filter the data for the specified year and season
  torts_year_w <- torts %>%
    filter(Year == year) %>%
    filter(season == "Wet") %>%
    droplevels()
  
  torts_year_d <- torts %>%
    filter(Year == year) %>%
    filter(season == "Dry") %>%
    droplevels()
  
  # Get the number of observations in the year
  n_ind_w <- nrow(torts_year_w)
  n_aug_w <- n_aug - n_ind_w
  
  n_ind_d <- nrow(torts_year_d)
  n_aug_d <- n_aug - n_ind_d
  
  # Create a vector of observed distances
  x_w <- c(torts_year_w$DistFromTransect, rep(NA, n_aug_w))
  x_d <- c(torts_year_d$DistFromTransect, rep(NA, n_aug_d))
  x <- cbind(x_w, x_d)
  
  # Create a vector of binary indicators for real or augmented data
  y_w <- c(rep(1, n_ind_w), rep(0, n_aug_w))
  y_d <- c(rep(1, n_ind_d), rep(0, n_aug_d))
  y <- cbind(y_w, y_d)
  
  # Create a vector of section codes
  site_w <- c(torts_year_w$SectionCode, rep(NA, n_aug_w))
  site_d <- c(torts_year_d$SectionCode, rep(NA, n_aug_d))
  site <- cbind(x_w, x_d)
  
  # Return the data in a list
  return(list(x = x, y = y, site = site))
}

# Create a function to cycle over years and create 3-dimensional arrays for x, y, and site
create_array <- function(torts, n_aug) {
  
  # Create a list to store the 3-dimensional arrays
  arrays <- list()
  
  # Cycle over years
  for (year in 2018:2022) {
    # Split the data for the specified year
    data <- data_aug(torts, year, n_aug)
    
    # Create 3-dimensional arrays for x, y, and site
    array_x <- abind(data$x, along = 2)
    array_y <- abind(data$y, along = 2)
    array_site <- abind(data$site, along = 2)
    
    # Store the 3-dimensional arrays in the list
    arrays[[year]] <- list(x = array_x, y = array_y, site = array_site)
  }
  
  # Return the list of 3-dimensional arrays
  return(arrays[2018:2022])
}

# Create the 3-dimensional arrays
arrays <- create_array(torts, 5000)

# Function to convert a list of matrices to a 3D array
convert_to_3D_array <- function(matrices_list) {
  # Bind matrices along the third dimension to create a 3D array
  return(abind(matrices_list, along = 3))
}

## Extract and convert the 'x', 'y', and 'site' variables for each year and season
x <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$x)))
y <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$y)))
site <- do.call(abind, lapply(arrays, function(year_list) convert_to_3D_array(year_list$site)))

## Adjust vegetation types
levels(torts$VegType)
summary(torts$VegType)

# changes to implement
# CRB --> EXP
# CSB --> SND
# MSC --> SMX

levels(torts$VegType) <- sapply(levels(torts$VegType), function(x) {
  if (x == "CRB") "EXP" else if (x == "CSB") "SND" else if (x == "MSC") "SMX" else x
})

## See variation in veg type for each site...
# Group by site and count unique veg types
veg_types_per_site <- torts %>%
  group_by(SectionCode) %>%
  summarise(number_of_veg_types = n_distinct(VegType)) %>%
  filter(number_of_veg_types > 1)
# So each site only ever has one veg type assigned to it

torts <- droplevels(torts)
## Additional required data
B <- 50
nsites <- length(levels(torts$SectionCode))
VegTypeM <- torts %>%
  select(SectionCode, VegType) %>%
  distinct() %>%
  as.data.frame()
VegTypeM <- model.matrix(~ -1 + VegType, data = VegTypeM)
nVegTypes <- length(levels(torts$VegType))
nseason <- length(levels(as.factor(torts$season)))
nind <- 5000 # augmented and real for each sampling occasion
nyears <- max(torts$Year) - min(torts$Year) + 1

#tidy
rm(arrays)
rm(veg_types_per_site)
rm(joined_df)
rm(fdates)


total_occasions <- n_distinct(torts$occasion)
section_codes_in_all_occasions <- torts %>%
  group_by(SectionCode) %>%
  summarise(OccasionCount = n_distinct(occasion)) %>%
  filter(OccasionCount == total_occasions)
# only 8 sections are visited in every season/year

## have a look at distances
hist(torts$DistFromTransect)
## negative exponential looks good

## save data
save.image("tortsReady.RData")
