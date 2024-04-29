library(AHMbook)
# Obtain a data set
set.seed(1236)

## set some parameter
B <- 50
nsites <- 397
nyears <- 5
K <- 2

## simulate some data
#Sim <- simHDSOPEN("line", nreps = 2, nyears = 5, B = B, nsites = 397, beta.trend = 0,
#                   mean.lam = 1, beta.lam = 0, mean.sig = 1, beta.sig = 0)

Sim <- simHDSopen("line", nreps = 2, nyears = 5, B = B, nsites = 397, beta.trend = 0,
                  mean.lam = 8, beta.lam = 0, mean.sig = 1, beta.sig = 0)

## take a look at true population size per year
apply(Sim$M.true, 2, sum)  

# Define distance class information
delta <- 1 # distance interval width
nD <- B %/% delta               # Number of distance classes
midpt <- seq(delta / 2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0, dim = c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- Sim$data[[yr]][[rep]]
    site <- data[, 1]
    dclass <- data[, "d"] %/% delta + 1
    ndclass <- B %/% delta
    dclass <- factor(dclass, levels=  1:ndclass)
    # ~~~~~ this cannot work ~~~~~~~~~~
    # y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
    # ~~~~~ use this instead ~~~~~~~~~~~
    ttt <- table(site, dclass)
    siteID <- as.numeric(rownames(ttt))
    y4d[siteID, 1:nD, rep, yr] <- ttt
  }
}