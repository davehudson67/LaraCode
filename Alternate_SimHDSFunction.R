simHDSOPEN <- function (type = c("line", "point"), nsites = 100, mean.lam = 2, 
          beta.lam = 0, mean.sig = 1, beta.sig = 0, B = 50, discard0 = TRUE, 
          nreps = 2, phi = 0.7, nyears = 5, beta.trend = 0) 
{
  nsites <- round(nsites[1])
  #stopifNegative(mean.lam, allowZero = FALSE)
  #stopifNegative(B, allowZero = FALSE)
  type <- match.arg(type)
  parmvec <- c(mean.lam, beta.lam, mean.sig, beta.sig, phi, 
               beta.trend)
  names(parmvec) <- c("mean.lam", "beta.lam", "mean.sig", "beta.sig", 
                      "phi", "beta.trend")
  habitat <- rnorm(nsites)
  M <- lambda <- matrix(NA, nrow = nsites, ncol = nyears)
  Na <- wind <- array(NA, dim = c(nsites, nreps, nyears))
  Na.real <- array(0, dim = c(nsites, nreps, nyears))
  for (i in 1:nyears) {
    lambda[, i] <- exp(log(mean.lam) + beta.lam * habitat + 
                         beta.trend * (i - nyears/2))
    M[, i] <- rpois(nsites, lambda[, i])
    Na[, , i] <- matrix(rbinom(nsites * nreps, M[, i], phi), 
                        nrow = nsites, byrow = FALSE)
    wind[, , i] <- runif(nsites * nreps, -2, 2)
  }
  sigma <- exp(log(mean.sig) + beta.sig * wind)
  outlist <- list()
  for (yr in 1:nyears) {
    list.yr <- list()
    for (rep in 1:nreps) {
      data <- NULL
      for (i in 1:nsites) {
        if (Na[i, rep, yr] == 0) {
          data <- rbind(data, c(i, NA, NA, NA, NA))
          next
        }
        if (type == "line") {
          d <- runif(Na[i, rep, yr], 0, B)
          Na.real[i, rep, yr] <- sum(d <= B)
          p <- exp(-d * d/(2 * (sigma[i, rep, yr]^2)))
          y <- rbinom(Na[i, rep, yr], 1, p)
          u1 <- u2 <- rep(NA, Na[i, rep, yr])
          d <- d[y == 1]
          u1 <- u1[y == 1]
          u2 <- u2[y == 1]
          y <- y[y == 1]
        }
        if (type == "point") {
          angle <- runif(Na[i, rep, yr], 0, 2 * pi)
          dd <- B * sqrt(runif(Na[i, rep, yr], 0, 1))
          u1 <- dd * cos(angle) + (B)
          u2 <- dd * sin(angle) + (B)
          d <- sqrt((u1 - B)^2 + (u2 - B)^2)
          Na.real[i, rep, yr] <- sum(d <= B)
          p <- exp(-d * d/(2 * (sigma[i, rep, yr]^2)))
          pp <- ifelse(d < B, 1, 0) * p
          y <- rbinom(Na[i, rep, yr], 1, pp)
          u1 <- u1[y == 1]
          u2 <- u2[y == 1]
          d <- d[y == 1]
          y <- y[y == 1]
        }
        if (sum(y) > 0) {
          data <- rbind(data, cbind(rep(i, sum(y)), y, 
                                    u1, u2, d))
        }
        else {
          data <- rbind(data, c(i, NA, NA, NA, NA))
        }
      }
      colnames(data) <- c("site", "y", "u1", "u2", "d")
      if (discard0) 
        data <- data[!is.na(data[, 2]), ]
      list.yr[[rep]] <- data
    }
    outlist[[yr]] <- list.yr
  }
  list(data = outlist, B = B, nsites = nsites, habitat = habitat, 
       wind = wind, M.true = M, K = nreps, nyears = nyears, 
       Na = Na, Na.real = Na.real, mean.lam = mean.lam, beta.lam = beta.lam, 
       mean.sig = mean.sig, beta.sig = beta.sig, phi = phi, 
       beta.trend = beta.trend, parms = parmvec)
}
