# x = distance
# gs = observed group size

## covariates
# t = year
# d = distance from shore
# h = altitude
# c = cue


model{
  
  # Likelihood:
  
  for(i in 1:M){M = superpopulation size
  for(j in 1:J){J = years
  
  w[i,j] ∼ dbern(psi[i, j])
  
  logit(psi[i,j]) < - alfa0[t[i,j]] + alfa1*d[i,j] + alfa2*d[i,j]*d[i,j]
  
  d[i,j] ∼ dnorm(0, 0.7)
  x[i,j] ∼ dunif(0,1) # distances
  
  np[i,j] < - -((x[i,j] * x[i,j]) / (2 * sigma2[i,j]))
  #sigma2[i,j] < - sigma[i,j] * sigma[i,j]
  
  log(sigma[i,j]) < - beta0[h[i,j] * c[i,j]]
  
  c[i,j] ∼ dcat(x[])
  p[i,j] < - exp(np[i,j])
  
  mu[i,j] < - w[i,j] * p[i,j]
  
  y[i,j] ∼ dbern(mu[i, j]) # detection yes/no
  gs[i,j] ∼ dpois(g[j])
  wgs[i,j] < - w[i,j] * (gs[i,j] + 1)
  
  } j
  
  } i
  N2008 < - sum(wgs[1:M,1])
  N2011 < - sum(wgs[1:M,2])
  N2015 < - sum(wgs[1:M,3])
  
}
#