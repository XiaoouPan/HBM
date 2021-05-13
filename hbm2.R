library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(psych)
library(mvtnorm)

rm(list = ls())

## generate data
geneData = function(N, C, p0, mu0, mu1, mu2, rho, ninter) {
  # N: number of subgroups
  # p0, mu0: scalars, null values in the null hypothesis, may be changed to vectors later
  # mu1, mu2: N-vectors
  all = rep(ninter, N) # number of subjects in each subgroup during first stage
  response = activity = matrix(0, N, ninter)
  Sigma = matrix(c(1, rho, rho, 1), 2, 2)
  for (i in 1:N) {
    dat = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(dat[, 1] > 0)
    activity[i, ] = dat[, 2]
  }
  cutoff = qnorm(p0)
  cutoff2 = mu0
  return(
    list(
      response = response,
      activity = activity,
      all = all,
      N = N,
      C = C,
      cutoff = cutoff,
      cutoff2 = cutoff2
    )
  )
}

posterior_simu = function (dat, iter = 1000) {
  thismodel = jags.model(file = "basket.txt", 
                         data = dat, 
                         inits = list(mu = cbind(rep(0, dat$N), rep(1, dat$N)),
                                      mumix = c(-1, -1, 1),
                                      muprec = c(1, 1, 1),
                                      mumix2 = c(0, 5, 0),
                                      muprec2 = c(1, 1, 1)),
                         n.adapt = iter)
  res.bugs = jags.samples(thismodel, 
                          variable.names = c('mu', 'mumix', 'muprec', 'mumix2', 'muprec2'),
                          n.iter = iter)
  if (length(names(res.bugs)) == 0) {
    return(res.bugs)
  }
  return (list(mu1 = matrix(res.bugs$mu[, 1], nrow = dat$N),
               mu2 = matrix(res.bugs$mu[, 2], nrow = dat$N),
               mumix = matrix(res.bugs$mumix, nrow = 3),
               muprec = matrix(res.bugs$muprec, nrow = 3),
               mumix2 = matrix(res.bugs$mumix2, nrow = 3),
               muprec2 = matrix(res.bugs$muprec2, nrow = 3)))
}

summary_posterior = function (dataVal, mcmcVal) {
  response = dataVal$response
  activity = dataVal$activity
  all = dataVal$all
  N = dataVal$N
  C = dataVal$C
  group = dataVal$group
  
  #parm
  mu1 = mcmcVal$mu1
  mu2 = mcmcVal$mu2
  mumix = mcmcVal$mumix
  muprec = mcmcVal$muprec
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  res1 = res2 = res3 = 0
  for (n in 1:N) {
    p1 = dbinom(x = response[n], size = all[n], prob = pnorm(0, mu1[n, ]), log = T)
    p2 = 
    p2 = dnorm(x = theta[n, ], mean = mumix[group[n], ], sd = 1 / sqrt(muprec[group[n], ]), log = T)
    
    p3 = dnorm(x = mu[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    
    res1 <- res1 + mean(p1, na.rm = T)
    res2 <- res2 + mean(p2, na.rm = T)
    res3 = res3 + mean(p3, na.rm = T)
    
  }
  return (res1 + res2 + res3)
}



ninter = 10
N = 4
C = 3
M = 20
post_prob_all = post_acti_all = matrix(0, N, M)
all_cluster <- permutations(n = C, r = N, repeats.allowed = T)

p0 = 0.3
mu0 = 3
prob = c(0.45, 0.45, 0.45, 0.45)
mu1 = qnorm(prob)
mu2 = c(1, 1, 5, 5)
rho = 0.5
trial = geneData(N, C, p0, mu0, mu1, mu2, rho, ninter)
true_p = prob
true_mu = mu2






n = 100
mu = c(0, 3)
Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
dat = mvrnorm(n, mu, Sigma)
X = as.numeric(dat[, 1] > 0)
Y = dat[, 2]
biserial(Y, X)
