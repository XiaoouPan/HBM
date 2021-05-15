library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(psych)
library(mvtnorm)

rm(list = ls())

posterior_simu = function (dat, C, iter = 2000) {
  thismodel = try(jags.model(file = "trial.txt", 
                             data = dat, 
                             inits = list(mu1 = rep(-0.5, dat$N),
                                          mu2 = rep(1, dat$N),
                                          mumix = c(-1, -1, 1),
                                          muprec = c(1, 1, 1),
                                          mumix2 = c(0, 5, 0),
                                          muprec2 = c(1, 1, 1)),
                             n.adapt = iter), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu1', 'mu2', 'mumix', 'muprec', 'mumix2', 'muprec2'),
                              n.iter = iter), silent = TRUE)
  if (length(names(res.bugs)) == 0) {
    return(res.bugs)
  }
  return (list(mu1 = matrix(res.bugs$mu1, nrow = dat$N),
               mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               mumix = matrix(res.bugs$mumix, nrow = C),
               muprec = matrix(res.bugs$muprec, nrow = C),
               mumix2 = matrix(res.bugs$mumix2, nrow = C),
               muprec2 = matrix(res.bugs$muprec2, nrow = C)))
}

summary_posterior = function (dataVal, mcmcVal) {
  response = dataVal$response
  activity = dataVal$activity
  all = dataVal$all
  N = dataVal$N
  rhoEst = dat$rhoEst
  ninter = dataVal$ninter
  group = dataVal$group
  
  #parm
  mu1 = mcmcVal$mu1
  mu2 = mcmcVal$mu2
  mumix = mcmcVal$mumix
  muprec = mcmcVal$muprec
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  rst = 0
  for (n in 1:N) {
    for (j in 1:ninter) {
      p1 = dnorm(activity[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1[n, ] + rhoEst * (activity[n, j] - mu2[n, ]), sd = 1, log.p = TRUE)
      if (response[n, j] == 1) {
        p2 = 1 - p2
      }
      p3 = dnorm(mu1[n, ], mean = mumix[group[n], ], sd = 1 / sqrt(muprec[group[n], ]), log = T)
      p4 = dnorm(mu2[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
      rst = rst + mean(p1, na.rm = T) + mean(p2, na.rm = T) + mean(p3, na.rm = T) + mean(p4, na.rm = T)
    }
  }
  return (rst)
}



ninter = 20
N = 4
C = 3
M = 10

p0 = 0.3
mu0 = 3
prob = c(0.15, 0.15, 0.15, 0.45) ## true p
mu1 = qnorm(prob)
mu2 = c(1, 5, 5, 5) ## true mu
rho = 0.5

all = rep(ninter, N) # number of subjects in each subgroup during first stage
response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = array(0, c(N, ninter, 2)) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho, rho, 1), 2, 2)
rhoEst = rep(0, N)
cutoff = qnorm(p0)
cutoff2 = mu0
all_cluster = permutations(n = C, r = N, repeats.allowed = T)
post_prob_all = post_prob_upper_all = post_prob_lower_all = post_acti_all = post_cluster_all = matrix(0, N, M)

for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[i, , 1] > 0)
    activity[i, ] = Z[i, , 2]
    if (length(unique(response[i, ])) == 1) {
      rhoEst[i] = 0
    } else {
      rhoEst[i] = as.numeric(cor.test(activity[i, ], response[i, ])$estimate)
    }
  }
  bayes_cluster = NULL
  prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  mu_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    dat = list(response = response,
               activity = activity,
               rhoEst = rhoEst,
               N = N,
               ninter = ninter,
               group = group,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    this_posterior = posterior_simu(dat, C)
    mcmcVal = list(mu1 = this_posterior$mu1,
                   mu2 = this_posterior$mu2,
                   mumix = this_posterior$mumix,
                   muprec = this_posterior$muprec,
                   mumix2 = this_posterior$mumix2,
                   muprec2 = this_posterior$muprec2)
    
    this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    prob_rec = rbind(prob_rec, rowMeans(this_prob))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    mu_rec = rbind(mu_rec, rowMeans(this_posterior$mu2))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior(dat, mcmcVal)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[, m] = all_cluster[index, ]
  post_prob_all[, m] = prob_rec[index, ]
  post_prob_upper_all[, m] = prob_upper_rec[index, ]
  post_prob_lower_all[, m] = prob_lower_rec[index, ]
  post_acti_all[, m] = mu_rec[index, ]
}

rowMeans(post_cluster_all == c(1, 2, 2, 3))
rowMeans(post_prob_all)
rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob)
rowMeans(post_acti_all)





#cbind(post_cluster_all, prob, post_prob_all, post_prob_lower_all, post_prob_upper_all, mu2, post_acti_all)
