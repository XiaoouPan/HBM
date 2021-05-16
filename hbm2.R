library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)

rm(list = ls())

posterior_simu = function (dat, C, iter = 4000) {
  thismodel = try(jags.model(file = "trial.txt", 
                             data = dat, 
                             inits = list(mu1 = rep(-0.5, dat$N),
                                          mu2 = rep(1, dat$N),
                                          rho = rep(0.5, dat$N),
                                          mumix = c(-2, -2, 0),
                                          muprec = c(1, 1, 1),
                                          mumix2 = c(1, 5, 1),
                                          muprec2 = c(1, 1, 1)),
                             n.adapt = iter), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu1', 'mu2', 'rho', 'mumix', 'muprec', 'mumix2', 'muprec2'),
                              n.iter = iter), silent = TRUE)
  return (list(mu1 = matrix(res.bugs$mu1, nrow = dat$N),
               mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               rho = matrix(res.bugs$rho, nrow = dat$N),
               mumix = matrix(res.bugs$mumix, nrow = C),
               muprec = matrix(res.bugs$muprec, nrow = C),
               mumix2 = matrix(res.bugs$mumix2, nrow = C),
               muprec2 = matrix(res.bugs$muprec2, nrow = C)))
}

summary_posterior = function (dataVal, mcmcVal) {
  response = dataVal$response
  activity = dataVal$activity
  N = dataVal$N
  ninter = dataVal$ninter
  group = dataVal$group
  
  #parm
  mu1 = mcmcVal$mu1
  mu2 = mcmcVal$mu2
  rho = mcmcVal$rho
  mumix = mcmcVal$mumix
  muprec = mcmcVal$muprec
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  rst = 0
  for (n in 1:N) {
    p3 = dnorm(mu1[n, ], mean = mumix[group[n], ], sd = 1 / sqrt(muprec[group[n], ]), log = T)
    p4 = dnorm(mu2[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    rst = rst + ninter * (mean(p3, na.rm = T) + mean(p4, na.rm = T))
    for (j in 1:ninter) {
      p1 = dnorm(activity[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity[n, j] - mu2[n, ]), sd = 1, log.p = TRUE)
      if (response[n, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(p2, na.rm = T)
    }
  }
  return (rst)
}



ninter = 20
N = 4
C = 3
M = 20

p0 = 0.15
mu0 = 3
prob = c(0.15, 0.15, 0.15, 0.15) ## true p
mu1 = qnorm(prob)
mu2 = c(3, 3, 3, 3) ## true mu
rho0 = 0.5
cluster = c(1, 1, 1, 1)

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = array(0, c(N, ninter, 2)) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0)
cutoff2 = mu0
all_cluster = permutations(n = C, r = N, repeats.allowed = T)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = post_cluster_all = matrix(0, N, M)

for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[i, , 1] > 0)
    activity[i, ] = Z[i, , 2]
  }
  bayes_cluster = NULL
  prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  mu_rec = mu_upper_rec = mu_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    dat = list(response = response,
               activity = activity,
               N = N,
               ninter = ninter,
               group = group,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    this_posterior = posterior_simu(dat, C)
    
    this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    prob_rec = rbind(prob_rec, rowMeans(this_prob))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    this_mu = this_posterior$mu2
    mu_rec = rbind(mu_rec, rowMeans(this_mu))
    mu_upper_rec = rbind(mu_upper_rec, apply(this_mu, 1, quantile, 0.975))
    mu_lower_rec = rbind(mu_lower_rec, apply(this_mu, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[, m] = all_cluster[index, ]
  post_prob_all[, m] = prob_rec[index, ]
  post_prob_upper_all[, m] = prob_upper_rec[index, ]
  post_prob_lower_all[, m] = prob_lower_rec[index, ]
  post_acti_all[, m] = mu_rec[index, ]
  post_acti_upper_all[, m] = mu_upper_rec[index, ]
  post_acti_lower_all[, m] = mu_lower_rec[index, ]
}


report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3),
               prob,
               rowMeans(post_prob_all),
               rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob),
               mu2,
               rowMeans(post_acti_all),
               rowMeans(post_acti_lower_all < mu2 & post_acti_upper_all > mu2))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3", "true_p", "p_hat", "p_CI", "true_mu", "mu_hat", "mu_CI")
report


