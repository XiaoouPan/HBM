library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)

rm(list = ls())

source('hbm_cont.R')

ninter = 19
n1 = 10
N = 4
C = 3
M = 20
n.adapt = 1000
n.burn = 1000
n.iter = 5000

epsilon_p = 0.05
epsilon_mu = 0.5
epsilon_1 = 0.2 ## buffer for the first stage
p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
mu0 = c(3, 3, 3, 3) ## null activity level
rho0 = 0.5
reject_rate = 0.9 ## For hypothesis testing
prob = c(0.15, 0.15, 0.15, 0.15) ## true p
acti = c(3, 3, 3, 3)  ## true activity
mu1 = qnorm(prob) - qnorm(p0)
mu2 = acti - mu0
cluster = c(1, 1, 1, 1) ## true cluster structure

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0[1] + epsilon_p) - qnorm(p0[1])
cutoff2 = epsilon_mu
cutoff_int = epsilon_1
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
early_stop = matrix(0, N, M)
reject_prob = reject_acti = matrix(0, N, M)
#post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
#post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(0, N, M)
index = rep(0, M)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  #set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = Z[, 2]
  }
  
  ## stage 1 with only sctivity
  activity_s1 = activity[, 1:n1]
  bayes_cluster = NULL
  #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  #mu_rec = mu_upper_rec = mu_lower_rec = NULL 
  #mu_rec = NULL
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    res = post_s1(activity_s1, n1, group, cutoff_int)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index[m] = which.max(bayes_cluster)
  setTxtProgressBar(pb, m / M)
}

index

## Cluster report
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3),
               rowMeans(early_stop), 
               rowMeans(reject_prob | reject_acti),
               rowMeans(reject_prob & reject_acti))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3", "early", "weak", "strong")
report




### Details of estimates
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


