library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)

rm(list = ls())

source('hbm.R')

ninter = 19
n1 = 10
N = 4
C = 3
M = 20

epsilon_p = 0.05
epsilon_mu = 0.5
p0 = 0.15 ## null response rate
mu0 = 3 ## null activity level
prob = c(0.15, 0.15, 0.15, 0.45) ## true p
mu1 = qnorm(prob)
mu2 = c(3, 3, 3, 5)  ## true mu
rho0 = 0.5
cluster = c(1, 1, 1, 3) ## true cluster structure

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = array(0, c(N, ninter, 2)) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0 + epsilon_p)
cutoff2 = mu0 + epsilon_mu
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
#post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
#post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(0, N, M)


for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[i, , 1] > 0)
    activity[i, ] = Z[i, , 2]
  }
  
  ## stage 1 with only sctivity
  activity_s1 = activity[, 1:n1]
  bayes_cluster = NULL
  #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  #mu_rec = mu_upper_rec = mu_lower_rec = NULL 
  #mu_rec = NULL
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    dat = list(activity = activity_s1,
               N = N,
               ninter = n1,
               group = group,
               cutoff2 = cutoff2)
    this_posterior = posterior_simu_s1(dat, C)
    
    #this_mu = this_posterior$mu2
    #mu_rec = rbind(mu_rec, rowMeans(this_mu))
    #mu_upper_rec = rbind(mu_upper_rec, apply(this_mu, 1, quantile, 0.975))
    #mu_lower_rec = rbind(mu_lower_rec, apply(this_mu, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior_s1(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  if (sum(s1_cluster[index, ] == 1) == 4) {
    post_cluster_all[, m] = c(1, 1, 1, 1)
    next
  }
  #post_cluster_all[, m] = all_cluster[index, ]
  #post_prob_all[, m] = prob_rec[index, ]
  #post_prob_upper_all[, m] = prob_upper_rec[index, ]
  #post_prob_lower_all[, m] = prob_lower_rec[index, ]
  #post_acti_all[, m] = mu_rec[index, ]
  #post_acti_upper_all[, m] = mu_upper_rec[index, ]
  #post_acti_lower_all[, m] = mu_lower_rec[index, ]
  
  ## stage 2 with both response and activity
  bayes_cluster = NULL
  arm_remain = which(s1_cluster[index, ] == 2)
  N_remain = length(arm_remain)
  response_remain = response[arm_remain, , drop = FALSE]
  activity_remain = activity[arm_remain, , drop = FALSE]
  all_cluster = permutations(n = C, r = N_remain, repeats.allowed = T)
  #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  #acti_rec = acti_upper_rec = acti_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    dat = list(response = response_remain,
               activity = activity_remain,
               N = N_remain,
               ninter = ninter,
               group = group,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    this_posterior = posterior_simu(dat, C)
    
    #this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    #prob_rec = rbind(prob_rec, rowMeans(this_prob))
    #prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    #prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    #this_acti = pnorm(0, mean = this_posterior$mu2, sd = 1, lower.tail = FALSE)
    #acti_rec = rbind(acti_rec, rowMeans(this_acti))
    #acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 0.975))
    #acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[-arm_remain, m] = rep(1, N - N_remain)
  post_cluster_all[arm_remain, m] = all_cluster[index, ]
  #post_prob_all[, m] = prob_rec[index, ]
  #post_prob_upper_all[, m] = prob_upper_rec[index, ]
  #post_prob_lower_all[, m] = prob_lower_rec[index, ]
  #post_acti_all[, m] = acti_rec[index, ]
  #post_acti_upper_all[, m] = acti_upper_rec[index, ]
  #post_acti_lower_all[, m] = acti_lower_rec[index, ]
}


## Cluster report
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3")
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


