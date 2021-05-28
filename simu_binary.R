library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)

rm(list = ls())

source('hbm.R')

ninter = 19
n1 = 10
N = 4
C = 3
M = 50

epsilon_p = 0.1
p0 = c(0.15, 0.4)
prob1 = c(0.15, 0.15, 0.15, 0.45) ## true p1
mu1 = qnorm(prob1)
prob2 = c(0.4, 0.4, 0.4, 0.7) ## true p2
mu2 = qnorm(prob2)
rho0 = 0.5
cluster = c(1, 1, 1, 3)

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = array(0, dim = c(N, ninter, 2))
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0 + epsilon_p)
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
#post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
#post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(0, N, M)

for (m in 1:M) {
  #set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[i, , 1] > 0)
    activity[i, ] = as.numeric(Z[i, , 2] > 0)
  }
  
  ## stage 1 with only sctivity
  activity_s1 = activity[, 1:n1]
  bayes_cluster = NULL
  #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  #acti_rec = acti_upper_rec = acti_lower_rec = NULL 
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    dat = list(activity = activity_s1,
               N = N,
               ninter = n1,
               group = group,
               cutoff = cutoff[2])
    this_posterior = posterior_bi_simu_s1(dat, C)
    
    #this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    #prob_rec = rbind(prob_rec, rowMeans(this_prob))
    #prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    #prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    #this_acti = pnorm(0, mean = this_posterior$mu2, sd = 1, lower.tail = FALSE)
    #acti_rec = rbind(acti_rec, rowMeans(this_acti))
    #acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 0.975))
    #acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior_bi_s1(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  if (sum(s1_cluster[index, ] == 1) == 4) {
    post_cluster_all[, m] = c(1, 1, 1, 1)
    next
  }
  
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
               cutoff = cutoff)
    this_posterior = posterior_bi_simu(dat, C)
    
    #this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    #prob_rec = rbind(prob_rec, rowMeans(this_prob))
    #prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    #prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    #this_acti = pnorm(0, mean = this_posterior$mu2, sd = 1, lower.tail = FALSE)
    #acti_rec = rbind(acti_rec, rowMeans(this_acti))
    #acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 0.975))
    #acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior_bi(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[-arm_remain, m] = rep(1, N - N_remain)
  post_cluster_all[arm_remain, m] = all_cluster[index, ]
}


## Cluster report
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3")
report







#### Details
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3),
               prob1,
               rowMeans(post_prob_all),
               rowMeans(post_prob_lower_all < prob1 & post_prob_upper_all > prob1),
               prob2,
               rowMeans(post_acti_all),
               rowMeans(post_acti_lower_all < prob2 & post_acti_upper_all > prob2))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3", "true_p1", "p1_hat", "p1_CI", "true_p2", "p2_hat", "p2_CI")
report


