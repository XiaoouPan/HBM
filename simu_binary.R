library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)

rm(list = ls())

source('hbm.R')

ninter = 20
N = 4
C = 3
M = 1

epsilon_p = 0.05
p0 = c(0.15, 0.4)
prob1 = c(0.15, 0.15, 0.15, 0.15) ## true p1
mu1 = qnorm(prob1)
prob2 = c(0.4, 0.4, 0.4, 0.4) ## true p2
mu2 = qnorm(prob2)
rho0 = 0.5
cluster = c(1, 1, 1, 1)

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = array(0, dim = c(N, ninter, 2))
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0 + epsilon_p)
all_cluster = permutations(n = C, r = N, repeats.allowed = T)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = post_cluster_all = matrix(0, N, M)

for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[i, , 1] > 0)
    activity[i, ] = as.numeric(Z[i, , 2] > 0)
  }
  bayes_cluster = NULL
  prob_rec = prob_upper_rec = prob_lower_rec = NULL 
  acti_rec = acti_upper_rec = acti_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    dat = list(response = response,
               activity = activity,
               N = N,
               ninter = ninter,
               group = group,
               cutoff = cutoff)
    this_posterior = posterior_bi_simu(dat, C)
    
    this_prob = pnorm(0, mean = this_posterior$mu1, sd = 1, lower.tail = FALSE)
    prob_rec = rbind(prob_rec, rowMeans(this_prob))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    this_acti = pnorm(0, mean = this_posterior$mu2, sd = 1, lower.tail = FALSE)
    acti_rec = rbind(acti_rec, rowMeans(this_acti))
    acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 0.975))
    acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, 0.025))
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res = summary_posterior_bi(dat, this_posterior)
    bayes_cluster = c(bayes_cluster, res)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[, m] = all_cluster[index, ]
  post_prob_all[, m] = prob_rec[index, ]
  post_prob_upper_all[, m] = prob_upper_rec[index, ]
  post_prob_lower_all[, m] = prob_lower_rec[index, ]
  post_acti_all[, m] = acti_rec[index, ]
  post_acti_upper_all[, m] = acti_upper_rec[index, ]
  post_acti_lower_all[, m] = acti_lower_rec[index, ]
}


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


