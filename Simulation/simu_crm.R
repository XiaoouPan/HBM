library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)

rm(list = ls())

source('sas_bina.R')

ninter = 22
n1 = 11
N = 4
C = 3
M = 1
n.adapt = 1000
n.burn = 1000
n.iter = 5000

epsilon_p = 0.1
epsilon_a = 0.1
#epsilon_1 = 0.05
#epsilon_2 = 0.2 ## buffer for the first stage
p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
a0 = c(0.15, 0.15, 0.15, 0.15) ## null activity level
rho0 = 0.5
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

prob = c(0.15, 0.15, 0.15, 0.15) ## true p
acti = c(0.15, 0.15, 0.15, 0.15)  ## true activity
cluster = c(1, 1, 1, 1) ## true cluster structure

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
outcome = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)

## calculate corresponding cutoffs
mu1 = qnorm(p0 + epsilon_p)
mu2 = qnorm(a0 + epsilon_a)
p_c0 = pbivnorm(-mu1, -mu2, rho0)
p_c1 = pnorm(0, mean = mu1) - p_c0
p_c2 = pnorm(0, mean = mu1, lower.tail = FALSE)
cutoff = log(p_c1 / p_c0)
cutoff2 = log(p_c2 / (1 - p_c2))

## null parameters
mu1 = qnorm(p0)
mu2 = qnorm(a0)
p_c0 = pbivnorm(-mu1, -mu2, rho0)
p_c1 = pnorm(0, mean = mu1) - p_c0
p_c2 = pnorm(0, mean = mu1, lower.tail = FALSE)

#cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
#cutoff_int2 = epsilon_2
#s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
#early_stop = matrix(0, N, M)
reject_weak = reject_strong = matrix(0, N, M)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(NA, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(NA, N, M)
all_cluster = permutations(n = C, r = N, repeats.allowed = T)


pb = txtProgressBar(style = 3)
for (m in 1:M) {
  #set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = as.numeric(Z[, 2] > 0)
  }
  
  bayes_cluster = NULL
  prob_rec = prob_est = prob_upper_rec = prob_lower_rec = NULL 
  acti_rec = acti_est = acti_upper_rec = acti_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    res = post(response, activity, ninter, group, cutoff, cutoff2, n.adapt, n.burn, n.iter)
    this_prob = pnorm(0, mean = qnorm(p0) + res$mu1_rec, sd = 1, lower.tail = FALSE)
    prob_est = rbind(prob_est, as.numeric(rowMeans(this_prob)))
    prob_rec = rbind(prob_rec, as.numeric(rowMeans(this_prob > p0) > reject_rate))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 0.975))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, 0.025))
    this_acti = pnorm(0, mean = qnorm(a0) + res$mu2_rec, sd = 1, lower.tail = FALSE)
    acti_est = rbind(acti_est, as.numeric(rowMeans(this_acti)))
    acti_rec = rbind(acti_rec, as.numeric(rowMeans(this_acti > a0) > reject_rate))
    acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 0.975))
    acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, 0.025))
    bayes_cluster = c(bayes_cluster, res$factor)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[, m] = all_cluster[index, ]
  reject_prob[, m] = prob_rec[index, ]
  reject_acti[, m] = acti_rec[index, ]
  post_prob_all[, m] = prob_est[index, ]
  post_prob_upper_all[, m] = prob_upper_rec[index, ]
  post_prob_lower_all[, m] = prob_lower_rec[index, ]
  post_acti_all[, m] = acti_est[index, ]
  post_acti_upper_all[, m] = acti_upper_rec[index, ]
  post_acti_lower_all[, m] = acti_lower_rec[index, ]
  
  setTxtProgressBar(pb, m / M)
}


## report
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3),
               prob,
               rowMeans(post_prob_all, na.rm = TRUE),
               rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob, na.rm = TRUE),
               acti,
               rowMeans(post_acti_all, na.rm = TRUE),
               rowMeans(post_acti_lower_all < acti & post_acti_upper_all > acti, na.rm = TRUE),
               rowMeans(reject_prob | reject_acti, na.rm = TRUE),
               rowMeans(reject_prob & reject_acti, na.rm = TRUE))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3", "true_p", "p_hat", "p_CI", "true_a", "a_hat", "a_CI", "weak", "strong")
report




