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
prob = c(0.15, 0.15, 0.15, 0.15)  ## true activity
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
cutoff = log(p_c1 / p_c0)[1]
cutoff2 = log(p_c2 / (1 - p_c2))[1]

## true parameters for triCRM
mu1 = qnorm(prob)
mu2 = qnorm(prob)
p_c0 = pbivnorm(-mu1, -mu2, rho0)
p_c1 = pnorm(0, mean = mu1) - p_c0
p_c2 = pnorm(0, mean = mu1, lower.tail = FALSE)

#cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
#cutoff_int2 = epsilon_2
#s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
#early_stop = matrix(0, N, M)
post_p_c0_all = post_p_c0_upper_all = post_p_c0_lower_all = matrix(NA, N, M)
post_p_c1_all = post_p_c1_upper_all = post_p_c1_lower_all = matrix(NA, N, M)
post_p_c2_all = post_p_c2_upper_all = post_p_c2_lower_all = matrix(NA, N, M)
all_cluster = permutations(n = C, r = N, repeats.allowed = T)


pb = txtProgressBar(style = 3)
for (m in 1:M) {
  #set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = as.numeric(Z[, 2] > 0)
    for (j in 1:ninter) {
      if (response[i, j] == 0 & activity[i, j] == 0) {
        outcome[i, j] = 1
      } else if (response[i, j] == 0 & activity[i, j] == 1) {
        outcome[i, j] = 2
      } else {
        outcome[i, j] = 3
      }
    }
  }
  
  bayes_cluster = NULL
  p_c0_est = p_c0_upper_rec = p_c0_lower_rec = NULL 
  p_c1_est = p_c1_upper_rec = p_c1_lower_rec = NULL
  p_c2_est = p_c2_upper_rec = p_c2_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    res = post_crm(outcome, ninter, group, cutoff, cutoff2, n.adapt, n.burn, n.iter)
    p_c0_est = rbind(p_c0_est, as.numeric(rowMeans(res$p_c0_rec)))
    p_c0_upper_rec = rbind(p_c0_upper_rec, apply(p_c0_est, 1, quantile, 0.975))
    p_c0_lower_rec = rbind(p_c0_lower_rec, apply(p_c0_est, 1, quantile, 0.025))
    p_c1_est = rbind(p_c1_est, as.numeric(rowMeans(res$p_c1_rec)))
    p_c1_upper_rec = rbind(p_c1_upper_rec, apply(p_c1_est, 1, quantile, 0.975))
    p_c1_lower_rec = rbind(p_c1_lower_rec, apply(p_c1_est, 1, quantile, 0.025))
    p_c2_est = rbind(p_c2_est, as.numeric(rowMeans(res$p_c2_rec)))
    p_c2_upper_rec = rbind(p_c2_upper_rec, apply(p_c2_est, 1, quantile, 0.975))
    p_c2_lower_rec = rbind(p_c2_lower_rec, apply(p_c2_est, 1, quantile, 0.025))
    bayes_cluster = c(bayes_cluster, res$factor)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[, m] = all_cluster[index, ]
  post_p_c0_all[, m] = p_c0_est[index, ]
  post_p_c0_upper_all[, m] = p_c0_upper_rec[index, ]
  post_p_c0_lower_all[, m] = p_c0_lower_rec[index, ]
  post_p_c1_all[, m] = p_c1_est[index, ]
  post_p_c1_upper_all[, m] = p_c1_upper_rec[index, ]
  post_p_c1_lower_all[, m] = p_c1_lower_rec[index, ]
  post_p_c2_all[, m] = p_c2_est[index, ]
  post_p_c2_upper_all[, m] = p_c2_upper_rec[index, ]
  post_p_c2_lower_all[, m] = p_c2_lower_rec[index, ]
  
  setTxtProgressBar(pb, m / M)
}


## report
report = cbind(cluster,
               rowMeans(post_cluster_all == 1),
               rowMeans(post_cluster_all == 2),
               rowMeans(post_cluster_all == 3),
               p_c0,
               rowMeans(post_p_c0_all, na.rm = TRUE),
               rowMeans(post_p_c0_lower_all < p_c0 & post_p_c0_upper_all > p_c0, na.rm = TRUE),
               p_c1,
               rowMeans(post_p_c1_all, na.rm = TRUE),
               rowMeans(post_p_c1_lower_all < p_c1 & post_p_c1_upper_all > p_c1, na.rm = TRUE),
               p_c2,
               rowMeans(post_p_c2_all, na.rm = TRUE),
               rowMeans(post_p_c2_lower_all < p_c2 & post_p_c2_upper_all > p_c2, na.rm = TRUE))
report = as.data.frame(report)
colnames(report) = c("cluster", "C1", "C2", "C3", "p_c0", "p_c0_hat", "p_c0_CI", "p_c1", "p_c1_hat", "p_c1_CI", "p_c2", "p_c2_hat", "p_c2_CI")
report




