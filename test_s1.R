library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)

rm(list = ls())

source('hbm_cont.R')

ninter = 22
n1 = 11
N = 4
C = 3
M = 50
n.adapt = 1000
n.burn = 1000
n.iter = 5000
alpha = 0.026
reject_rate = 1 - alpha

p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
mu0 = c(3, 3, 3, 3) ## null activity level
rho0 = 0.75
prob = c(0.15, 0.15, 0.15, 0.15) ## true p
acti = c(3, 3, 3, 3)  ## true activity
mu1 = qnorm(prob) - qnorm(p0)
mu2 = acti - mu0
cluster = c(1, 1, 1, 1) ## true cluster structure

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
reject_prob = reject_acti = matrix(0, N, M)
post_cluster_all = matrix(0, N, M)
#post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(0, N, M)
#post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(0, N, M)


######## activity
epsilon_seq = seq(0, 0.9, by = 0.05)
l = length(epsilon_seq)
early1 = early2 = matrix(0, N, l)
pb = txtProgressBar(style = 3)
for (j in 1:l) {
  epsilon_2 = epsilon_seq[j]
  #cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
  cutoff_int2 = epsilon_2
  for (m in 1:M) {
    #set.seed(m)
    ## Data generation
    for (i in 1:N) {
      Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
      response[i, ] = as.numeric(Z[, 1] > 0)
      activity[i, ] = Z[, 2]
    }
    
    activity_s1 = activity[, 1:n1]
    bayes_cluster = NULL
    #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
    #mu_rec = mu_upper_rec = mu_lower_rec = NULL 
    mu_rec = NULL
    for (i in 1:nrow(s1_cluster)) {
      group = s1_cluster[i, ]
      res = post_s1_acti(activity_s1, n1, group, cutoff_int2)
      bayes_cluster = c(bayes_cluster, res$factor)
      this_acti = res$mu2_rec
      mu_rec = rbind(mu_rec, as.numeric(rowMeans(this_acti > 0) > reject_rate))
    }
    index = which.max(bayes_cluster)
    post_cluster_all[, m] = s1_cluster[index, ]
    reject_acti[, m] = mu_rec[index, ]
    
    setTxtProgressBar(pb, ((j - 1) * M +  m) / (l * M))
  }
  early1[, j] = rowMeans(post_cluster_all == 1)
  early2[, j] = rowMeans(reject_acti)
}

early1
early2



