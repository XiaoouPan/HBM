library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(tikzDevice)
library(ggplot2)

rm(list = ls())

source('hbm_cont.R')

ninter = 22
N = 4
n1 = 11
C = 3
M = 100
n.adapt = 1000
n.burn = 1000
n.iter = 5000

epsilon_p = 0.15 ## buffer for the response in the final stage
epsilon_mu_seq = c(0.5, 0.6, 0.7, 0.8, 0.9) ## buffer for the activity in the final stage
epsilon_1 = 0.02
epsilon_2 = 0.05  ## buffer for the first stage
p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
mu0 = c(3, 3, 3, 3) ## null activity level
rho0 = 0.75
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

prob = c(0.15, 0.15, 0.15, 0.15) ## true p
acti = c(3, 3, 3, 3)  ## true activity

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0[1] + epsilon_p) - qnorm(p0[1])
cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
cutoff_int2 = epsilon_2
weak = matrix(0, 4, 5)
strong = matrix(0, 4, 5)
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(qnorm(prob)[i], acti[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = Z[, 2]
  }
  
  ## Estimate correlation as preliminary analysis
  cor_est = get_cor(response[, 1:n1], activity[, 1:n1], N, n1, p0, mu0, n.adapt, n.burn, n.iter)
  
  ## Interim stage with only one outcome
  bayes_cluster = NULL
  activity_s1 = activity[, 1:n1]
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    res = post_s1_acti(activity_s1, n1, group, cutoff_int2, mu0, n.adapt, n.burn, n.iter)
    bayes_cluster = c(bayes_cluster, res$factor)
  }
  index2 = which.max(bayes_cluster)

  bayes_cluster = NULL
  response_s1 = response[, 1:n1]
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    res = post_s1_resp(response_s1, n1, group, cutoff_int1, p0, n.adapt, n.burn, n.iter)
    bayes_cluster = c(bayes_cluster, res$factor)
  }
  index1 = which.max(bayes_cluster)
  index = ifelse(mean(cor_est > 0.5) > 0.9, index2, index1)

  if (sum(s1_cluster[index, ] == 1) == 4) {
    setTxtProgressBar(pb, m / M)
    next
  }
  
  arm_remain = which(s1_cluster[index, ] == 2)
  N_remain = length(arm_remain)
  response_remain = response[arm_remain, , drop = FALSE]
  activity_remain = activity[arm_remain, , drop = FALSE]
  all_cluster = permutations(n = C, r = N_remain, repeats.allowed = T)
  
  for (j in 1:5) {
    cutoff2 = epsilon_mu_seq[j]
    bayes_cluster = NULL
    weak_rec = strong_rec = NULL 
    for (k in 1:nrow(all_cluster)) {
      group = all_cluster[k, ]
      res = post(response_remain, activity_remain, ninter, group, cutoff, cutoff2, p0[arm_remain], mu0[arm_remain], n.adapt, n.burn, n.iter)
      this_prob = pnorm(0, mean = qnorm(p0[arm_remain]) + res$mu1_rec, sd = 1, lower.tail = FALSE)
      this_acti = mu0[arm_remain] + res$mu2_rec
      weak_rec = rbind(weak_rec, as.numeric(rowMeans(this_prob > p0[arm_remain] | this_acti > mu0[arm_remain]) > reject_rate))
      strong_rec = rbind(strong_rec, as.numeric(rowMeans(this_prob > p0[arm_remain] & this_acti > mu0[arm_remain]) > reject_rate))
      bayes_cluster = c(bayes_cluster, res$factor)# this the result vector of BF after iterating thru every permutation
    }
    index = which.max(bayes_cluster)
    weak[, j] = weak[, j] + weak_rec[index, ]
    strong[, j] = strong[, j] + strong_rec[index, ]
  }
  setTxtProgressBar(pb, m / M)
}

setwd("~/Dropbox/Mayo-intern/HBM_Simulation/Results/Test_s2/s2_0act")
weak = cbind(as.matrix(read.csv("weak_015.csv")[, -1]), 
             as.matrix(read.csv("weak_018.csv")[, -1]), 
             as.matrix(read.csv("weak_021.csv")[, -1]),
             as.matrix(read.csv("weak_024.csv")[, -1]),
             as.matrix(read.csv("weak_027.csv")[, -1]))
strong = cbind(as.matrix(read.csv("strong_015.csv")[, -1]), 
             as.matrix(read.csv("strong_018.csv")[, -1]), 
             as.matrix(read.csv("strong_021.csv")[, -1]),
             as.matrix(read.csv("strong_024.csv")[, -1]),
             as.matrix(read.csv("strong_027.csv")[, -1]))

## Weak
k = 2
rst1 = c(rep(0.15, 5), rep(0.18, 5), rep(0.21, 5), rep(0.24, 5), rep(0.27, 5))
rst2 = rep(c(0.5, 0.6, 0.7, 0.8, 0.9), 5)
rst3 = weak[k, ] / 100
dat = as.data.frame(cbind(rst1, rst2, rst3))
colnames(dat) = c("Response", "Activity", "Power")

setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
ggplot(dat, aes(x = Response, y = Activity, z = Power)) + geom_contour_filled() + 
  ylim(0.5, 0.9) + 
  scale_x_continuous(breaks = seq(0.15, 0.27, by = 0.03)) + 
  xlab("Buffer for response") + ylab("Buffer for activity") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

### Strong
k = 2
rst1 = c(rep(0.15, 5), rep(0.18, 5), rep(0.21, 5), rep(0.24, 5), rep(0.27, 5))
rst2 = rep(c(0.5, 0.6, 0.7, 0.8, 0.9), 5)
rst3 = strong[k, ] / 100
dat = as.data.frame(cbind(rst1, rst2, rst3))
colnames(dat) = c("Response", "Activity", "Power")

setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
ggplot(dat, aes(x = Response, y = Activity, z = Power)) + geom_contour_filled(aes(z = Power)) + 
  ylim(0.5, 0.9) + 
  scale_x_continuous(breaks = seq(0.15, 0.27, by = 0.03)) + 
  xlab("Buffer for response") + ylab("Buffer for activity") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


