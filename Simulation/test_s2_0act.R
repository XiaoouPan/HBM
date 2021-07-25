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
C = 3
M = 50
n.adapt = 1000
n.burn = 1000
n.iter = 5000

epsilon_p = 0.15 ## buffer for the response in the final stage
epsilon_mu_seq = c(0.5, 0.6, 0.7, 0.8, 0.9) ## buffer for the activity in the final stage
p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
mu0 = c(3, 3, 3, 3) ## null activity level
rho0 = 0.75
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

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
weak = matrix(0, 4, 5)
strong = matrix(0, 4, 5)
all_cluster = permutations(n = C, r = N, repeats.allowed = T)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = Z[, 2]
  }
  
  for (j in 1:5) {
    epsilon_mu = epsilon_mu_seq[j]
    cutoff2 = epsilon_mu
    bayes_cluster = NULL
    prob_rec = acti_rec = NULL 
    for (k in 1:nrow(all_cluster)) {
      group = all_cluster[k, ]
      res = post(response, activity, ninter, group, cutoff, cutoff2, n.adapt, n.burn, n.iter)
      this_prob = pnorm(0, mean = qnorm(p0) + res$mu1_rec, sd = 1, lower.tail = FALSE)
      prob_rec = rbind(prob_rec, as.numeric(rowMeans(this_prob > p0) > reject_rate))
      this_acti = mu0 + res$mu2_rec
      acti_rec = rbind(acti_rec, as.numeric(rowMeans(this_acti > mu0) > reject_rate))
      bayes_cluster = c(bayes_cluster, res$factor)# this the result vector of BF after iterating thru every permutation
    }
    index = which.max(bayes_cluster)
    reject_prob = prob_rec[index, ]
    reject_acti = acti_rec[index, ]
    weak[, j] = weak[, j] + (reject_prob | reject_acti)
    strong[, j] = strong[, j] + (reject_prob & reject_acti)
  }
  setTxtProgressBar(pb, m / M)
}

setwd("~/Dropbox/Mayo-intern/HBM_Simulation/Results/Test_s2/4act")
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
k = 1
rst1 = c(rep(0.15, 5), rep(0.18, 5), rep(0.21, 5), rep(0.24, 5), rep(0.27, 5))
rst2 = rep(c(0.5, 0.6, 0.7, 0.8, 0.9), 5)
rst3 = weak[k, ] / 50
dat = as.data.frame(cbind(rst1, rst2, rst3))
colnames(dat) = c("Response", "Activity", "Power")

setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = Response, y = Activity, z = Power)) + geom_contour_filled() + 
  xlim(0.15, 0.27) + ylim(0.5, 0.9) + 
  scale_x_continuous(breaks = seq(0.15, 0.27, by = 0.03)) + 
  xlab("Buffer for response") + ylab("Buffer for activity") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

### Strong
k = 3
rst1 = c(rep(0.15, 5), rep(0.18, 5), rep(0.21, 5), rep(0.24, 5), rep(0.27, 5))
rst2 = rep(c(0.5, 0.6, 0.7, 0.8, 0.9), 5)
rst3 = strong[k, ] / 50
dat = as.data.frame(cbind(rst1, rst2, rst3))
colnames(dat) = c("Response", "Activity", "Power")

setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(dat, aes(x = Response, y = Activity, z = Power)) + geom_contour_filled(aes(z = Power)) + 
  xlim(0.15, 0.27) + ylim(0.5, 0.9) + 
  scale_x_continuous(breaks = seq(0.15, 0.27, by = 0.03)) + 
  xlab("Buffer for response") + ylab("Buffer for activity") +
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


