library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)
library(tikzDevice)
library(ggplot2)


rm(list = ls())

source('sas_bina.R')

ninter = 22
n1 = 11
N = 4
M = 1
n.adapt = 1000
n.burn = 1000
n.iter = 5000

p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
a0 = c(0.15, 0.15, 0.15, 0.15) ## null activity level
rho0 = 0.5
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

prob = c(0.15, 0.15, 0.15, 0.15) ## true p
acti = c(0.15, 0.15, 0.15, 0.15)  ## true activity
cluster = c(1, 1, 1, 1) ## true cluster structure
mu1 = qnorm(prob) - qnorm(p0)
mu2 = qnorm(acti) - qnorm(a0)

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)

#early_stop = matrix(0, N, M)
reject_prob = reject_acti = matrix(0, N, M)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(NA, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(NA, N, M)
cluster = getCluster(N)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  #set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = as.numeric(Z[, 2] > 0)
  }
  
  prob_rec = prob_est = prob_upper_rec = prob_lower_rec = NULL 
  acti_rec = acti_est = acti_upper_rec = acti_lower_rec = NULL 
  
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    res = post_crm(outcome, ninter, group, cutoff, cutoff2, n.adapt, n.burn, n.iter)
    this_p_c0 = res$p_c0_rec
    p_c0_est = rbind(p_c0_est, as.numeric(rowMeans(this_p_c0)))
    p_c0_upper_rec = rbind(p_c0_upper_rec, apply(this_p_c0, 1, quantile, 0.975))
    p_c0_lower_rec = rbind(p_c0_lower_rec, apply(this_p_c0, 1, quantile, 0.025))
    this_p_c1 = res$p_c1_rec
    p_c1_est = rbind(p_c1_est, as.numeric(rowMeans(this_p_c1)))
    p_c1_upper_rec = rbind(p_c1_upper_rec, apply(this_p_c1, 1, quantile, 0.975))
    p_c1_lower_rec = rbind(p_c1_lower_rec, apply(this_p_c1, 1, quantile, 0.025))
    this_p_c2 = res$p_c2_rec
    p_c2_est = rbind(p_c2_est, as.numeric(rowMeans(this_p_c2)))
    p_c2_upper_rec = rbind(p_c2_upper_rec, apply(this_p_c2, 1, quantile, 0.975))
    p_c2_lower_rec = rbind(p_c2_lower_rec, apply(this_p_c2, 1, quantile, 0.025))
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



resp = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/Results/triCRM/resp_mix.csv")[, -1])
acti = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/Results/triCRM/acti_mix.csv")[, -1])

M = 50
response = as.numeric(c(resp[1, ], resp[2, ], resp[3, ], resp[4, ]))
activity = as.numeric(c(acti[1, ], acti[2, ], acti[3, ], acti[4, ]))
arm = c(rep("arm 1", M), rep("arm 2", M), rep("arm 3", M), rep("arm 4", M))
arm = factor(arm, levels = c("arm 1", "arm 2", "arm 3", "arm 4"))
rst = data.frame("response" = response, "activity" = activity, "arm" = arm)


setwd("~/Dropbox/Mayo-intern/Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
ggplot(rst, aes(x = arm, y = response, fill = arm)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
ggplot(rst, aes(x = arm, y = activity, fill = arm)) + 
  geom_boxplot(alpha = 1, width = 0.7, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 1) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



