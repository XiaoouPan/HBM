library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)
library(tikzDevice)
library(ggplot2)
library(xtable)

rm(list = ls())

source('hbm_bina.R')

ninter = 22
n1 = 11
N = 4
C = 3
M = 500
n.adapt = 1000
n.burn = 1000
n.adapt.c3 = 1000
n.iter = 5000

epsilon_p = 0.2
epsilon_a = 0.2
epsilon_1 = 0.02
epsilon_2 = 0.02 ## buffer for the first stage
p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
a0 = c(0.15, 0.15, 0.15, 0.15) ## null activity level
rho0 = 0.75
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

prob = c(0.15, 0.15, 0.15, 0.45) ## true p
acti = c(0.15, 0.15, 0.45, 0.45)  ## true activity
mu1 = qnorm(prob) - qnorm(p0)
mu2 = qnorm(acti) - qnorm(a0)
cluster = c(1, 1, 2, 3) ## true cluster structure

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
cutoff = qnorm(p0[1] + epsilon_p) - qnorm(p0[1])
cutoff2 = qnorm(a0[1] + epsilon_a) - qnorm(a0[1])
cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
cutoff_int2 = qnorm(a0[1] + epsilon_2) - qnorm(a0[1])
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)
early_stop = matrix(0, N, M)
reject_prob = reject_acti = matrix(0, N, M)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(NA, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(NA, N, M)


pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = as.numeric(Z[, 2] > 0)
  }
  
  ## Estimate correlation as preliminary analysis
  cor_est = get_cor(response[, 1:n1], activity[, 1:n1], N, n1, n.adapt, n.burn, n.iter)
  
  ## Interim stage with only one outcome
  prob_est = prob_upper_rec = prob_lower_rec = NULL 
  acti_est = acti_upper_rec = acti_lower_rec = NULL
  
  bayes_cluster = NULL
  activity_s1 = activity[, 1:n1]
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    res = post_s1_acti(activity_s1, n1, group, cutoff_int2, n.adapt, n.burn, n.iter)
    bayes_cluster = c(bayes_cluster, res$factor)
    this_acti = pnorm(0, mean = qnorm(a0) + res$mu2_rec, sd = 1, lower.tail = FALSE)
    acti_est = rbind(acti_est, as.numeric(rowMeans(this_acti)))
    acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 1 - alpha / 2))
    acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, alpha / 2))
  }
  index2 = which.max(bayes_cluster)
  post_acti_all[, m] = acti_est[index2, ]
  post_acti_upper_all[, m] = acti_upper_rec[index2, ]
  post_acti_lower_all[, m] = acti_lower_rec[index2, ]
  
  bayes_cluster = NULL
  response_s1 = response[, 1:n1]
  for (i in 1:nrow(s1_cluster)) {
    group = s1_cluster[i, ]
    res = post_s1_resp(response_s1, n1, group, cutoff_int1, n.adapt, n.burn, n.iter)
    bayes_cluster = c(bayes_cluster, res$factor)
    this_prob = pnorm(0, mean = qnorm(p0) + res$mu1_rec, sd = 1, lower.tail = FALSE)
    prob_est = rbind(prob_est, as.numeric(rowMeans(this_prob)))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 1 - alpha / 2))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, alpha / 2))
  }
  index1 = which.max(bayes_cluster)
  post_prob_all[, m] = prob_est[index1, ]
  post_prob_upper_all[, m] = prob_upper_rec[index1, ]
  post_prob_lower_all[, m] = prob_lower_rec[index1, ]
  index = ifelse(mean(cor_est > 0.5) > 0.9, index2, index1)
  
  if (sum(s1_cluster[index, ] == 1) == 4) {
    post_cluster_all[, m] = c(1, 1, 1, 1)
    early_stop[, m] = c(1, 1, 1, 1)
    setTxtProgressBar(pb, m / M)
    next
  }
  
  ## stage 2 with both response and activity
  arm_remain = which(s1_cluster[index, ] == 2)
  N_remain = length(arm_remain)
  early_stop[-arm_remain, m] = rep(1, N - N_remain)
  response_remain = response[arm_remain, , drop = FALSE]
  activity_remain = activity[arm_remain, , drop = FALSE]
  all_cluster = permutations(n = C, r = N_remain, repeats.allowed = T)
  bayes_cluster = NULL
  prob_rec = prob_est = prob_upper_rec = prob_lower_rec = NULL 
  acti_rec = acti_est = acti_upper_rec = acti_lower_rec = NULL 
  for (i in 1:nrow(all_cluster)) {
    group = all_cluster[i, ]
    res = post(response_remain, activity_remain, ninter, group, cutoff, cutoff2, n.adapt.c3, n.burn, n.iter)
    this_prob = pnorm(0, mean = qnorm(p0[arm_remain]) + res$mu1_rec, sd = 1, lower.tail = FALSE)
    prob_est = rbind(prob_est, as.numeric(rowMeans(this_prob)))
    prob_rec = rbind(prob_rec, as.numeric(rowMeans(this_prob > p0[arm_remain]) > reject_rate))
    prob_upper_rec = rbind(prob_upper_rec, apply(this_prob, 1, quantile, 1 - alpha / 2))
    prob_lower_rec = rbind(prob_lower_rec, apply(this_prob, 1, quantile, alpha / 2))
    this_acti = pnorm(0, mean = qnorm(a0[arm_remain]) + res$mu2_rec, sd = 1, lower.tail = FALSE)
    acti_est = rbind(acti_est, as.numeric(rowMeans(this_acti)))
    acti_rec = rbind(acti_rec, as.numeric(rowMeans(this_acti > a0[arm_remain]) > reject_rate))
    acti_upper_rec = rbind(acti_upper_rec, apply(this_acti, 1, quantile, 1 - alpha / 2))
    acti_lower_rec = rbind(acti_lower_rec, apply(this_acti, 1, quantile, alpha / 2))
    bayes_cluster = c(bayes_cluster, res$factor)# this the result vector of BF after iterating thru every permutation
  }
  index = which.max(bayes_cluster)
  post_cluster_all[-arm_remain, m] = rep(1, N - N_remain)
  post_cluster_all[arm_remain, m] = all_cluster[index, ]
  reject_prob[arm_remain, m] = prob_rec[index, ]
  reject_acti[arm_remain, m] = acti_rec[index, ]
  post_prob_all[arm_remain, m] = prob_est[index, ]
  post_prob_upper_all[arm_remain, m] = prob_upper_rec[index, ]
  post_prob_lower_all[arm_remain, m] = prob_lower_rec[index, ]
  post_acti_all[arm_remain, m] = acti_est[index, ]
  post_acti_upper_all[arm_remain, m] = acti_upper_rec[index, ]
  post_acti_lower_all[arm_remain, m] = acti_lower_rec[index, ]
  
  setTxtProgressBar(pb, m / M)
}


setwd("~/Dropbox/Mayo-intern/HBM_Simulation/Results/500trials/binary/4act")
#prob = c(0.15, 0.15, 0.15, 0.45) ## true p
#acti = c(0.15, 0.15, 0.45, 0.45)  ## true activity
post_cluster_all = as.matrix(read.csv("cluster.csv")[, -1])
early_stop = as.matrix(read.csv("early.csv")[, -1])
post_prob_all = as.matrix(read.csv("prob.csv")[, -1])
post_prob_lower_all = as.matrix(read.csv("prob_lower.csv")[, -1])
post_prob_upper_all = as.matrix(read.csv("prob_upper.csv")[, -1])
post_acti_all = as.matrix(read.csv("acti.csv")[, -1])
post_acti_lower_all = as.matrix(read.csv("acti_lower.csv")[, -1])
post_acti_upper_all = as.matrix(read.csv("acti_upper.csv")[, -1])
reject_prob = as.matrix(read.csv("rej_prob.csv")[, -1])
reject_acti = as.matrix(read.csv("rej_acti.csv")[, -1])



setwd("~/Dropbox/Mayo-intern/HBM_Simulation/Results/500trials/binary/4act")
post_cluster_all = cbind(as.matrix(read.csv("2_cluster.csv")[, 2:201]), as.matrix(read.csv("4_cluster.csv")[, 202:401]), as.matrix(read.csv("5_cluster.csv")[, 402:501]))
early_stop = cbind(as.matrix(read.csv("2_early.csv")[, 2:201]), as.matrix(read.csv("4_early.csv")[, 202:401]), as.matrix(read.csv("5_early.csv")[, 402:501]))
post_prob_all = cbind(as.matrix(read.csv("2_prob.csv")[, 2:201]), as.matrix(read.csv("4_prob.csv")[, 202:401]), as.matrix(read.csv("5_prob.csv")[, 402:501]))
post_prob_lower_all = cbind(as.matrix(read.csv("2_prob_lower.csv")[, 2:201]), as.matrix(read.csv("4_prob_lower.csv")[, 202:401]), as.matrix(read.csv("5_prob_lower.csv")[, 402:501]))
post_prob_upper_all = cbind(as.matrix(read.csv("2_prob_upper.csv")[, 2:201]), as.matrix(read.csv("4_prob_upper.csv")[, 202:401]), as.matrix(read.csv("5_prob_upper.csv")[, 402:501]))
post_acti_all = cbind(as.matrix(read.csv("2_acti.csv")[, 2:201]), as.matrix(read.csv("4_acti.csv")[, 202:401]), as.matrix(read.csv("5_acti.csv")[, 402:501]))
post_acti_lower_all = cbind(as.matrix(read.csv("2_acti_lower.csv")[, 2:201]), as.matrix(read.csv("4_acti_lower.csv")[, 202:401]), as.matrix(read.csv("5_acti_lower.csv")[, 402:501]))
post_acti_upper_all = cbind(as.matrix(read.csv("2_acti_upper.csv")[, 2:201]), as.matrix(read.csv("4_acti_upper.csv")[, 202:401]), as.matrix(read.csv("5_acti_upper.csv")[, 402:501]))
reject_prob = cbind(as.matrix(read.csv("2_rej_prob.csv")[, 2:201]), as.matrix(read.csv("4_rej_prob.csv")[, 202:401]), as.matrix(read.csv("5_rej_prob.csv")[, 402:501]))
reject_acti = cbind(as.matrix(read.csv("2_rej_acti.csv")[, 2:201]), as.matrix(read.csv("4_rej_acti.csv")[, 202:401]), as.matrix(read.csv("5_rej_acti.csv")[, 402:501]))



## report
report = cbind(rowMeans(post_cluster_all == 1) * 100,
               rowMeans(post_cluster_all == 2) * 100,
               rowMeans(post_cluster_all == 3) * 100,
               rowMeans(early_stop) * 100, 
               rowMeans(post_prob_all, na.rm = TRUE),
               rowMeans(post_prob_lower_all, na.rm = TRUE),
               rowMeans(post_prob_upper_all, na.rm = TRUE),
               #rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob, na.rm = TRUE) * 100,
               rowMeans(post_acti_all, na.rm = TRUE),
               rowMeans(post_acti_lower_all, na.rm = TRUE),
               rowMeans(post_acti_upper_all, na.rm = TRUE),
               #rowMeans(post_acti_lower_all < acti & post_acti_upper_all > acti, na.rm = TRUE) * 100,
               rowMeans(reject_prob | reject_acti, na.rm = TRUE) * 100,
               rowMeans(reject_prob & reject_acti, na.rm = TRUE) * 100)
report = as.data.frame(report)
colnames(report) = c("C1", "C2", "C3", "early", "p_hat", "CI_l", "CI_u", "mu_hat", "CI_l", "CI_u", "weak", "strong")
report


xtable(report, digits = c(1, rep(1, 4), rep(2, 6), 1, 1))


### Plots of estimators

rst1 = c(post_prob_all[1, ], post_prob_all[2, ], post_prob_all[3, ], post_prob_all[4, ])
rst2 = c(post_acti_all[1, ], post_acti_all[2, ], post_acti_all[3, ], post_acti_all[4, ])
estimator = c(rep("Response", 2000), rep("Activity", 2000))
estimator = factor(estimator, levels = c("Response", "Activity"))
arm = rep(rep(c("Arm 1", "Arm 2", "Arm 3", "Arm 4"), each = 500), 2)
rst = data.frame("value" = c(rst1, rst2), "estimator" = estimator, "arm" = arm)

setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 5, height = 5)
ggplot(rst, aes(x = arm, y = value, fill = estimator)) + 
  geom_boxplot(lwd = 0.5, alpha = 1, width = 0.9, outlier.colour = "red", outlier.fill = "red", outlier.size = 2, outlier.alpha = 0) + 
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Parameters estimation") + 
  #scale_y_continuous(breaks = seq(0, 0.8, 0.2)) + 
  ylim(0, 0.8) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  theme(legend.position = "none")
dev.off()
tools::texi2dvi("plot.tex", pdf = T)

