library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)
library(dplyr)
library(tikzDevice)
library(ggplot2)
library(xtable)


rm(list = ls())

source('sas_v4.R')

ninter = 22
n1 = 11
N = 4
M = 5
n.adapt = 5000
n.burn = 5000
n.iter = 10000

p0 = c(0.15, 0.15, 0.15, 0.15) ## null response rate
a0 = c(0.15, 0.15, 0.15, 0.15) ## null activity level
rho0 = 0.75
alpha = 0.026
reject_rate = 1 - alpha ## For hypothesis testing

prob = c(0.15, 0.15, 0.15, 0.45) ## true p
acti = c(0.15, 0.15, 0.15, 0.45)  ## true activity
mu1 = qnorm(prob) - qnorm(p0)
mu2 = qnorm(acti) - qnorm(a0)

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)

#early_stop = matrix(0, N, M)
early_stop = reject_weak = reject_strong = matrix(0, N, M)
post_prob_all = post_prob_upper_all = post_prob_lower_all = matrix(NA, N, M)
post_acti_all = post_acti_upper_all = post_acti_lower_all = matrix(NA, N, M)
#cluster = getCluster(N) ## 15 or 41 possibilities in total
#trans = getTrans(cluster)

pb = txtProgressBar(style = 3)
for (m in 1:M) {
  set.seed(m)
  ## Data generation
  for (i in 1:N) {
    Z = mvrnorm(ninter, c(qnorm(prob)[i], qnorm(acti)[i]), Sigma)
    response[i, ] = as.numeric(Z[, 1] > 0)
    activity[i, ] = as.numeric(Z[, 2] > 0)
  }
  
  ## Final stage
  arm_remain = c(1, 2, 3, 4)
  #arm_remain = which(reject_acti[, m] == 1)
  N_remain = length(arm_remain)
  #early_stop[-arm_remain, m] = rep(1, N - N_remain)
  response_remain = response[arm_remain, , drop = FALSE]
  activity_remain = activity[arm_remain, , drop = FALSE]
  
  res = post_sas(response_remain, activity_remain, N_remain, ninter, p0, mu1_h0 = qnorm(p0)[arm_remain], a0, mu2_h0 = qnorm(a0)[arm_remain], n.adapt, n.burn, n.iter)
  this_prob = pnorm(0, mean = qnorm(p0)[arm_remain] + res$mu1_rec, sd = 1, lower.tail = FALSE)
  post_prob_all[arm_remain, m] = as.numeric(rowMeans(this_prob))
  post_prob_upper_all[arm_remain, m] = apply(this_prob, 1, quantile, 1 - alpha / 2)
  post_prob_lower_all[arm_remain, m] = apply(this_prob, 1, quantile, alpha / 2)
  this_acti = pnorm(0, mean = qnorm(a0)[arm_remain] + res$mu2_rec, sd = 1, lower.tail = FALSE)
  post_acti_all[arm_remain, m] = as.numeric(rowMeans(this_acti))
  post_acti_upper_all[arm_remain, m] = apply(this_acti, 1, quantile, 1 - alpha / 2)
  post_acti_lower_all[arm_remain, m] = apply(this_acti, 1, quantile, alpha / 2)
  reject_weak[arm_remain, m] = as.numeric(rowMeans(this_prob > p0[arm_remain] | this_acti > a0[arm_remain]) > reject_rate)
  reject_strong[arm_remain, m] = as.numeric(rowMeans(this_prob > p0[arm_remain] & this_acti > a0[arm_remain]) > reject_rate)
  
  setTxtProgressBar(pb, m / M)
}


## report
report = cbind(rowMeans(post_prob_all, na.rm = TRUE),
               rowMeans(post_prob_lower_all, na.rm = TRUE),
               rowMeans(post_prob_upper_all, na.rm = TRUE),
               #rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob, na.rm = TRUE) * 100,
               rowMeans(post_acti_all, na.rm = TRUE),
               rowMeans(post_acti_lower_all, na.rm = TRUE),
               rowMeans(post_acti_upper_all, na.rm = TRUE),
               #rowMeans(post_acti_lower_all < acti & post_acti_upper_all > acti, na.rm = TRUE) * 100,
               rowMeans(reject_weak, na.rm = TRUE) * 100,
               rowMeans(reject_strong, na.rm = TRUE) * 100)
report = as.data.frame(report)
colnames(report) = c("p_hat", "CI_l", "CI_u", "mu_hat", "CI_l", "CI_u", "weak", "strong")
report




#### Results
setwd("~/Dropbox/Mayo-intern/HBM_Simulation/Results/sas/mode")
#prob = c(0.15, 0.15, 0.15, 0.45) ## true p
#acti = c(3, 3, 4, 4)  ## true activity
post_prob_all = as.matrix(read.csv("prob.csv")[, -1])
post_prob_lower_all = as.matrix(read.csv("prob_lower.csv")[, -1])
post_prob_upper_all = as.matrix(read.csv("prob_upper.csv")[, -1])
post_acti_all = as.matrix(read.csv("acti.csv")[, -1])
post_acti_lower_all = as.matrix(read.csv("acti_lower.csv")[, -1])
post_acti_upper_all = as.matrix(read.csv("acti_upper.csv")[, -1])
reject_weak = as.matrix(read.csv("rej_weak.csv")[, -1])
reject_strong = as.matrix(read.csv("rej_strong.csv")[, -1])

report = cbind(rowMeans(post_prob_all, na.rm = TRUE),
               rowMeans(post_prob_lower_all, na.rm = TRUE),
               rowMeans(post_prob_upper_all, na.rm = TRUE),
               #rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob, na.rm = TRUE) * 100,
               rowMeans(post_acti_all, na.rm = TRUE),
               rowMeans(post_acti_lower_all, na.rm = TRUE),
               rowMeans(post_acti_upper_all, na.rm = TRUE),
               #rowMeans(post_acti_lower_all < acti & post_acti_upper_all > acti, na.rm = TRUE) * 100,
               rowMeans(reject_weak, na.rm = TRUE) * 100,
               rowMeans(reject_strong, na.rm = TRUE) * 100)
report = as.data.frame(report)
colnames(report) = c("p_hat", "CI_l", "CI_u", "mu_hat", "CI_l", "CI_u", "weak", "strong")
report

xtable(report, digits = c(1, rep(2, 6), 1, 1))


M = 100
rst1 = c(post_prob_all[1, ], post_prob_all[2, ], post_prob_all[3, ], post_prob_all[4, ])
rst2 = c(post_acti_all[1, ], post_acti_all[2, ], post_acti_all[3, ], post_acti_all[4, ])
estimator = c(rep("Response", 4 * M), rep("Activity", 4 * M))
estimator = factor(estimator, levels = c("Response", "Activity"))
arm = rep(rep(c("Arm 1", "Arm 2", "Arm 3", "Arm 4"), each = M), 2)
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



