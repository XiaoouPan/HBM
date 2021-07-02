library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(pbivnorm)
library(tikzDevice)
library(ggplot2)


rm(list = ls())

source('sas_v2.R')

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

  res = post_sas(response, activity, N, ninter, cluster, n.adapt, n.burn, n.iter)
  this_prob = pnorm(0, mean = qnorm(p0) + res$mu1_rec, sd = 1, lower.tail = FALSE)
  post_prob_all[, m] = as.numeric(rowMeans(this_prob))
  reject_prob[, m] = as.numeric(rowMeans(this_prob > p0) > reject_rate)
  post_prob_upper_all[, m] = apply(this_prob, 1, quantile, 0.975)
  post_prob_lower_all[, m] = apply(this_prob, 1, quantile, 0.025)
  this_acti = pnorm(0, mean = qnorm(a0) + res$mu2_rec, sd = 1, lower.tail = FALSE)
  post_acti_all[, m] = as.numeric(rowMeans(this_acti))
  reject_acti[, m] = as.numeric(rowMeans(this_acti > a0) > reject_rate)
  post_acti_upper_all[, m] = apply(this_acti, 1, quantile, 0.975)
  post_acti_lower_all[, m] = apply(this_acti, 1, quantile, 0.025)
  
  setTxtProgressBar(pb, m / M)
}


## report
report = cbind(prob,
               rowMeans(post_prob_all, na.rm = TRUE),
               rowMeans(post_prob_lower_all < prob & post_prob_upper_all > prob, na.rm = TRUE),
               acti,
               rowMeans(post_acti_all, na.rm = TRUE),
               rowMeans(post_acti_lower_all < acti & post_acti_upper_all > acti, na.rm = TRUE),
               rowMeans(reject_prob | reject_acti, na.rm = TRUE),
               rowMeans(reject_prob & reject_acti, na.rm = TRUE))
report = as.data.frame(report)
colnames(report) = c("true_p", "p_hat", "p_CI", "true_a", "a_hat", "a_CI", "weak", "strong")
report



#resp = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/Results/triCRM/resp_mix.csv")[, -1])
#acti = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/Results/triCRM/acti_mix.csv")[, -1])

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



