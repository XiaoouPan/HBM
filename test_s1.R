library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(mvtnorm)
library(tikzDevice)

rm(list = ls())

source('hbm_cont.R')

ninter = 22
n1 = 11
N = 4
C = 3
M = 5
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
epsilon_seq = seq(0, 1, by = 0.05)
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




######## response
epsilon_seq = seq(0, 0.3, by = 0.02)
l = length(epsilon_seq)
early1 = early2 = matrix(0, N, l)
pb = txtProgressBar(style = 3)
for (j in 1:l) {
  epsilon_1 = epsilon_seq[j]
  cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
  #cutoff_int2 = epsilon_2
  for (m in 1:M) {
    #set.seed(m)
    ## Data generation
    for (i in 1:N) {
      Z = mvrnorm(ninter, c(mu1[i], mu2[i]), Sigma)
      response[i, ] = as.numeric(Z[, 1] > 0)
      activity[i, ] = Z[, 2]
    }
    
    response_s1 = response[, 1:n1]
    bayes_cluster = NULL
    #prob_rec = prob_upper_rec = prob_lower_rec = NULL 
    #mu_rec = mu_upper_rec = mu_lower_rec = NULL 
    prob_rec = NULL
    for (i in 1:nrow(s1_cluster)) {
      group = s1_cluster[i, ]
      res = post_s1_resp(response_s1, n1, group, cutoff_int1)
      bayes_cluster = c(bayes_cluster, res$factor)
      this_prob = pnorm(0, mean = qnorm(p0) + res$mu1_rec, sd = 1, lower.tail = FALSE)
      prob_rec = rbind(prob_rec, as.numeric(rowMeans(this_prob > p0) > reject_rate))
    }
    index = which.max(bayes_cluster)
    post_cluster_all[, m] = s1_cluster[index, ]
    reject_prob[, m] = prob_rec[index, ]
    
    setTxtProgressBar(pb, ((j - 1) * M +  m) / (l * M))
  }
  early1[, j] = rowMeans(post_cluster_all == 1)
  early2[, j] = rowMeans(reject_prob)
}

early1
early2


early1 = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/early1_acti.csv")[, -1])
early2 = as.matrix(read.csv("~/Dropbox/Mayo-intern/Simulation/early2_acti.csv")[, -1])

setwd("~/Dropbox/Mayo-intern/Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
plot(epsilon_seq, early1[1, ], type = "b", pch = 20, lwd = 5, cex = 1, axes = FALSE, ylim = c(0.45, 1), xlab = "", ylab = "")
lines(epsilon_seq, early1[2, ], type = "b", pch = 0, lwd = 5, cex = 1, col = "darkorange")
lines(epsilon_seq, early1[3, ], type = "b", pch = 4, lwd = 5, cex = 1, col = "forestgreen")
lines(epsilon_seq, early1[4, ], type = "b", pch = 5, lwd = 5, cex = 1, col = "blue")
color = c("black", "darkorange", "forestgreen", "blue")
labels = c("\\texttt{arm 1}", "\\texttt{arm 2}", "\\texttt{arm 3}", "\\texttt{arm 4}")
pch = c(20, 0, 4, 5)
legend("bottomright", labels, col = color, pch = pch, lwd = 5, cex = 1.7, box.lwd = 1, bg = "white")
axis(1, epsilon_seq[c(1, 6, 11, 16, 21)], line = 0, cex.axis = 1.5)
axis(2, c(0.5, 0.65, 0.8, 0.95), line = 0, cex.axis = 1.5)
box()
abline(h = c(0.5, 0.65, 0.8, 0.95), v = epsilon_seq[c(1, 6, 11, 16, 21)], col = "gray", lty = 2)
title(xlab = "Buffer $\\epsilon$", line = 3, cex.lab = 2)
title(ylab = "Early stopping rate", line = 2.2, cex.lab = 1.6)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
plot(epsilon_seq, 1 - early2[1, ], type = "b", pch = 20, lwd = 5, cex = 1, axes = FALSE, ylim = c(0.45, 1), xlab = "", ylab = "")
lines(epsilon_seq, 1 - early2[2, ], type = "b", pch = 0, lwd = 5, cex = 1, col = "darkorange")
lines(epsilon_seq, 1 - early2[3, ], type = "b", pch = 4, lwd = 5, cex = 1, col = "forestgreen")
lines(epsilon_seq, 1 - early2[4, ], type = "b", pch = 5, lwd = 5, cex = 1, col = "blue")
color = c("black", "darkorange", "forestgreen", "blue")
labels = c("\\texttt{arm 1}", "\\texttt{arm 2}", "\\texttt{arm 3}", "\\texttt{arm 4}")
pch = c(20, 0, 4, 5)
legend("bottomright", labels, col = color, pch = pch, lwd = 5, cex = 1.7, box.lwd = 1, bg = "white")
axis(1, epsilon_seq[c(1, 6, 11, 16, 21)], line = 0, cex.axis = 1.5)
axis(2, c(0.5, 0.65, 0.7, 0.85, 1), line = 0, cex.axis = 1.5)
box()
abline(h = c(0.5, 0.65, 0.7, 0.85, 1), v = epsilon_seq[c(1, 6, 11, 16, 21)], col = "gray", lty = 2)
title(xlab = "Buffer $\\epsilon$", line = 3, cex.lab = 2)
title(ylab = "Early stopping rate", line = 2.2, cex.lab = 1.6)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



