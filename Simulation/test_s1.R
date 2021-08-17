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

response = matrix(0, N, ninter)
activity = matrix(0, N, ninter)
Z = matrix(0, ninter, 2) ## underlying bivariate normal, one of them is unobservable
Sigma = matrix(c(1, rho0, rho0, 1), 2, 2)
s1_cluster = permutations(n = 2, r = N, repeats.allowed = T)
post_cluster_all = matrix(0, N, M)


######## activity
epsilon_seq = seq(0, 1, by = 0.05)
l = length(epsilon_seq)
early2 = matrix(0, N, l)
pb = txtProgressBar(style = 3)
for (j in 1:l) {
  cutoff_int2 = epsilon_seq[j]
  for (m in 1:M) {
    set.seed(m)
    ## Data generation
    for (i in 1:N) {
      Z = mvrnorm(ninter, c(qnorm(prob)[i], acti[i]), Sigma)
      response[i, ] = as.numeric(Z[, 1] > 0)
      activity[i, ] = Z[, 2]
    }
    
    activity_s1 = activity[, 1:n1]
    bayes_cluster = NULL
    for (i in 1:nrow(s1_cluster)) {
      group = s1_cluster[i, ]
      res = post_s1_acti(activity_s1, n1, group, cutoff_int2, mu0, n.adapt, n.burn, n.iter)
      bayes_cluster = c(bayes_cluster, res$factor)
    }
    index = which.max(bayes_cluster)
    post_cluster_all[, m] = s1_cluster[index, ]
    
    setTxtProgressBar(pb, ((j - 1) * M +  m) / (l * M))
  }
  early2[, j] = rowMeans(post_cluster_all == 1)
}


######## response
epsilon_seq = seq(0, 0.3, by = 0.02)
l = length(epsilon_seq)
early1 = matrix(0, N, l)
pb = txtProgressBar(style = 3)
for (j in 1:l) {
  epsilon_1 = epsilon_seq[j]
  cutoff_int1 = qnorm(p0[1] + epsilon_1) - qnorm(p0[1])
  for (m in 1:M) {
    set.seed(m)
    ## Data generation
    for (i in 1:N) {
      Z = mvrnorm(ninter, c(qnorm(prob)[i], acti[i]), Sigma)
      response[i, ] = as.numeric(Z[, 1] > 0)
      activity[i, ] = Z[, 2]
    }
    
    response_s1 = response[, 1:n1]
    bayes_cluster = NULL
    for (i in 1:nrow(s1_cluster)) {
      group = s1_cluster[i, ]
      res = post_s1_resp(response_s1, n1, group, cutoff_int1, p0, n.adapt, n.burn, n.iter)
      bayes_cluster = c(bayes_cluster, res$factor)
    }
    index = which.max(bayes_cluster)
    post_cluster_all[, m] = s1_cluster[index, ]
    
    setTxtProgressBar(pb, ((j - 1) * M +  m) / (l * M))
  }
  early1[, j] = rowMeans(post_cluster_all == 1)
}

early1
early2


early1 = as.matrix(read.csv("~/Dropbox/Mayo-intern/HBM_Simulation/Results/Test_s2/s1_4act/early1.csv")[, -1])

### response
epsilon_seq = seq(0, 0.3, by = 0.02)
setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
plot(epsilon_seq, early1[1, ], type = "b", pch = 20, lwd = 5, cex = 1, axes = FALSE, ylim = c(0.1, 0.55), xlab = "", ylab = "")
lines(epsilon_seq, early1[2, ], type = "b", pch = 0, lwd = 5, cex = 1, col = "darkorange")
lines(epsilon_seq, early1[3, ], type = "b", pch = 4, lwd = 5, cex = 1, col = "forestgreen")
lines(epsilon_seq, early1[4, ], type = "b", pch = 5, lwd = 5, cex = 1, col = "blue")
#color = c("black", "darkorange", "forestgreen", "blue")
#labels = c("\\texttt{arm 1}", "\\texttt{arm 2}", "\\texttt{arm 3}", "\\texttt{arm 4}")
#pch = c(20, 0, 4, 5)
#legend("bottomright", labels, col = color, pch = pch, lwd = 5, cex = 2, box.lwd = 1, bg = "white")
axis(1, epsilon_seq[c(1, 6, 11, 16)], line = 0, cex.axis = 1.5)
axis(2, c(0.1, 0.25, 0.4, 0.55), line = 0, cex.axis = 1.5)
box()
abline(h = c(0.1, 0.25, 0.4, 0.55), v = epsilon_seq[c(1, 6, 11, 16)], col = "gray", lty = 2)
title(xlab = "Buffer $\\epsilon$", line = 2.5, cex.lab = 2)
title(ylab = "Early stopping rate", line = 2.5, cex.lab = 1.6)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



early2 = as.matrix(read.csv("~/Dropbox/Mayo-intern/HBM_Simulation/Results/Test_s2/s1_4act/early2.csv")[, -1])

### activity
epsilon_seq = seq(0, 1, by = 0.05)
setwd("~/Dropbox/Mayo-intern/HBM_Simulation")
tikz("plot.tex", standAlone = TRUE, width = 6, height = 5)
plot(epsilon_seq, early2[1, ], type = "b", pch = 20, lwd = 5, cex = 1, axes = FALSE, ylim = c(0, 0.55), xlab = "", ylab = "")
lines(epsilon_seq, early2[2, ], type = "b", pch = 0, lwd = 5, cex = 1, col = "darkorange")
lines(epsilon_seq, early2[3, ], type = "b", pch = 4, lwd = 5, cex = 1, col = "forestgreen")
lines(epsilon_seq, early2[4, ], type = "b", pch = 5, lwd = 5, cex = 1, col = "blue")
#color = c("black", "darkorange", "forestgreen", "blue")
#labels = c("\\texttt{arm 1}", "\\texttt{arm 2}", "\\texttt{arm 3}", "\\texttt{arm 4}")
#pch = c(20, 0, 4, 5)
#legend("bottomright", labels, col = color, pch = pch, lwd = 5, cex = 2, box.lwd = 1, bg = "white")
axis(1, epsilon_seq[c(1, 6, 11, 16, 21)], line = 0, cex.axis = 1.5)
axis(2, c(0, 0.15, 0.3, 0.45), line = 0, cex.axis = 1.5)
box()
abline(h = c(0, 0.15, 0.3, 0.45), v = epsilon_seq[c(1, 6, 11, 16, 21)], col = "gray", lty = 2)
title(xlab = "Buffer $\\epsilon$", line = 3, cex.lab = 2)
title(ylab = "Early stopping rate", line = 2.2, cex.lab = 1.6)
dev.off()
tools::texi2dvi("plot.tex", pdf = T)



