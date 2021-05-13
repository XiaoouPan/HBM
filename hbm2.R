library(MASS)
library(coda)
library(rjags) 
library(gtools)
library(psych)

rm(list = ls())

n = 100
mu = c(0, 3)
Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
dat = mvrnorm(n, mu, Sigma)
X = as.numeric(dat[, 1] > 0)
Y = dat[, 2]
biserial(Y, X)
