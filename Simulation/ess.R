rm(list = ls())

n = 1000
p0 = 0.15
m = 500
var.pri = seq(1, 50, length.out = m)
ess = rep(0, m)
pb = txtProgressBar(style = 3)
for (i in 1:m) {
  w = rbinom(n, 1, 0.5)
  mu = w * rnorm(1000, 1, var.pri[i])
  prob = pnorm(qnorm(p0) + mu)
  E = mean(prob)
  V = var(prob)
  a = (1 - E) * E^2 / V - E
  b = a * (1 - E) / E
  ess[i] = a + b
  setTxtProgressBar(pb, i / m)
}
idx = which(ess > 0.55 & ess < 0.6)
var.pri[idx]
plot(var.pri, ess, type = "l")
