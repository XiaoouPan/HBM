library(ggplot2)
library(tikzDevice)

rm(list = ls())

n = 5000
m = 100
J = 4
var.pri = seq(2, 50, length.out = m)
p0 = c(0.15, 0.15, 0.18, 0.18)
pa = p0 + 0.3
ess.all = matrix(0, J, m)
ess1 = rep(0, m)
for (i in 1:m) {
  u = rnorm(n, 1, var.pri[i])
  for (j in 1:J) {
    w = rbinom(n, 1, 0.5)
    mu = w * u
    prob = pnorm(qnorm(p0[j]) + mu)
    E = mean(prob)
    V = var(prob)
    a = (1 - E) * E^2 / V - E
    b = a * (1 - E) / E
    ess.all[j, i] = a + b
  }
}
ess1 = colMeans(ess.all)
p0 = c(0.12, 0.12, 0.15, 0.15)
pa = p0 + 0.3
ess.all = matrix(0, J, m)
ess2 = rep(0, m)
for (i in 1:m) {
  u = rnorm(n, 1, var.pri[i])
  for (j in 1:J) {
    w = rbinom(n, 1, 0.5)
    mu = w * u
    prob = pnorm(qnorm(p0[j]) + mu)
    E = mean(prob)
    V = var(prob)
    a = (1 - E) * E^2 / V - E
    b = a * (1 - E) / E
    ess.all[j, i] = a + b
  }
}
ess2 = colMeans(ess.all)

var.pri = rep(var.pri, 4)
ess = c(ess1, ess2, rep(0.4, m), rep(0.6, m))
dat = as.data.frame(cbind(var.pri, ess))
colnames(dat) = c("var", "ess")
dat$type = c(rep("ess1", m), rep("ess2", m), rep("lower", m),  rep("upper", m))
dat$type = factor(dat$type, levels = c("ess1", "ess2", "lower", "upper"))

setwd("~/Dropbox/Mayo-intern/Spike-and-slab")
tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
ggplot(dat, aes(x = var, y = ess, color = type)) +
  geom_line(aes(y = ess, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dashed")) +
  theme_bw() + xlab("Candidates of $\\sigma^2$") + 
  ylab("Approximate effective sample size") + 
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


