library(ggplot2)
library(tikzDevice)

rm(list = ls())

n = 5000
p0 = 0.15
m = 100
var.pri = seq(2, 50, length.out = m)
ess = rep(0, m)
pb = txtProgressBar(style = 3)
for (i in 1:m) {
  w = rbinom(n, 1, 0.5)
  mu = w * rnorm(n, 1, var.pri[i])
  prob = pnorm(qnorm(p0) + mu)
  E = mean(prob)
  V = var(prob)
  a = (1 - E) * E^2 / V - E
  b = a * (1 - E) / E
  ess[i] = a + b
  setTxtProgressBar(pb, i / m)
}

var.pri = rep(var.pri, 3)
ess = c(ess, rep(0.5, m), rep(0.6, m))
dat = as.data.frame(cbind(var.pri, ess))
colnames(dat) = c("var", "ess")
dat$type = c(rep("ess", m), rep("lower", m),  rep("upper", m))
dat$type = factor(dat$type, levels = c("ess", "lower", "upper"))

setwd("~/Dropbox/Mayo-intern/Spike-and-slab")
tikz("plot.tex", standAlone = TRUE, width = 7, height = 5)
ggplot(dat, aes(x = var, y = ess, color = type)) +
  geom_line(aes(y = ess, color = type, linetype = type), size = 3) + 
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  theme_bw() + xlab("Candidates of $\\sigma^2$") + 
  ylab("Effective sample size") + 
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
dev.off()
tools::texi2dvi("plot.tex", pdf = T)


