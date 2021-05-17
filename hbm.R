## Generate parameters based on posterior distributons
posterior_simu = function (dat, C, iter = 2000) {
  thismodel = try(jags.model(file = "trial.txt", 
                             data = dat, 
                             inits = list(mu1 = rep(-0.5, dat$N),
                                          mu2 = rep(1, dat$N),
                                          rho = rep(0.5, dat$N),
                                          mumix = c(-2, -2, 0),
                                          muprec = c(1, 1, 1),
                                          mumix2 = c(1, 5, 1),
                                          muprec2 = c(1, 1, 1)),
                             n.adapt = iter), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu1', 'mu2', 'rho', 'mumix', 'muprec', 'mumix2', 'muprec2'),
                              n.iter = iter), silent = TRUE)
  return (list(mu1 = matrix(res.bugs$mu1, nrow = dat$N),
               mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               rho = matrix(res.bugs$rho, nrow = dat$N),
               mumix = matrix(res.bugs$mumix, nrow = C),
               muprec = matrix(res.bugs$muprec, nrow = C),
               mumix2 = matrix(res.bugs$mumix2, nrow = C),
               muprec2 = matrix(res.bugs$muprec2, nrow = C)))
}

# Calculate Bayesian factors for each clustering permutation
summary_posterior = function (dataVal, mcmcVal) {
  response = dataVal$response
  activity = dataVal$activity
  N = dataVal$N
  ninter = dataVal$ninter
  group = dataVal$group
  
  #parm
  mu1 = mcmcVal$mu1
  mu2 = mcmcVal$mu2
  rho = mcmcVal$rho
  mumix = mcmcVal$mumix
  muprec = mcmcVal$muprec
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  rst = 0
  for (n in 1:N) {
    p3 = dnorm(mu1[n, ], mean = mumix[group[n], ], sd = 1 / sqrt(muprec[group[n], ]), log = T)
    p4 = dnorm(mu2[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    rst = rst + ninter * (mean(p3, na.rm = T) + mean(p4, na.rm = T))
    for (j in 1:ninter) {
      p1 = dnorm(activity[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity[n, j] - mu2[n, ]), sd = 1, log.p = TRUE)
      if (response[n, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(p2, na.rm = T)
    }
  }
  return (rst)
}

