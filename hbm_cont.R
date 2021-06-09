#### MCMC Sampling and likelihood for the interim stage
post_s1 = function(activity, n1, group, cutoff_int, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    dat = list(activity = activity[ind, ],
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1_uni.txt", 
                               data = dat, 
                               inits = list(mu2 = 0, 
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(p2, na.rm = T) + mean(p3, na.rm = T)
    for (j in 1:n1) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1 / sqrt(prec), log = TRUE)
      rst = rst + mean(p1, na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = activity_sub,
               N = N,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1.txt", 
                               data = dat, 
                               inits = list(mu2 = rep(0, N),
                                            prec = 1,
                                            mumix2 = 0,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "prec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T)
    for (n in 1:N) {
      p2 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        rst = rst + mean(p1, na.rm = T)
      }
    }
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N == 1) {
    dat = list(activity = activity[ind, ],
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2_uni.txt", 
                               data = dat, 
                               inits = list(mu2 = 5),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    p2 = dnorm(mu2, mean = 5, sd = 1 / sqrt(0.00001))
    rst = rst + mean(log(p2), na.rm = T)
    for (j in 1:n1) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1, log = TRUE)
      rst = rst + mean(p1, na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = activity_sub,
               N = N,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2.txt", 
                               data = dat, 
                               inits = list(mu2 = rep(5, N),
                                            mumix2 = 5,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c('mu2', 'mumix2', 'muprec2'),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    p3 = dnorm(mumix2, mean = 5, sd = 1 / sqrt(0.00001))
    p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T)
    for (n in 1:N) {
      p4 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(log(p4), na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
        rst = rst + mean(p1, na.rm = T)
      }
    }
  }
  return (rst)
}


#### MCMC Sampling and calculate likelihood for the final stage
post = function(response, activity, ninter, group, cutoff, cutoff2, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = -0.5,
                                            mu2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2), sd = 1)
      if (response[ind, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    response_sub = response[ind, ]
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1.txt", 
                               data = dat, 
                               inits = list(rho = rep(0.5, N),
                                            mu1 = rep(-2, N),
                                            mu2 = rep(1, N),
                                            mumix = -2,
                                            muprec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = N)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec)) / pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2)) / pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(log(p7), na.rm = T) + mean(log(p8), na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity_sub[n, j] - mu2[n, ]), sd = 1)
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N == 1) {
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = -0.5,
                                            mu2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2), sd = 1)
      if (response[ind, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    response_sub = response[ind, ]
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1.txt", 
                               data = dat, 
                               inits = list(rho = rep(0.5, N),
                                            mu1 = rep(-2, N),
                                            mu2 = rep(1, N),
                                            mumix = -2,
                                            muprec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = N)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec)) / pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2)) / pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(log(p7), na.rm = T) + mean(log(p8), na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity_sub[n, j] - mu2[n, ]), sd = 1)
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  
  ## Cluster 3
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = -0.5,
                                            mu2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1, log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2), sd = 1)
      if (response[ind, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    response_sub = response[ind, ]
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1.txt", 
                               data = dat, 
                               inits = list(rho = rep(0.5, N),
                                            mu1 = rep(-2, N),
                                            mu2 = rep(1, N),
                                            mumix = -2,
                                            muprec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = N)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec)) / pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2)) / pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(log(p7), na.rm = T) + mean(log(p8), na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity_sub[n, j] - mu2[n, ]), sd = 1)
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  return (rst)
}





## Generate parameters based on posterior distributons
## Activity is continuous, stage 2
posterior_simu = function (dat, C, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  thismodel = try(jags.model(file = "trial_ct.txt", 
                             data = dat, 
                             inits = list(mu1 = rep(-0.5, dat$N),
                                          mu2 = rep(1, dat$N),
                                          rho = rep(0.5, dat$N),
                                          mumix = c(-2, -2, 0),
                                          muprec = c(1, 1, 1),
                                          mumix2 = c(1, 5, 1),
                                          muprec2 = c(1, 1, 1)),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  update(thismodel, n.burn, progress.bar = "none")
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu1', 'mu2', 'rho', 'mumix', 'muprec', 'mumix2', 'muprec2'),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  return (list(mu1 = matrix(res.bugs$mu1, nrow = dat$N),
               mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               rho = matrix(res.bugs$rho, nrow = dat$N),
               mumix = matrix(res.bugs$mumix, nrow = C),
               muprec = matrix(res.bugs$muprec, nrow = C),
               mumix2 = matrix(res.bugs$mumix2, nrow = C),
               muprec2 = matrix(res.bugs$muprec2, nrow = C)))
}

## Calculate Bayesian factors for each clustering permutation
## Activity is continuous, stage 2
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
      p2 = pnorm(0, mean = mu1[n, ] + rho[n, ] * (activity[n, j] - mu2[n, ]), sd = 1)
      if (response[n, j] == 1) {
        p2 = 1 - p2
      }
      rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
    }
  }
  return (rst)
}


