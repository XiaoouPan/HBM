#### Estimate the correlation as our first step
get_cor = function(response, activity, N, ninter, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter)
  thismodel = try(jags.model(file = "bugs/binary/cor.txt", 
                             data = dat, 
                             inits = list(mu1 = rep(0, N),
                                          mu2 = rep(0, N),
                                          prec = 1,
                                          rho = 0.5),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  update(thismodel, n.burn, progress.bar = "none") 
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "prec", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  return (as.vector(res.bugs$rho))
}


posterior_simu_s1 = function (dat, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  thismodel = try(jags.model(file = "trial_ct_s1.txt", 
                             data = dat, 
                             inits = list(mu2 = rep(1, dat$N),
                                          mumix2 = c(1, 5),
                                          muprec2 = c(1, 1)),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  update(thismodel, n.burn, progress.bar = "none") 
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu2', 'mumix2', 'muprec2'),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  return (list(mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               mumix2 = matrix(res.bugs$mumix2, nrow = 2),
               muprec2 = matrix(res.bugs$muprec2, nrow = 2)))
}

## Calculate Bayesian factors for each clustering permutation
## Activity is continuous, stage 1 based on activity
summary_posterior_s1 = function (dataVal, mcmcVal) {
  activity = dataVal$activity
  N = dataVal$N
  ninter = dataVal$ninter
  group = dataVal$group
  
  #parm
  mu2 = mcmcVal$mu2
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  rst = 0
  for (n in 1:N) {
    p4 = dnorm(mu2[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    rst = rst + ninter * mean(p4, na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[n, j], mean = mu2[n, ], sd = 1, log = TRUE)
      rst = rst + mean(p1, na.rm = T)
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



#####################################################
####  Biomarker activity is binary
#####################################################

## Generate parameters based on posterior distributons
## Activity is continuous, stage 1 based on activity
posterior_bi_simu_s1 = function (dat, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  thismodel = try(jags.model(file = "trial_bi_s1.txt", 
                             data = dat, 
                             inits = list(mu2 = rep(0, dat$N),
                                          mumix2 = c(-0.5, 0.5),
                                          muprec2 = c(1, 1)),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  update(thismodel, n.burn, progress.bar = "none")
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c('mu2', 'mumix2', 'muprec2'),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  return (list(mu2 = matrix(res.bugs$mu2, nrow = dat$N),
               mumix2 = matrix(res.bugs$mumix2, nrow = 2),
               muprec2 = matrix(res.bugs$muprec2, nrow = 2)))
}

## Calculate Bayesian factors for each clustering permutation
## Activity is continuous, stage 1 based on activity
summary_posterior_bi_s1 = function (dataVal, mcmcVal) {
  activity = dataVal$activity
  N = dataVal$N
  ninter = dataVal$ninter
  group = dataVal$group
  
  #parm
  mu2 = mcmcVal$mu2
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  rst = 0
  for (n in 1:N) {
    p4 = dnorm(mu2[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    rst = rst + ninter * mean(p4, na.rm = T)
    for (j in 1:ninter) {
      p1 = pnorm(0, mean = mu2[n, ])
      if (response[n, j] == 1) {
        p1 = 1 - p1
      }
      rst = rst + mean(log(p1), na.rm = T)
    }
  }
  return (rst)
}


## Generate parameters based on posterior distributons
## Activity is binary
posterior_bi_simu = function (dat, C, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  thismodel = try(jags.model(file = "trial_bi.txt", 
                             data = dat, 
                             inits = list(Z = array(c(dat$response, dat$activity), dim = c(dat$N, dat$ninter, 2)),
                                          mu1 = rep(-0.5, dat$N),
                                          mu2 = rep(0, dat$N),
                                          rho = rep(0.5, dat$N),
                                          mumix = c(-2, -2, 0),
                                          muprec = c(1, 1, 1),
                                          mumix2 = c(-0.5, 0.5, 0.5),
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
## Activity is binary
summary_posterior_bi = function (dataVal, mcmcVal) {
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
    p00 = pbivnorm(-mu1[n, ], -mu2[n, ], rho[n, ])
    p01 = pnorm(0, mean = mu1[n, ]) - p00
    p10 = pnorm(0, mean = mu2[n, ]) - p00
    p11 = 1 - p00 - p01 - p10
    for (j in 1:ninter) {
      if (response[n, j] == 0 & activity[n, j] == 0) {
        rst = rst + mean(log(p00), na.rm = T)
      } else if (response[n, j] == 0 & activity[n, j] == 1) {
        rst = rst + mean(log(p01), na.rm = T)
      } else if (response[n, j] == 1 & activity[n, j] == 0) {
        rst = rst + mean(log(p10), na.rm = T)
      } else {
        rst = rst + mean(log(p11), na.rm = T)
      }
    }
  }
  return (rst)
}

