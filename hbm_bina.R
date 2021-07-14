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
                                          rho = 0.5),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  update(thismodel, n.burn, progress.bar = "none") 
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  return (as.vector(res.bugs$rho))
}

#### MCMC Sampling and likelihood for the interim stage, using response (sum of response)
post_s1_resp = function(response, n1, group, cutoff_int, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu1_rec = matrix(0, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    r1 = sum(response[ind, ])
    dat = list(response = r1,
               ninter = n1,
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c1_uni_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = 0),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu1_rec[ind, ] = mu1
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p0 = pnorm(0, mean = mu1, sd = 1)
    rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
  } else if (N >= 2) {
    response_sub = response[ind, ]
    dat = list(response = rowSums(response_sub),
               N = N,
               ninter = n1,
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c1_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(0, N),
                                            mumix1 = 0,
                                            muprec1 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mumix1", "muprec1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu1_rec[ind, ] = mu1
    mumix1 = matrix(res.bugs$mumix1, nrow = 1)
    muprec1 = matrix(res.bugs$muprec1, nrow = 1)
    #p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu1[n, ], mean = mumix1, sd = 1 / sqrt(muprec1), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      r1 = sum(response_sub[n, ])
      p0 = pnorm(0, mean = mu1[n, ], sd = 1)
      rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
    }
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N == 1) {
    r1 = sum(response[ind, ])
    dat = list(response = r1,
               ninter = n1,
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c2_uni_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu1_rec[ind, ] = mu1
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p0 = pnorm(0, mean = mu1, sd = 1)
    rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
  } else if (N >= 2) {
    response_sub = response[ind, ]
    dat = list(response = rowSums(response_sub),
               N = N,
               ninter = n1,
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c2_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(1, N),
                                            mumix1 = 1,
                                            muprec1 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mumix1", "muprec1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu1_rec[ind, ] = mu1
    mumix1 = matrix(res.bugs$mumix1, nrow = 1)
    muprec1 = matrix(res.bugs$muprec1, nrow = 1)
    #p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu1[n, ], mean = mumix1, sd = 1 / sqrt(muprec1), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      r1 = sum(response_sub[n, ])
      p0 = pnorm(0, mean = mu1[n, ], sd = 1)
      rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
    }
  }
  return (list("factor" = rst, "mu1_rec" = mu1_rec))
}


#### MCMC Sampling and likelihood for the interim stage, using activity (sum of activity)
post_s1_acti = function(activity, n1, group, cutoff_int, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu2_rec = matrix(0, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    r1 = sum(activity[ind, ])
    dat = list(activity = r1,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c1_uni_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = 0),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    mu2_rec[ind, ] = mu2
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p0 = pnorm(0, mean = mu2, sd = 1)
    rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = rowSums(activity_sub),
               N = N,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c1_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = rep(0, N),
                                            mumix2 = 0,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu1, nrow = N)
    mu2_rec[ind, ] = mu2
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    #p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      r1 = sum(activity_sub[n, ])
      p0 = pnorm(0, mean = mu2[n, ], sd = 1)
      rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
    }
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N == 1) {
    r1 = sum(activity[ind, ])
    dat = list(activity = r1,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c2_uni_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    mu2_rec[ind, ] = mu2
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p0 = pnorm(0, mean = mu2, sd = 1)
    rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = rowSums(activity_sub),
               N = N,
               ninter = n1,
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/binary/int_c2_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = rep(1, N),
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    mu2_rec[ind, ] = mu2
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    #p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      r1 = sum(activity_sub[n, ])
      p0 = pnorm(0, mean = mu2[n, ], sd = 1)
      rst = rst + (n1 - r1) * mean(log(p0), na.rm = T) + r1 * mean(log(1 - p0), na.rm = T)
    }
  }
  return (list("factor" = rst, "mu2_rec" = mu2_rec))
}

















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

