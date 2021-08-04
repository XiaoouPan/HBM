#### Estimate the correlation as our first step
get_cor = function(response, activity, N, ninter, p0, mu0, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter,
             p_h0 = qnorm(p0),
             mu_h0 = mu0)
  thismodel = try(jags.model(file = "bugs/cont/cor.txt", 
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


#### MCMC Sampling and likelihood for the interim stage, using response (sum of response)
post_s1_resp = function(response, n1, group, cutoff_int, p0, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu1_rec = matrix(0, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    r1 = sum(response[ind, ])
    dat = list(response = r1,
               ninter = n1,
               p_h0 = qnorm(p0)[ind],
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1_uni_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = 0),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu1_rec[ind, ] = mu1
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p1 = pnorm(0, mean = qnorm(p0)[ind] + mu1, sd = 1)
    if (r1 < n1) {
      rst = rst + (n1 - r1) * mean(log(p1), na.rm = T)
    }
    if (r1 > 0) {
      rst = rst + r1 * mean(log(1 - p1), na.rm = T)
    }
  } else if (N >= 2) {
    response_sub = response[ind, ]
    dat = list(response = rowSums(response_sub),
               N = N,
               ninter = n1,
               p_h0 = qnorm(p0)[ind],
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1_resp.txt", 
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
      p1 = pnorm(0, mean = qnorm(p0)[ind][n] + mu1[n, ], sd = 1)
      if (r1 < n1) {
        rst = rst + (n1 - r1) * mean(log(p1), na.rm = T)
      }
      if (r1 > 0) {
        rst = rst + r1 * mean(log(1 - p1), na.rm = T)
      }
    }
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N == 1) {
    r1 = sum(response[ind, ])
    dat = list(response = r1,
               ninter = n1,
               p_h0 = qnorm(p0)[ind],
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2_uni_resp.txt", 
                               data = dat, 
                               inits = list(mu1 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu1_rec[ind, ] = mu1
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    p1 = pnorm(0, mean =  qnorm(p0)[ind] + mu1, sd = 1)
    if (r1 < n1) {
      rst = rst + (n1 - r1) * mean(log(p1), na.rm = T)
    }
    if (r1 > 0) {
      rst = rst + r1 * mean(log(1 - p1), na.rm = T)
    }
  } else if (N >= 2) {
    response_sub = response[ind, ]
    dat = list(response = rowSums(response_sub),
               N = N,
               ninter = n1,
               p_h0 = qnorm(p0)[ind],
               cutoff1 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2_resp.txt", 
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
      p3 = pnorm(0, mean = qnorm(p0)[ind][n] + mu1[n, ], sd = 1)
      if (r1 < n1) {
        rst = rst + (n1 - r1) * mean(log(p3), na.rm = T)
      }
      if (r1 > 0) {
        rst = rst + r1 * mean(log(1 - p3), na.rm = T)
      }
    }
  }
  return (list("factor" = rst, "mu1_rec" = mu1_rec))
}


#### MCMC Sampling and likelihood for the interim stage, using activity
post_s1_acti = function(activity, n1, group, cutoff_int, mu0, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu2_rec = matrix(0, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    dat = list(activity = activity[ind, ],
               ninter = n1,
               mu_h0 = mu0[ind],
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1_uni_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = 0, 
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = as.numeric(res.bugs$mu2)
    mu2_rec[ind, ] = mu2
    prec = matrix(res.bugs$prec, nrow = 1)
    #p2 = dnorm(mu2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    for (j in 1:n1) {
      p1 = dnorm(activity[ind, j], mean = mu0[ind] + mu2, sd = 1 / sqrt(prec), log = TRUE)
      rst = rst + mean(p1, na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = activity_sub,
               N = N,
               ninter = n1,
               mu_h0 = mu0[ind],
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c1_acti.txt", 
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
    mu2_rec[ind, ] = mu2
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    #p3 = dnorm(mumix2, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 0, sd = 1 / sqrt(0.00001), log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu0[ind][n] + mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
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
               mu_h0 = mu0[ind],
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2_uni_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = 1,
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = as.numeric(res.bugs$mu2)
    mu2_rec[ind, ] = mu2
    prec = matrix(res.bugs$prec, nrow = 1)
    #p2 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 1, sd = 1 / sqrt(0.00001), lower.tail = FALSE, log = TRUE)
    #p3 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p2, na.rm = T) + mean(p3, na.rm = T))
    for (j in 1:n1) {
      p1 = dnorm(activity[ind, j], mean = mu0[ind] + mu2, sd = 1 / sqrt(prec), log = TRUE)
      rst = rst + mean(p1, na.rm = T)
    }
  } else if (N >= 2) {
    activity_sub = activity[ind, ]
    dat = list(activity = activity_sub,
               N = N,
               ninter = n1,
               mu_h0 = mu0[ind],
               cutoff2 = cutoff_int)
    thismodel = try(jags.model(file = "bugs/cont/int_c2_acti.txt", 
                               data = dat, 
                               inits = list(mu2 = rep(1, N),
                                            prec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu2", "prec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    mu2_rec[ind, ] = mu2
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    #p3 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001), log = TRUE) - pnorm(cutoff_int, mean = 1, sd = 1 / sqrt(0.00001), lower.tail = FALSE, log = TRUE)
    #p4 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dgamma(prec, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + (mean(p3, na.rm = T) + mean(p4, na.rm = T) + mean(p5, na.rm = T))
    for (n in 1:N) {
      p2 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #- pnorm(cutoff_int, mean = mumix2, sd = 1 / sqrt(muprec2), lower.tail = FALSE, log = TRUE)
      rst = rst + mean(p2, na.rm = T)
      for (j in 1:n1) {
        p1 = dnorm(activity_sub[n, j], mean = mu0[ind][n] + mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        rst = rst + mean(p1, na.rm = T)
      }
    }
  }
  return (list("factor" = rst, "mu2_rec" = mu2_rec))
}


#### MCMC Sampling and calculate likelihood for the final stage
post = function(response, activity, ninter, group, cutoff, cutoff2, p0, mu0, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu1_rec = mu2_rec = matrix(0, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N == 1) {
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               ninter = ninter,
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = 0,
                                            mu2 = 0,
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu0[ind] + mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = qnorm(p0)[ind] + mu1 + rho * (activity[ind, j] - mu0[ind] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c1.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(0, N),
                                            mu2 = rep(0, N),
                                            prec = 1,
                                            rho = rep(0.5, N),
                                            mumix = 0,
                                            muprec = 1,
                                            mumix2 = 0,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec), log = TRUE) #/ pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #/ pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(p7, na.rm = T) + mean(p8, na.rm = T)
      for (j in 1:ninter) {
        p1 = dnorm(activity_sub[n, j], mean = mu0[ind][n] + mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = qnorm(p0)[ind][n] + mu1[n, ] + rho * (activity_sub[n, j] - mu0[ind][n] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
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
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c2_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = 0,
                                            mu2 = 1,
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu0[ind] + mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = qnorm(p0)[ind] + mu1 + rho * (activity[ind, j] - mu0[ind] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/cont/c2.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(0, N),
                                            mu2 = rep(1, N),
                                            prec = 1,
                                            rho = rep(0.5, N),
                                            mumix = 0,
                                            muprec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec), log = TRUE) #/ pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #/ pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(p7, na.rm = T) + mean(p8, na.rm = T)
      for (j in 1:ninter) {
        p1 = dnorm(activity_sub[n, j], mean = mu0[ind][n] + mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = qnorm(p0)[ind][n] + mu1[n, ] + rho * (activity_sub[n, j] - mu0[ind][n] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  
  ## Cluster 3
  ind = which(group == 3)
  N = length(ind)
  if (N == 1) {
    dat = list(response = response[ind, ],
               activity = activity[ind, ],
               ninter = ninter,
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff)
    thismodel = try(jags.model(file = "bugs/cont/c3_uni.txt", 
                               data = dat, 
                               inits = list(rho = 0.5,
                                            mu1 = 1,
                                            mu2 = 1,
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu0[ind] + mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = qnorm(p0)[ind] + mu1 + rho * (activity[ind, j] - mu0[ind] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
               p_h0 = qnorm(p0)[ind],
               mu_h0 = mu0[ind], 
               cutoff = cutoff)
    thismodel = try(jags.model(file = "bugs/cont/c3.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(1, N),
                                            mu2 = rep(1, N),
                                            prec = 1,
                                            rho = rep(0.5, N),
                                            mumix = 1,
                                            muprec = 1,
                                            mumix2 = 1,
                                            muprec2 = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec", "mumix", "muprec", "mumix2", "muprec2"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = N)
    mu2 = matrix(res.bugs$mu2, nrow = N)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mumix = matrix(res.bugs$mumix, nrow = 1)
    muprec = matrix(res.bugs$muprec, nrow = 1)
    mumix2 = matrix(res.bugs$mumix2, nrow = 1)
    muprec2 = matrix(res.bugs$muprec2, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mumix, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dgamma(muprec, shape = 0.001, rate = 0.001, log = TRUE)
    #p5 = dnorm(mumix2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #p6 = dgamma(muprec2, shape = 0.001, rate = 0.001, log = TRUE)
    #rst = rst + mean(log(p3), na.rm = T) + mean(p4, na.rm = T) + mean(log(p5), na.rm = T) + mean(p6, na.rm = T)
    for (n in 1:N) {
      p7 = dnorm(mu1[n, ], mean = mumix, sd = 1 / sqrt(muprec), log = TRUE) #/ pnorm(cutoff, mean = mumix, sd = 1 / sqrt(muprec))
      p8 = dnorm(mu2[n, ], mean = mumix2, sd = 1 / sqrt(muprec2), log = TRUE) #/ pnorm(cutoff2, mean = mumix2, sd = 1 / sqrt(muprec2))
      rst = rst + mean(p7, na.rm = T) + mean(p8, na.rm = T)
      for (j in 1:ninter) {
        p1 = dnorm(activity_sub[n, j], mean = mu0[ind][n] + mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = qnorm(p0)[ind][n] + mu1[n, ] + rho * (activity_sub[n, j] - mu0[ind][n] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  return (list("factor" = rst, "mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec))
}





