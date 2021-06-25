
#### MCMC Sampling and calculate likelihood for the final stage
post = function(response, activity, ninter, group, cutoff, cutoff2, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu1_rec = mu2_rec = matrix(0, length(group), n.iter)
  
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
                                            mu1 = 0,
                                            mu2 = 0,
                                            prec = 1),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    update(thismodel, n.burn, progress.bar = "none") 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho", "prec"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
                               inits = list(mu1 = rep(0, N),
                                            mu2 = rep(0, N),
                                            prec = 1,
                                            rho = 0.5,
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
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho * (activity_sub[n, j] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
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
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
    thismodel = try(jags.model(file = "bugs/cont/c2.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(0, N),
                                            mu2 = rep(1, N),
                                            prec = 1,
                                            rho = 0.5,
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
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho * (activity_sub[n, j] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
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
    mu1 = matrix(res.bugs$mu1, nrow = 1)
    mu2 = matrix(res.bugs$mu2, nrow = 1)
    rho = matrix(res.bugs$rho, nrow = 1)
    prec = matrix(res.bugs$prec, nrow = 1)
    mu1_rec[ind, ] = mu1
    mu2_rec[ind, ] = mu2
    #p3 = dnorm(mu1, mean = -2, sd = 1 / sqrt(0.00001)) / pnorm(cutoff, mean = -2, sd = 1 / sqrt(0.00001))
    #p4 = dnorm(mu2, mean = 1, sd = 1 / sqrt(0.00001)) / pnorm(cutoff2, mean = 1, sd = 1 / sqrt(0.00001))
    #rst = rst + mean(log(p3), na.rm = T) + mean(log(p4), na.rm = T)
    for (j in 1:ninter) {
      p1 = dnorm(activity[ind, j], mean = mu2, sd = 1 / sqrt(prec), log = TRUE)
      p2 = pnorm(0, mean = mu1 + rho * (activity[ind, j] - mu2) * sqrt(prec), sd = sqrt(1 - rho^2))
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
               cutoff = cutoff)
    thismodel = try(jags.model(file = "bugs/cont/c3.txt", 
                               data = dat, 
                               inits = list(mu1 = rep(1, N),
                                            mu2 = rep(1, N),
                                            prec = 1,
                                            rho = 0.5,
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
        p1 = dnorm(activity_sub[n, j], mean = mu2[n, ], sd = 1 / sqrt(prec), log = TRUE)
        p2 = pnorm(0, mean = mu1[n, ] + rho * (activity_sub[n, j] - mu2[n, ]) * sqrt(prec), sd = sqrt(1 - rho^2))
        if (response_sub[n, j] == 1) {
          p2 = 1 - p2
        }
        rst = rst + mean(p1, na.rm = T) + mean(log(p2), na.rm = T)
      }
    }
  }
  return (list("factor" = rst, "mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec))
}





