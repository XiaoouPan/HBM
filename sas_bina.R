
#### MCMC Sampling and calculate likelihood for the final stage
post = function(response, activity, ninter, group, cutoff, cutoff2, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  mu1_rec = mu2_rec = matrix(NA, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N >= 1) {
    activity_sub = activity[ind, , drop = FALSE]
    response_sub = response[ind, , drop = FALSE]
    dat = list(response = response_sub,
               activity = activity_sub,
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/sas_binary/c1.txt", 
                               data = dat, 
                               inits = list(Z = array(c(dat$response, dat$activity), dim = c(dat$N, dat$ninter, 2)),
                                            mu1 = 0,
                                            mu2 = 0,
                                            rho = 0.5),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = as.numeric(res.bugs$rho)
    mu1_rec[ind, ] = matrix(rep(mu1, N), nrow = N, byrow = TRUE)
    mu2_rec[ind, ] = matrix(rep(mu2, N), nrow = N, byrow = TRUE)
    p00 = pbivnorm(-mu1, -mu2, rho)
    p01 = pnorm(0, mean = mu1) - p00
    p10 = pnorm(0, mean = mu2) - p00
    p11 = 1 - p00 - p01 - p10
    s00 = sum(response_sub == 0 & activity_sub == 0)
    s01 = sum(response_sub == 0 & activity_sub == 1)
    s10 = sum(response_sub == 1 & activity_sub == 0)
    s11 = sum(response_sub == 1 & activity_sub == 1)
    rst = rst + s00 * mean(log(p00), na.rm = T) + s01 * mean(log(p01), na.rm = T) + s10 * mean(log(p10), na.rm = T) + s11 * mean(log(p11), na.rm = T)
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N >= 1) {
    activity_sub = activity[ind, , drop = FALSE]
    response_sub = response[ind, , drop = FALSE]
    dat = list(response = response_sub,
               activity = activity_sub,
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/sas_binary/c2.txt", 
                               data = dat, 
                               inits = list(Z = array(c(dat$response, dat$activity), dim = c(dat$N, dat$ninter, 2)),
                                            mu1 = 0,
                                            mu2 = 1,
                                            rho = 0.5),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = as.numeric(res.bugs$rho)
    mu1_rec[ind, ] = matrix(rep(mu1, N), nrow = N, byrow = TRUE)
    mu2_rec[ind, ] = matrix(rep(mu2, N), nrow = N, byrow = TRUE)
    p00 = pbivnorm(-mu1, -mu2, rho)
    p01 = pnorm(0, mean = mu1) - p00
    p10 = pnorm(0, mean = mu2) - p00
    p11 = 1 - p00 - p01 - p10
    s00 = sum(response_sub == 0 & activity_sub == 0)
    s01 = sum(response_sub == 0 & activity_sub == 1)
    s10 = sum(response_sub == 1 & activity_sub == 0)
    s11 = sum(response_sub == 1 & activity_sub == 1)
    rst = rst + s00 * mean(log(p00), na.rm = T) + s01 * mean(log(p01), na.rm = T) + s10 * mean(log(p10), na.rm = T) + s11 * mean(log(p11), na.rm = T)
  }
  
  ## Cluster 3
  ind = which(group == 3)
  N = length(ind)
  if (N >= 1) {
    activity_sub = activity[ind, , drop = FALSE]
    response_sub = response[ind, , drop = FALSE]
    dat = list(response = response_sub,
               activity = activity_sub,
               N = N,
               ninter = ninter,
               cutoff = cutoff)
    thismodel = try(jags.model(file = "bugs/sas_binary/c3.txt", 
                               data = dat, 
                               inits = list(Z = array(c(dat$response, dat$activity), dim = c(dat$N, dat$ninter, 2)),
                                            mu1 = 1,
                                            mu2 = 1,
                                            rho = 0.5),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE) 
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("mu1", "mu2", "rho"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    mu1 = as.numeric(res.bugs$mu1)
    mu2 = as.numeric(res.bugs$mu2)
    rho = as.numeric(res.bugs$rho)
    mu1_rec[ind, ] = matrix(rep(mu1, N), nrow = N, byrow = TRUE)
    mu2_rec[ind, ] = matrix(rep(mu2, N), nrow = N, byrow = TRUE)
    p00 = pbivnorm(-mu1, -mu2, rho)
    p01 = pnorm(0, mean = mu1) - p00
    p10 = pnorm(0, mean = mu2) - p00
    p11 = 1 - p00 - p01 - p10
    s00 = sum(response_sub == 0 & activity_sub == 0)
    s01 = sum(response_sub == 0 & activity_sub == 1)
    s10 = sum(response_sub == 1 & activity_sub == 0)
    s11 = sum(response_sub == 1 & activity_sub == 1)
    rst = rst + s00 * mean(log(p00), na.rm = T) + s01 * mean(log(p01), na.rm = T) + s10 * mean(log(p10), na.rm = T) + s11 * mean(log(p11), na.rm = T)
  }
  
  return (list("factor" = rst, "mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec))
}





