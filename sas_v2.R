#### All possibilities of clustering
getCluster = function(num = 4) {
  rst = NULL
  perm = permutations(n = num, r = num, repeats.allowed = T)
  for (i in 1:dim(perm)[1]) {
    cluster = perm[i, ]
    add = TRUE
    for (j in 1:num) {
      if (cluster[cluster[j]] != cluster[j]) {
        add = FALSE
        break
      }
    }
    if (add) {
      rst = rbind(rst, cluster)
    }
  }
  return (as.matrix(rst))
}


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




#### triCRM: MCMC Sampling and calculate likelihood for the final stage
post_crm = function(outcome, ninter, group, cutoff, cutoff2, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  rst = 0  ## Bayesian factor for this group
  p_c0_rec = p_c1_rec = p_c2_rec = matrix(NA, length(group), n.iter)
  
  ## Cluster 1
  ind = which(group == 1)
  N = length(ind)
  if (N >= 1) {
    outcome_sub = outcome[ind, , drop = FALSE]
    dat = list(outcome = outcome_sub,
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/sas_binary/c1_crm.txt", 
                               data = dat, 
                               inits = list(theta = -2.5,
                                            eta = -2),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("theta", "eta"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    theta = as.numeric(res.bugs$theta)
    eta = as.numeric(res.bugs$eta)
    p_c0 = 1.0 / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c1 = exp(theta) / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c2 = exp(eta) / (1 + exp(eta))
    p_c0_rec[ind, ] = matrix(rep(p_c0, N), nrow = N, byrow = TRUE)
    p_c1_rec[ind, ] = matrix(rep(p_c1, N), nrow = N, byrow = TRUE)
    p_c2_rec[ind, ] = matrix(rep(p_c2, N), nrow = N, byrow = TRUE)
    s_c0 = sum(outcome_sub == 1)
    s_c1 = sum(outcome_sub == 2)
    s_c2 = sum(outcome_sub == 3)
    rst = rst + s_c0 * mean(log(p_c0), na.rm = T) + s_c1 * mean(log(p_c1), na.rm = T) + s_c2 * mean(log(p_c2), na.rm = T)
  }
  
  ## Cluster 2
  ind = which(group == 2)
  N = length(ind)
  if (N >= 1) {
    outcome_sub = outcome[ind, , drop = FALSE]
    dat = list(outcome = outcome_sub,
               N = N,
               ninter = ninter,
               cutoff = cutoff,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/sas_binary/c2_crm.txt", 
                               data = dat, 
                               inits = list(theta = -0.5,
                                            eta = -2),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("theta", "eta"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    theta = as.numeric(res.bugs$theta)
    eta = as.numeric(res.bugs$eta)
    p_c0 = 1.0 / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c1 = exp(theta) / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c2 = exp(eta) / (1 + exp(eta))
    p_c0_rec[ind, ] = matrix(rep(p_c0, N), nrow = N, byrow = TRUE)
    p_c1_rec[ind, ] = matrix(rep(p_c1, N), nrow = N, byrow = TRUE)
    p_c2_rec[ind, ] = matrix(rep(p_c2, N), nrow = N, byrow = TRUE)
    s_c0 = sum(outcome_sub == 1)
    s_c1 = sum(outcome_sub == 2)
    s_c2 = sum(outcome_sub == 3)
    rst = rst + s_c0 * mean(log(p_c0), na.rm = T) + s_c1 * mean(log(p_c1), na.rm = T) + s_c2 * mean(log(p_c2), na.rm = T)
  }
  
  ## Cluster 3
  ind = which(group == 3)
  N = length(ind)
  if (N >= 1) {
    outcome_sub = outcome[ind, , drop = FALSE]
    dat = list(outcome = outcome_sub,
               N = N,
               ninter = ninter,
               cutoff2 = cutoff2)
    thismodel = try(jags.model(file = "bugs/sas_binary/c3_crm.txt", 
                               data = dat, 
                               inits = list(theta = -0.5,
                                            eta = 0),
                               n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
    try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
    res.bugs = try(jags.samples(thismodel, 
                                variable.names = c("theta", "eta"),
                                n.iter = n.iter, progress.bar = "none"), silent = TRUE)
    theta = as.numeric(res.bugs$theta)
    eta = as.numeric(res.bugs$eta)
    p_c0 = 1.0 / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c1 = exp(theta) / ((1 + exp(theta)) * (1 + exp(eta)))
    p_c2 = exp(eta) / (1 + exp(eta))
    p_c0_rec[ind, ] = matrix(rep(p_c0, N), nrow = N, byrow = TRUE)
    p_c1_rec[ind, ] = matrix(rep(p_c1, N), nrow = N, byrow = TRUE)
    p_c2_rec[ind, ] = matrix(rep(p_c2, N), nrow = N, byrow = TRUE)
    s_c0 = sum(outcome_sub == 1)
    s_c1 = sum(outcome_sub == 2)
    s_c2 = sum(outcome_sub == 3)
    rst = rst + s_c0 * mean(log(p_c0), na.rm = T) + s_c1 * mean(log(p_c1), na.rm = T) + s_c2 * mean(log(p_c2), na.rm = T)
  }
  
  return (list("factor" = rst, "p_c0_rec" = p_c0_rec, "p_c1_rec" = p_c1_rec, "p_c2_rec" = p_c2_rec))
}




