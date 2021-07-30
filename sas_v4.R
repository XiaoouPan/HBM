#### Spike and slab code: pairwise difference   ####

s1_acti_sas = function(activity, N, n1, mu2_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  r1 = sum(activity)
  dat = list(activity = r1,
             N = N,
             ninter = n1,
             mu2_h0 = mu2_h0)
  thismodel = try(jags.model(file = "bugs/sas_binary/sas_s1_acti.txt", 
                             data = dat, 
                             inits = list(diff2 = 1,
                                          mu21 = 0,
                                          ss2 = rep(0, N - 1)),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu21", "mu22", "mu23", "mu24"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu2_rec = rbind(as.numeric(res.bugs$mu21), as.numeric(res.bugs$mu22), as.numeric(res.bugs$mu23), as.numeric(res.bugs$mu24))
  return (list("mu2_rec" = mu2_rec))
}


post_sas = function(response, activity, N, ninter, mu1_h0, mu2_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter,
             mu1_h0 = mu1_h0, 
             mu2_h0 = mu2_h0)
  Z = array(0, dim = c(N, ninter, 2))
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1_h0[i], mu2_h0[i]), matrix(c(1, 0.5, 0.5, 1), 2, 2))
    for (j in 1:ninter) {
      Z[i, j, 1] = ifelse(dat$response[i, j] == 1, abs(Z[i, j, 1]), -abs(Z[i, j, 1]))
      Z[i, j, 2] = ifelse(dat$activity[i, j] == 1, abs(Z[i, j, 2]), -abs(Z[i, j, 2]))
    }
  }
  thismodel = try(jags.model(file = "bugs/sas_binary/sas_v4.txt", 
                             data = dat, 
                             inits = list(Z = Z,
                                          diff1 = 1,
                                          diff2 = 1,
                                          mu11 = 0,
                                          mu21 = 0,
                                          ss1 = rep(0, N - 1),
                                          ss2 = rep(0, N - 1),
                                          rho = 0.5),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu11", "mu12", "mu13", "mu14", "mu21", "mu22", "mu23", "mu24", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = rbind(as.numeric(res.bugs$mu11), as.numeric(res.bugs$mu12), as.numeric(res.bugs$mu13), as.numeric(res.bugs$mu14))
  mu2_rec = rbind(as.numeric(res.bugs$mu21), as.numeric(res.bugs$mu22), as.numeric(res.bugs$mu23), as.numeric(res.bugs$mu24))
  rho = as.numeric(res.bugs$rho)
  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec, "rho" = rho))
}
