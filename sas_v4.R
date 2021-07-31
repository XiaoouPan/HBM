#### Spike and slab code: pairwise difference   ####

s1_acti_sas = function(activity, N, n1, mu2_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  r1 = rowSums(activity)
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


post_sas = function(response, activity, N, ninter, p0, mu1_h0, a0, mu2_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
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
  ss10 = as.numeric(rowMeans(response) > p0 + 0.15)
  ss20 = as.numeric(rowMeans(activity) > a0 + 0.15)
  p10 = 0.25 + 0.5 * as.numeric(rowMeans(response) > p0 + 0.15)
  p20 = 0.25 + 0.5 * as.numeric(rowMeans(activity) > a0 + 0.15)
  thismodel = try(jags.model(file = "bugs/sas_binary/sas_v4.txt", 
                             data = dat, 
                             inits = list(Z = Z,
                                          diff1 = 1,
                                          diff2 = 1,
                                          ss1 = ss10,
                                          ss2 = ss20,
                                          p1 = p10,
                                          p2 = p20,
                                          rho = 0.5),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "p1", "p2", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = matrix(res.bugs$mu1, N, n.iter)
  mu2_rec = matrix(res.bugs$mu2, N, n.iter)
  p1 = matrix(res.bugs$p1, N, n.iter)
  p2 = matrix(res.bugs$p2, N, n.iter)
  rho = as.numeric(res.bugs$rho)
  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec, "p1" = p1, "p2" = p2, "rho" = rho))
}
