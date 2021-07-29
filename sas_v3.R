#### Updated spike and slab codes with only 15 clusters
#### All possibilities of clustering
getCluster = function(num = 4) {
  rst = rbind(c(1, 1, 1, 1),
              c(1, 2, 2, 2),
              c(1, 2, 1, 1),
              c(1, 1, 3, 1),
              c(1, 1, 1, 4),
              c(1, 1, 3, 3),
              c(1, 2, 1, 2),
              c(1, 2, 2, 1),
              c(1, 2, 3, 3),
              c(1, 2, 3, 2),
              c(1, 2, 2, 4),
              c(1, 2, 3, 1),
              c(1, 2, 1, 4),
              c(1, 1, 3, 4),
              c(1, 2, 3, 4))
  return (as.matrix(rst))
}

## This function only works when we have 4 groups
getTrans = function(cluster) {
  trans = matrix(0, 15, 15)
  trans[1, 15] = 1 / 2
  for (j in 2:5) {
    trans[1, j] = 1 / 8
  } 
  for (i in 2:5) {
    trans[i, 1] = 1 / 2
    for (j in 6:8) {
      trans[i, j] = 1 / 6
    }
  }
  for (i in 6:8) {
    for (j in 2:5) {
      trans[i, j] = 1 / 8
    }
    for (j in 9:14) {
      trans[i, j] = 1 / 12
    }
  }
  for (i in 9:14) {
    for (j in 6:8) {
      trans[i, j] = 1 / 6
    }
    trans[i, 15] = 1 / 2
  }
  trans[15, 1] = 1 / 2
  for (j in 9:14) {
    trans[15, j] = 1 / 12
  }
  return (trans)
}

#### MCMC Sampling and calculate likelihood for the final stage
post_sas = function(response, activity, N, ninter, cluster, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter,
             cluster = cluster,
             cluster_prob = rep(1 / 15, 15))
  Z = array(0, dim = c(N, ninter, 2))
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2, 2))
    for (j in 1:ninter) {
      Z[i, j, 1] = ifelse(dat$response[i, j] == 1, abs(Z[i, j, 1]), -abs(Z[i, j, 1]))
      Z[i, j, 2] = ifelse(dat$activity[i, j] == 1, abs(Z[i, j, 2]), -abs(Z[i, j, 2]))
    }
  }
  thismodel = try(jags.model(file = "bugs/sas_binary/sas.txt", 
                             data = dat, 
                             inits = list(Z = Z,
                                          nor1 = rep(0, N),
                                          nor2 = rep(0, N),
                                          rho = 0.5,
                                          cluster_ind = 15),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "rho", "cluster_ind"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = matrix(res.bugs$mu1, nrow = N, ncol = n.iter)
  mu2_rec = matrix(res.bugs$mu2, nrow = N, ncol = n.iter)
  rho = as.numeric(res.bugs$rho)
  cluster_ind = as.numeric(res.bugs$cluster_ind)
  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec))
}
