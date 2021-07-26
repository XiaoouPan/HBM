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

## This function only works when we have 4 groups
getTrans = function(cluster) {
  trans = matrix(0, 41, 41)
  for (i in 1:41) {
    cur = cluster[i, ]
    uni = unique(cur)
    if (length(uni) == 4) {
      for (j in 1:41) {
        if (length(unique(cluster[j, ])) == 1) {
          trans[i, j] = 1 / 8
        } else if (length(unique(cluster[j, ])) == 3) {
          trans[i, j] = 1 / 24
        }
      }
    } else if (length(uni) == 1) {
      for (j in 1:41) {
        if (length(unique(cluster[j, ])) == 4) {
          trans[i, j] = 1 / 2
        } else if (sum(cluster[j, ] == uni) == 3) {
          trans[i, j] = 1 / 6
        }
      }
    } else if (length(uni) == 2 & as.numeric(table(cur))[1] == 2) {
      for (j in 1:41) {
        if (sum(cluster[j, ] == cur) == 3) {
          trans[i, j] = 1 / 4
        }
      } 
    } else if (length(uni) == 2 & max(table(cur)) == 3) {
      for (j in 1:41) {
        if (sum(cluster[j, ] == cur) == 3 & length(unique(cluster[j, ])) <= 2) {
          trans[i, j] = 1 / 6
        } else if (sum(cluster[j, ] == cur) == 3 & length(unique(cluster[j, ])) == 3) {
          trans[i, j] = 1 / 4
        }
      }
    } else if (length(uni) == 3) {
      for (j in 1:41) {
        if (sum(cluster[j, ] == cur) == 3 & length(unique(cluster[j, ])) == 3) {
          trans[i, j] = 1 / 4
        } else if (length(unique(cluster[j, ])) == 4) {
          trans[i, j] = 1 / 2
        }
      }
    }
  } 
  return (trans)
}

#### MCMC Sampling and calculate likelihood for the final stage
post_sas = function(response, activity, N, ninter, cluster, n.adapt = 1000, n.burn = 1000, n.iter = 5000) {
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter,
             cluster = cluster,
             cluster_prob = rep(1 / dim(cluster)[1], dim(cluster)[1]))
  thismodel = try(jags.model(file = "bugs/sas_binary/sas.txt", 
                             data = dat, 
                             inits = list(Z = array(c(dat$response, dat$activity), dim = c(dat$N, dat$ninter, 2)),
                                          nor1 = rep(0, N),
                                          nor2 = rep(0, N),
                                          rho = 0.5,
                                          cluster_ind = 16),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "rho", "cluster_ind"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = matrix(res.bugs$mu1, nrow = N, ncol = n.iter)
  mu2_rec = matrix(res.bugs$mu2, nrow = N, ncol = n.iter)
  #rho = as.numeric(res.bugs$rho)
  #cluster_ind = as.numeric(res.bugs$cluster_ind)

  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec))
}
