library(MASS)
library(coda)
library(rjags) 
library(gtools)

GetSC <- function(nsub, N, p0, p1, e, e1, e2, prob, ninter) {
  # N = 5 # number of subgroups
  all = rep(ninter, N) # number of subjects in each subgroup during first stage
  all2 = rep(nsub, N) # number of subjects in each subgroup total
  
  # true probability of response
  # prob = c(0.15, 0.15, 0.15, 0.15, 0.15) 
  prob = prob
  
  # responses of those in the first stage
  response = rbinom(n = N, size = all, prob = prob)
  # responses of the remaining subjects plus initial responses
  response2 = rbinom(n = N,
                     size = all2 - all,
                     prob = prob) + response
  
  cutoff = log( (p0 + e) / (1 - (p0 + e)) ) #cutoff for the interim step;
  cutoff2 = log( (p0 + e1) / (1 - (p0 + e1)) ) #left cutoff for the final step;
  cutoff3 = log( (p1 - e2) / (1 - (p1 - e2)) ) #right cutoff for the final step;
  return(
    list(
      response = response,
      all = all,
      N = N,
      cutoff = cutoff,
      response2 = response2,
      all2 = all2,
      cutoff2 = cutoff2,
      cutoff3 = cutoff3,
      true_prob = prob
    )
  )
}

posterior_simu2 <- function (dat, iter = 1) {
  thismodel <-
    try(jags.model(
      file = "basket.txt",
      data = dat,
      # initial values: theta_j's at zero
      #                 mumix is length 2, 1st cluster centered at -10 and 2nd cluster at 0 
      #                 muprec is length 2, 1st cluster 1 and 2nd cluster 1
      inits = list( 
        theta = rep(0, dat$N),
        mumix = c(-10, rep(0, dat$C - 1)),
        muprec = rep(1, dat$C)
      ),
      n.adapt = 1000
    ),
    silent = T)
  ;
  res.bugs <-
    try(jags.samples(
      thismodel,
      variable.names = c('prob', 'theta', 'mumix', 'muprec'),
      n.iter = 4000
    ),
    silent = T)
  ;
  if (length(names(res.bugs)) == 0) {
    return(res.bugs)
  }
  return(list(
    prob = matrix(res.bugs$prob, nrow = dat$N),
    theta = matrix(res.bugs$theta, nrow = dat$N),
    mumix = matrix(res.bugs$mumix, nrow = dat$C),
    muprec = matrix(res.bugs$muprec, nrow = dat$C)
  ))
}

