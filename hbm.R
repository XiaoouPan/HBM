library(MASS)
library(coda)
library(rjags) 
library(gtools)

rm(list = ls())

GetSC <- function(N, p0, mu0, mu, prob, ninter) {
  # N: number of subgroups
  all = rep(ninter, N) # number of subjects in each subgroup during first stage
  # responses of those in the first stage
  response = rbinom(n = N, size = all, prob = prob)
  cutoff = log(p0 / (1 - p0)) #cutoff for the interim step
  cutoff2 = mu0
  activity = rnorm(n = N, mean = mu, sd = rep(1, N))
  return(
    list(
      response = response,
      activity = activity,
      all = all,
      N = N,
      cutoff = cutoff,
      cutoff2 = cutoff2,
      true_prob = prob,
      true_mu = mu
    )
  )
}

posterior_simu2 <- function (dat, iter = 1000) {
  thismodel <-
    try(jags.model(
      file = "basket.txt",
      data = dat,
      inits = list( 
        theta = rep(0, dat$N),
        mu = rep(0, dat$N),
        mumix = c(-10, -10, 0),
        muprec = c(1, 1, 1),
        mumix2 = c(0, 10, 5),
        muprec2 = c(1, 1, 1)
      ),
      n.adapt = iter
    ),
    silent = T)
  ;
  res.bugs <-
    try(jags.samples(
      thismodel,
      variable.names = c('prob', 'theta', 'mu', 'mumix', 'muprec', 'mumix2', 'muprec2'),
      n.iter = iter
    ),
    silent = T)
  ;
  if (length(names(res.bugs)) == 0) {
    return(res.bugs)
  }
  return(list(
    prob = matrix(res.bugs$prob, nrow = dat$N),
    theta = matrix(res.bugs$theta, nrow = dat$N),
    mu = matrix(res.bugs$mu, nrow = dat$N),
    mumix = matrix(res.bugs$mumix, nrow = 3),
    muprec = matrix(res.bugs$muprec, nrow = 3),
    mumix2 = matrix(res.bugs$mumix2, nrow = 3),
    muprec2 = matrix(res.bugs$muprec2, nrow = 3)
  ))
}

summary_posterior2 <- function (dataVal, mcmcVal) {
  #value
  response <- dataVal$response
  activity = dataVal$activity
  all <- dataVal$all
  
  N <- dataVal$N
  #number of patients
  C <- dataVal$C
  #number of clustering
  group <- dataVal$group
  
  #parm
  prob <- mcmcVal$prob
  theta <- mcmcVal$theta
  mu = mcmcVal$mu
  mumix <- mcmcVal$mumix
  muprec <- mcmcVal$muprec
  mumix2 = mcmcVal$mumix2
  muprec2 = mcmcVal$muprec2
  
  #hyper parm
  # # gamma <- 1/sqrt(0.000001);
  # # alpha <- 0.001; beta <- 0.001;
  #bayes factor
  #summarize over all simulations
  res1 <- 0
  res2 <- 0
  res3 = 0
  
  for (n in 1:N) {
    p1 <-
      dbinom(
        x = response[n],
        size = all[n],
        prob = prob[n, ],
        log = T
      )
    
    p2 <-
      dnorm(
        x = theta[n, ],
        mean = mumix[group[n], ],
        sd = 1 / sqrt(muprec[group[n], ]),
        log = T
      )
    
    p3 = dnorm(x = mu[n, ], mean = mumix2[group[n], ], sd = 1 / sqrt(muprec2[group[n], ]), log = T)
    
    res1 <- res1 + mean(p1, na.rm = T)
    res2 <- res2 + mean(p2, na.rm = T)
    res3 = res3 + mean(p3, na.rm = T)
    
  }
  return (res1 + res2 + res3)
}


ninter = 10
N = 4
C = 3
M = 20
post_prob_all = post_acti_all = matrix(0, N, M)
all_cluster <- permutations(n = C, r = N, repeats.allowed = T)

p0 = 0.3
prob = c(0.45, 0.45, 0.45, 0.45)
mu0 = 3
mu = c(1, 1, 5, 5)
SC = GetSC(N, p0, mu0, mu, prob, ninter)
true_prob = SC$true_prob
true_mu = SC$true_mu

for (m in 1:M) {
  set.seed(m)
  #interim analysis
  #initialize recorder for bayes factor of different clustering
  bayes_cluster2 <- NULL
  #to record the cohort clustering
  prob_rec <- NULL 
  mu_rec <- NULL 
  
  for (i in 1:dim(all_cluster)[1]) { # for every cluster permutation for the groups
    group <- all_cluster[i, ] # select a permutation of cluster assignments to each subgroupS
    #cluster with hypothesis
    dat <-
      list(
        response = SC$response,
        activity = SC$activity,
        all = SC$all,
        N = SC$N,
        C = C,
        group = group,
        cutoff = SC$cutoff,
        cutoff2 = SC$cutoff2
      )
    this_posterior <- posterior_simu2(dat) # generate draws from posterior using rjags
    # break the loop if posterior can't be generated
    if (length(names(this_posterior)) == 0) {
      break
    }
    mcmcVal = list(
      prob = this_posterior$prob,
      theta = this_posterior$theta,
      mu = this_posterior$mu,
      mumix = this_posterior$mumix,
      muprec = this_posterior$muprec,
      mumix2 = this_posterior$mumix2,
      muprec2 = this_posterior$muprec2
    )
    prob_rec <- rbind(prob_rec, apply(this_posterior$prob, 1, mean))
    mu_rec = rbind(mu_rec, apply(this_posterior$mu, 1, mean))
    
    
    # Calculate the Bayes Factors for the interim analysis cluster permutations
    res <- summary_posterior2(dataVal = dat, mcmcVal = mcmcVal)
    bayes_cluster2 <- c(bayes_cluster2, res)# this the result vector of BF after iterating thru every permutation
  }
  
  index = which.max(bayes_cluster2)
  post_prob_all[, m]  = prob_rec[index, ]
  post_acti_all[, m] = mu_rec[index, ]
}

max(abs(rowMeans(post_prob_all) - true_prob))
mean((rowMeans(post_acti_all) - true_mu)^2)

mean((post_prob_all > p0 | post_acti_all > mu0) == (true_prob > p0 | true_mu > mu0))
