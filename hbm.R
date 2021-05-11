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

posterior_simu2 <- function (dat, iter = 1000) {
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

summary_posterior2 <- function (dataVal, mcmcVal) {
  #value
  response <- dataVal$response
  all <- dataVal$all
  
  N <- dataVal$N
  #number of patients
  C <- dataVal$C
  #number of clustering (3 default)
  group <- dataVal$group
  
  #parm
  prob <- mcmcVal$prob
  theta <- mcmcVal$theta
  
  mumix <- mcmcVal$mumix
  muprec <- mcmcVal$muprec
  
  #hyper parm
  # # gamma <- 1/sqrt(0.000001);
  # # alpha <- 0.001; beta <- 0.001;
  #bayes factor
  #summarize over all simulations
  res1 <- 0
  res2 <- 0
  
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
    
    res1 <- res1 + mean(p1, na.rm = T)
    
    res2 <- res2 + mean(p2, na.rm = T)
    
  }
  res <- res1 + res2
  
  #print(res1);print(res2);print(res3);print(res4);
  #only use the first part as bayes factor #return(res) as all parts
  return(res)
}


nsub = 20
N = 2
p0 = 0.15
p1 = 0.45
e = e1 = e2 = 0.05
k = 4
prob = c(0.15, 0.15, 0.15, 0.15, 0.15)
ninter = 10

SC <- GetSC(nsub, N, p0, p1, e, e1, e2, prob, ninter)
true_prob  <- SC$true_prob
N <- SC$N # This is the number of subgroups

obs_rr_init_orig <- SC$response/SC$all
obs_rr_final_orig <- (SC$response2 - SC$response)/(SC$all2 - SC$all)
obs_rr_total_orig <- (SC$response2)/(SC$all2)

list_discard <- rep(FALSE, N)

#interim analysis
#initialize recorder for bayes factor of different clustering
bayes_cluster2 <- NULL
#to record the bayes factor for each cohort clustering
group_rec <- NULL
#to record the cohort clustering
prob_rec <- NULL 
prob_rec_up <- NULL
prob_rec_down <- NULL
#response rate for all clustering models;
response_rate <- rep(NA, N)
response_rate_up <- rep(NA, N)
response_rate_down <- rep(NA, N)
#
post_prob_p0 <- rep(NA, N)
post_prob_p0_p1 <- rep(NA, N)
post_prob_p1 <- rep(NA, N)

all_cluster <- permutations(n = 2,
                            r = SC$N,
                            repeats.allowed = T)
for (i in 1:dim(all_cluster)[1]) { # for every cluster permutation for the groups
  group <- all_cluster[i, ] # select a permutation of cluster assignments to each subgroupS
  
  group_rec <- rbind(group_rec, group)
  
  #cluster with hypothesis
  dat <-
    list(
      response = SC$response,
      all = SC$all,
      N = SC$N,
      C = 2,
      group = group,
      cutoff = SC$cutoff
    )
  
  this_posterior <- posterior_simu2(dat) # generate draws from posterior using rjags
  
  # break the loop if posterior can't be generated
  if (length(names(this_posterior)) == 0) {
    break
  }
  mcmcVal = list(
    prob = this_posterior$prob,
    theta = this_posterior$theta,
    mumix = this_posterior$mumix,
    muprec = this_posterior$muprec
  )
  
  prob_rec <- rbind(prob_rec, apply(this_posterior$prob, 1, mean))
  
  prob_rec_up <-
    rbind(prob_rec_up,
          apply(this_posterior$prob, 1, quantile, 0.975))
  
  prob_rec_down <-
    rbind(prob_rec_down,
          apply(this_posterior$prob, 1, quantile, 0.025))
  
  # Calculate the Bayes Factors for the interim analysis cluster permutations
  res <- summary_posterior2(dataVal = dat, mcmcVal = mcmcVal)
  
  bayes_cluster2 <- c(bayes_cluster2, res)# this the result vector of BF after iterating thru every permutation
  
}

best_cluster2 <- group_rec[which.max(bayes_cluster2), ]


this_response <- prob_rec[which.max(bayes_cluster2), ]

this_response_up <- prob_rec_up[which.max(bayes_cluster2), ]

this_response_down <- prob_rec_down[which.max(bayes_cluster2), ]

#get the cluster for unsorted original prob
best_cluster2 <- best_cluster2

list_discard[best_cluster2 == 1] <- TRUE # if in inactive cluster add to discard list

list_keep <- !list_discard


response_rate[list_discard] <- this_response[list_discard]

response_rate_up[list_discard] <- this_response_up[list_discard]

response_rate_down[list_discard] <- this_response_down[list_discard]

response_rate

post_prob_p0 <- as.numeric(response_rate <= p0)
post_prob_p0_p1 <- as.numeric((response_rate > p0) & (response_rate <= p1))
post_prob_p1 <- as.numeric(response_rate > p1) 
