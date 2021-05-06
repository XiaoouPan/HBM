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
