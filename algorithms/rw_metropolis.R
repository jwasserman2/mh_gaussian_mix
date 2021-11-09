source("constants.R")
source("data_sim.R")
#' Run Random Walk Metropolis Metropolis-Hastings Algorithm
#' Proposal dist q(x, y) for x, y in R_{d} at time t: x_{t-1} + epsilon
#' Epsilon ~ N(0, proposal_sigma)
#' @param n_dim int
#' @param mu numeric vector, component means for mixture distribution
#' @param sigma numeric vector, component variances for mixture distribution
#' @param mixture_probs numeric vector, mixture probabilities for mixture distribution
#' @param proposal_sigma float, variance for multivariate normal proposal distribution
#' @param n_iters int, number of steps in the chain
#' @param seed int
#' 
#' @return matrix of steps in the chain
rw_metropolis <- function(n_dim, mu, sigma, mixture_probs, proposal_sigma, n_iters, seed = NA) {
  if (is.na(seed)) {
    seed <- round(runif(1) * 1e7)
  }
  set.seed(seed)
  chain <- runif(n_dim)
  message(paste0("Running Random Walk Metropolis chain with ", n_iters, " steps and seed: ", seed))
  for (i in 1:n_iters) {
    if (i %% 1000 == 0) {
      message(paste0("Proposing step ", i))
    }
    proposal <- chain[i] + rnorm(n_dim, sd = sqrt(proposal_sigma))
    acceptance_prob <- min(1, calculate_gaussian_mixture_prob(proposal, mu, sigma, mixture_probs) /
                             calculate_gaussian_mixture_prob(chain[i], mu, sigma, mixture_probs))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      chain <- rbind(chain, proposal, deparse.level = 0)
    } else {
      chain <- rbind(chain, chain[i], deparse.level = 0)
    }
  }
  return(chain)
}

#' Make a spherically symmetric proposal
#' Propose a step in Random Walk Metropolis chain using multivariate normal with uncorrelated
#' standard normal components
#' @param pos_vec numeric vector
#' 
#' @return numeric vector of pos_vec + proposed steps
propose <- function(pos_vec) {
  n_covs <- length(pos_vec)
  steps <- rnorm(n_covs)
  # magnitude <- runif(1, min = MIN_RW_METROPOLIS_STEP, MAX_RW_METROPOLIS_STEP)
  # iterate_order <- sample(1:n_covs)
  # step_sum <- MAX_RW_METROPOLIS_STEP
  # while(step_sum > MAX_RW_METROPOLIS_STEP - MIN_RW_METROPOLIS_STEP) {
  #   steps <- runif(n_covs - 1, -MAX_RW_METROPOLIS_STEP / 2, MAX_RW_METROPOLIS_STEP / 2)
  #   step_sum <- sum(abs(steps))
  # }
  # steps[n_covs] <- ifelse(rbinom(1, 1, 0.5),
  #                         -(MAX_RW_METROPOLIS_STEP - step_sum),
  #                         MAX_RW_METROPOLIS_STEP - step_sum)
  # return(pos_vec + steps[iterate_order])
  return(pos_vec + steps)
}