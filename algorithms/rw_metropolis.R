library(dplyr)
library(tidyr)

source("./algorithms/constants.R")
source("./algorithms/data_sim.R")
source("./algorithms/posterior_calculations.R")
#' Run Random Walk Metropolis Metropolis-Hastings Algorithm
#' Proposal dist q(x, y) for x, y in R_{d} at time t: x_{t-1} + epsilon
#' Epsilon ~ N(0, proposal_sigma)
#' @param data matrix with dim = n x d
#' @param mu numeric vector, component means for mixture distribution
#' @param sigma numeric vector, component variances for mixture distribution
#' @param mixture_probs numeric vector, mixture probabilities for mixture distribution
#' @param proposal_sigma float, variance for multivariate normal proposal distribution
#' @param n_iters int, number of steps in the chain
#' @param seed int
#' 
#' @return matrix of steps in the chain
rw_metropolis <- function(data, priors, n_components, proposal_alpha, proposal_sd, n_iters, seed = NA) {
  if (n_iters < 2) {
    rlang::abort("Must run at least 2 iterations in the chain")
  }
  
  incorrect_priors <- setdiff(names(priors),
                              c("alpha", "k_0", "mu_0", "S_0", "v_0"))
  if (length(incorrect_priors) > 0) {
    msg <- paste0("Priors must have names \"alpha\", \"k_0\", \"mu_0\", \"S_0\", ",
                  "and \"v_0\"")
    rlang::abort(msg)
  }

  if (is.na(seed)) {
    seed <- round(runif(1) * 1e7)
  }

  message(paste0("Running Random Walk Metropolis chain with ", n_iters, " steps and seed: ", seed))
  set.seed(seed)

  # initial param values
  p_chain <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
  mu_chain <- rnorm(n_components, sd = proposal_sd)
  sigma_chain <- list(matrix(rWishart(1, n_components, diag(n_components)), ncol = n_components))
  
  # final pre-chain prep
  current_args <- priors
  current_args <- append(current_args, list(
    "n_components" = n_components,
    "data" = data,
    "p" = p_chain,
    "mu" = mu_chain,
    "sigma" = sigma_chain[[1]]))
  
  for (i in 2:n_iters) {
    if (i %% round(n_iters / 5) == 0) {
      message(paste0("Proposing step ", i))
    }
    # proposals
    # symmetric Dirichlet with concentration parameter alpha for mixture probs p
    p_proposal <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
    # symmetric Gaussian with mean 0 and variance proposal_sd^2 * I
    mu_proposal <- rnorm(n_components, mean = matrix(mu_chain, nrow = i - 1)[i - 1,],
                         sd = proposal_sd)
    # Wishart with scale matrix = I and df = n_components
    sigma_proposal <- matrix(rWishart(1, n_components, diag(n_components)), ncol = n_components)
    
    proposal_args <- current_args
    proposal_args[c("p", "mu", "sigma")] <- list(
      "p" = p_proposal,
      "mu" = mu_proposal,
      "sigma" = sigma_proposal
    )
    
    # log hastings ratio log((h(x)/c) / (h(y)/c)) = log(h(x)) - log(h(y))
    h_current <- suppressMessages(do.call(calculate_mixture_posterior, current_args))
    h_proposal <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
    log_hastings_ratio <- h_proposal - h_current
    
    acceptance_prob <- exp(min(0, log_hastings_ratio))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      p_chain <- rbind(p_chain, p_proposal, deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_proposal, deparse.level = 0)
      sigma_chain <- append(sigma_chain, sigma_proposal)
      current_args <- proposal_args
    } else {
      p_chain <- rbind(p_chain, p_chain[i - 1,], deparse.level = 0)
      mu_chain <- rbind(mu_chain, matrix(mu_chain, nrow = i - 1)[i - 1,], deparse.level = 0)
      sigma_chain <- append(sigma_chain, sigma_chain[[i - 1]])
    }
  }
  return(list(
    "p_chain" = p_chain,
    "mu_chain" = mu_chain,
    "sigma_chain" = sigma_chain))
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