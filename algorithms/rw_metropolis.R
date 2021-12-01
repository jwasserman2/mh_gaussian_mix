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
rw_metropolis <- function(data, priors, n_components, proposal_alpha, proposal_sd,
                          max_proposal_sigma, proposal_corr_0,
                          n_iters, seed = NA) {
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
  mu_chain <- matrix(ncol = n_components)
  mu_chain[1,] <- rnorm(n_components, sd = proposal_sd)
  init_sigma_variances <- runif(n_components, min = 0.001, max = max_proposal_sigma)
  init_sigma_cov <- runif(1, -min(abs(init_sigma_variances)), min(abs(init_sigma_variances))) / proposal_corr_0
  sigma_chain <- list(diag(init_sigma_variances - init_sigma_cov) + init_sigma_cov)
    
  # final pre-chain prep
  current_args <- priors
  current_args <- append(current_args, list(
    "n_components" = n_components,
    "data" = data,
    "p" = p_chain,
    "mu" = mu_chain[1,],
    "sigma" = sigma_chain[[1]]))

  acceptances <- c()
  for (i in 2:n_iters) {
    if (i %% round(n_iters / 5) == 0) {
      message(paste0("Proposing step ", i))
    }
    # proposals
    # symmetric Dirichlet with concentration parameter alpha for mixture probs p
    p_proposal <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
    # symmetric Gaussian with mean 0 and variance proposal_sd^2 * I
    mu_proposal <- rnorm(n_components, mean = mu_chain[i - 1,], sd = proposal_sd)
    # symmetric distributions for drawing variances and covariance independently
    proposal_variances <- runif(n_components, min = 0.001, max = max_proposal_sigma)
    proposal_cov <- runif(1, -min(abs(proposal_variances)), min(abs(proposal_variances))) / proposal_corr_0
    sigma_proposal <- diag(proposal_variances - proposal_cov) + proposal_cov

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
      sigma_chain[[i]] <- sigma_proposal
      current_args <- proposal_args
      q_current <- q_proposal
    } else {
      p_chain <- rbind(p_chain, p_chain[i - 1,], deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_chain[i - 1,], deparse.level = 0)
      sigma_chain[[i]] <- sigma_chain[[i - 1]]
    }
    acceptances <- c(acceptances, accept)
  }
  return(list(
    "p_chain" = p_chain,
    "mu_chain" = mu_chain,
    "sigma_chain" = sigma_chain,
    "acceptances" = acceptances))
}
