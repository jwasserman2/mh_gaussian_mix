library(dplyr)
library(tidyr)

source("./algorithms/constants.R")
source("./algorithms/data_sim.R")
source("./algorithms/posterior_calculations.R")
#' Run Metropolized Independent Sampler Metropolis Metropolis-Hastings Algorithm
#' Proposal dist q(x, y) indepdendent of current position x
#' @param data matrix with dim = n x d
#' @param mu numeric vector, component means for mixture distribution
#' @param sigma numeric vector, component variances for mixture distribution
#' @param mixture_probs numeric vector, mixture probabilities for mixture distribution
#' @param proposal_sigma float, variance for multivariate normal proposal distribution
#' @param n_iters int, number of steps in the chain
#' @param seed int
#' 
#' @return matrix of steps in the chain
run_metropolized_independent_sampler <- function(data,
                                                 priors,
                                                 n_components,
                                                 proposal_alpha,
                                                 proposal_sd,
                                                 proposal_df,
                                                 n_iters,
                                                 seed = NA) {
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
  
  message(paste0("Running Metropolized Independent Sampler chain with ", n_iters,
                 " steps and seed: ", seed))
  set.seed(seed)
  
  # initial param values
  p_chain <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
  mu_chain <- matrix(ncol = n_components)
  mu_chain[1,] <- rnorm(n_components, sd = proposal_sd)
  sigma_chain <- list(matrix(rWishart(1, proposal_df, diag(n_components) / proposal_df),
                             ncol = n_components))
  
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
    # Wishart with scale matrix = I and df = n_components
    sigma_proposal <- matrix(rWishart(1, proposal_df, diag(n_components) / proposal_df),
                             ncol = n_components)
    
    proposal_args <- current_args
    proposal_args[c("p", "mu", "sigma")] <- list(
      "p" = p_proposal,
      "mu" = mu_proposal,
      "sigma" = sigma_proposal
    )
    
    # log hastings ratio log(h(y) * g(x) / (h(x) * g(y))) = log(h(y)) + log(g(x)) - log(h(x)) - log(g(y))
    h_x <- suppressMessages(do.call(calculate_mixture_posterior, current_args))
    message(h_x)
    g_x <- as.vector(calculate_log_mis_proposal_prob(
      p = p_chain[i-1,],
      mu = mu_chain[i-1,],
      sigma = sigma_chain[[i - 1]],
      dir_parameter = rep(proposal_alpha, n_components),
      mvn_mu_parameter = mu_chain[i - 1,],
      mvr_sigma_parameter = diag(n_components) * proposal_sd^2,
      wish_df_parameter = proposal_df,
      wish_scale_mat_parameter = diag(n_components) / proposal_df))
    message(g_x)
    h_y <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
    message(h_y)
    g_y <- as.vector(calculate_log_mis_proposal_prob(
      p = p_proposal,
      mu = mu_proposal,
      sigma = sigma_proposal,
      dir_parameter = rep(proposal_alpha, n_components),
      mvn_mu_parameter = mu_chain[i - 1,],
      mvr_sigma_parameter = diag(n_components) * proposal_sd^2,
      wish_df_parameter = proposal_df,
      wish_scale_mat_parameter = diag(n_components) / proposal_df))
    message(g_y)
    log_hastings_ratio <- h_y + g_x - h_x - g_y
    
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

#' Calculate g(x) or g(y) for the Metropolized Independent Sampler
calculate_log_mis_proposal_prob <- function(p, mu, sigma,
                                            dir_parameter, mvn_mu_parameter, mvr_sigma_parameter,
                                            wish_df_parameter, wish_scale_mat_parameter) {
  return(
    log(pDirichlet(p, dir_parameter)) +
    log(pMultivariateNormal(mu, mvn_mu_parameter, mvr_sigma_parameter)) +
    log(pWishart(sigma, wish_df_parameter, wish_scale_mat_parameter))
  )
}
