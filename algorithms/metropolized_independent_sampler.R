library(gtools)
library(dplyr)
library(tidyr)
library(rlang)

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
                                                 # proposal_alpha,
                                                 proposal_param,
                                                 n_iters,
                                                 true_sigma,
                                                 seed = NA) {
  if (n_iters < 2) {
    rlang::abort("Must run at least 2 iterations in the chain")
  }
  
  # incorrect_priors <- setdiff(names(priors),
  #                             c("alpha", "k_0", "mu_0", "S_0", "v_0"))
  # if (length(incorrect_priors) > 0) {
  #   msg <- paste0("Priors must have names \"alpha\", \"k_0\", \"mu_0\", \"S_0\", ",
  #                 "and \"v_0\"")
  #   rlang::abort(msg)
  # }
  incorrect_priors <- setdiff(names(priors),
                              c("alpha", "mu_0", "tau_2"))
  if (length(incorrect_priors) > 0) {
    rlang::abort("Priors must have names \"alpha\", \"mu_0\", and \"tau_2\"")
  }
  
  if (is.na(seed)) {
    seed <- round(runif(1) * 1e7)
  }
  
  if (is.na(seed)) {
    seed <- round(runif(1) * 1e7)
  }
  
  message(paste0("Running Metropolized Independent Sampler chain with ", n_iters,
                 " steps and seed: ", seed))
  set.seed(seed)
  
  # initial param values
  p_chain <- gtools::rdirichlet(1, rep(4, n_components))
  mu_chain <- matrix(ncol = n_components)
  mu_chain[1,] <- rt(n_components, df = proposal_param)
  z <- sample(1:n_components, length(data), replace = T, prob = p_chain[1,])
  # sigma_chain <- list(matrix(rWishart(1, proposal_df, diag(n_components) / proposal_df),
  #                            ncol = n_components))
  
  # final pre-chain prep
  # current_args <- priors
  # current_args <- append(current_args, list(
  #   "n_components" = n_components,
  #   "data" = data,
  #   "p" = p_chain,
  #   "mu" = mu_chain[1,],
  #   "sigma" = sigma_chain[[1]]))
  
  acceptance_probs <- c()
  acceptances <- c()
  for (i in 2:n_iters) {
    if (i %% round(n_iters / 5) == 0) {
      message(paste0("Proposing step ", i))
    }
    # proposals
    # p_proposal <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
    p_chain <- rbind(p_chain, calculate_gibbs_p_updates(z, p_chain[i - 1,], n_components, priors$alpha))
    z <- calculate_gibbs_z_updates(data, p_chain[i,], mu_chain[i - 1,], true_sigma, n_components)
    if (length(unique(z)) == 1) {
      z[sample(1:length(z), 1)] <- sample(dplyr::setdiff(1:n_components, z), 1)
    }
    # sigma_proposal <- matrix(rWishart(1, proposal_df, diag(n_components) / proposal_df),
    #                          ncol = n_components)

    # proposal_args <- current_args
    # proposal_args[c("p", "mu", "sigma")] <- list(
    #   "p" = p_proposal,
    #   "mu" = mu_proposal,
    #   "sigma" = sigma_proposal
    # )
    
    # log hastings ratio log(h(y) * g(x) / (h(x) * g(y))) = log(h(y)) + log(g(x)) - log(h(x)) - log(g(y))
    # h_x <- suppressMessages(do.call(calculate_mixture_posterior, current_args))
    # g_x <- as.vector(calculate_log_mis_proposal_prob(
    #   p = p_chain[i-1,],
    #   mu = mu_chain[i-1,],
    #   sigma = sigma_chain[[i - 1]],
    #   dir_parameter = rep(proposal_alpha, n_components),
    #   mvn_mu_parameter = rep(0, n_components),
    #   mvr_sigma_parameter = diag(n_components) * proposal_sd^2,
    #   wish_df_parameter = proposal_df,
    #   wish_scale_mat_parameter = diag(n_components) / proposal_df))
    # h_y <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
    # g_y <- as.vector(calculate_log_mis_proposal_prob(
    #   p = p_proposal,
    #   mu = mu_proposal,
    #   sigma = sigma_proposal,
    #   dir_parameter = rep(proposal_alpha, n_components),
    #   mvn_mu_parameter = rep(0, n_components),
    #   mvr_sigma_parameter = diag(n_components) * proposal_sd^2,
    #   wish_df_parameter = proposal_df,
    #   wish_scale_mat_parameter = diag(n_components) / proposal_df))
    h_x <- max(calculate_mh_gibbs_mean_posterior(
      data, p_chain[i,], z, mu_chain[i-1,], true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    q_x <- as.vector(calculate_log_mis_proposal_prob(mu_chain[i-1,], proposal_param))
    
    mu_proposal <- rt(n_components, df = proposal_param)
    h_y <- max(calculate_mh_gibbs_mean_posterior(
      data, p_chain[i,], z, mu_proposal, true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    q_y <- as.vector(calculate_log_mis_proposal_prob(mu_proposal, proposal_param))
    
    log_hastings_ratio <- h_y + q_x - h_x - q_y

    acceptance_prob <- min(1, exp(log_hastings_ratio))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      # p_chain <- rbind(p_chain, p_proposal, deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_proposal, deparse.level = 0)
      # sigma_chain[[i]] <- sigma_proposal
      # current_args <- proposal_args
      # q_current <- q_proposal
    } else {
      # p_chain <- rbind(p_chain, p_chain[i - 1,], deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_chain[i - 1,], deparse.level = 0)
      # sigma_chain[[i]] <- sigma_chain[[i - 1]]
    }
    acceptances <- c(acceptances, accept)
  }
  return(list(
    "p_chain" = p_chain,
    "mu_chain" = mu_chain,
    "acceptance_probs" = acceptance_probs,
    # "sigma_chain" = sigma_chain,
    "acceptances" = acceptances))
}

#' Calculate g(x) or g(y) for the Metropolized Independent Sampler
calculate_log_mis_proposal_prob <- function(#p, mu, sigma,
                                            #dir_parameter, mvn_mu_parameter, mvr_sigma_parameter,
                                            #wish_df_parameter, wish_scale_mat_parameter) {
                                            mu, t_df_parameter) {
  return(
    # log(gtools::ddirichlet(p, dir_parameter)) +
    # log(dMultivariateNormal(mu, mvn_mu_parameter, mvr_sigma_parameter)) +
    # log(dWishart(sigma, wish_df_parameter, wish_scale_mat_parameter))
    sum(log(dt(mu, t_df_parameter)))
  )
}
