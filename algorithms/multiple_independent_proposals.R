library(dplyr)
library(tidyr)
library(MASS)

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
run_multiple_independent_proposals <- function(data, K, priors, n_components,
                                               # proposal_alpha,
                                               proposal_param,
                                               # max_proposal_sigma,
                                               # proposal_corr_0,
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
  
  message(paste0("Running Multiple Independent Proposals chain with ", n_iters, " steps and seed: ", seed))
  set.seed(seed)
  
  # initial param values
  p_chain <- gtools::rdirichlet(1, rep(4, n_components))
  mu_chain <- matrix(ncol = n_components)
  mu_chain[1,] <- rnorm(n_components, sd = proposal_param)
  z <- sample(1:n_components, length(data), replace = T, prob = p_chain[1,])
  # init_sigma_variances <- runif(n_components, min = 0.001, max = max_proposal_sigma)
  # init_sigma_cov <- runif(1, -min(abs(init_sigma_variances)), min(abs(init_sigma_variances))) / proposal_corr_0
  # sigma_chain <- list(diag(init_sigma_variances - init_sigma_cov) + init_sigma_cov)
  
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
    p_chain <- rbind(p_chain, calculate_gibbs_p_updates(z, p_chain[i - 1,], n_components, priors$alpha))
    z <- calculate_gibbs_z_updates(data, p_chain[i,], mu_chain[i - 1,], true_sigma, n_components)
    if (length(unique(z)) == 1) {
      z[sample(1:length(z), 1)] <- sample(dplyr::setdiff(1:n_components, z), 1)
    }
    # proposals
    # p_proposals <- gtools::rdirichlet(K, rep(proposal_alpha, n_components))
    # sigma_proposals <- rWishart(K, proposal_df, diag(n_components) / proposal_df)
    # proposal_variances <- matrix(runif(K * n_components, min = 0.001, max = max_proposal_sigma),
    #                              nrow = K)
    # cov_bounds <- apply(abs(proposal_variances), 1, min)
    # proposal_covs <- purrr::map_dbl(cov_bounds, ~ runif(1, -.x, .x) / proposal_corr_0)
    # sigma_proposals <- cbind(proposal_variances, proposal_covs, deparse.level = 0) %>%
    #   apply(., 1, function(x) {diag(x[1:(length(x) - 1)] - x[length(x)]) + x[length(x)]},
    #         simplify = F)
    
    mu_proposals <- MASS::mvrnorm(K, mu = mu_chain[i - 1,],
                                  Sigma = diag(n_components) * proposal_param^2)

    log_weights <- c()
    for (k in 1:K) {
      # proposal_args <- current_args
      # proposal_args[c("p", "mu", "sigma")] <- list(
      #   "p" = p_proposals[k,],
      #   "mu" = mu_proposals[k,],
      #   "sigma" = sigma_proposals[[k]]
      # )
      # h_k <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
      h_k <- max(calculate_mh_gibbs_mean_posterior(
        data, p_chain[i,], z, mu_proposals[k,], true_sigma, n_components, priors$mu_0, priors$tau_2),
        -1e7)
      q_x_y <- calculate_log_mip_proposal_prob(mu_proposals[k,], mu_chain[i - 1,], proposal_param)
      l_x_y <- -q_x_y
      log_weights <- c(log_weights, h_k)
    }
    
    # y_idx <- as.vector(rmultinom(1, 1, (-1 / log_weights) / sum(-1 / log_weights)))
    y_idx <- as.vector(rmultinom(1, 1, exp(log_weights) / sum(exp(log_weights))))
    y <- list(# "p" = p_proposals[which.max(y_idx), ],
              "mu" = mu_proposals[which.max(y_idx), ])
              # "sigma" = sigma_proposals[[which.max(y_idx)]])
    
    # p_star_proposals <- rbind(gtools::rdirichlet(K - 1, rep(proposal_alpha, n_components)),
    #                           p_chain[i - 1,])
    # mu_star_proposals <- rbind(matrix(rnorm((K - 1) * n_components, mean = y$mu,
    #                                         sd = proposal_param), nrow = K - 1),
    #                            mu_chain[i - 1,])
    mu_star_proposals <- rbind(MASS::mvrnorm(K - 1, mu = y$mu,
                                             Sigma = diag(n_components) * proposal_param^2),
                               mu_chain[i -1,])
    # sigma_proposals <- rWishart(K, proposal_df, diag(n_components) / proposal_df)
    # proposal_star_variances <- matrix(runif((K - 1) * n_components, min = 0.001,
    #                                         max = max_proposal_sigma), nrow = K - 1)
    # cov_star_bounds <- apply(abs(proposal_star_variances), 1, min)
    # proposal_star_covs <- purrr::map_dbl(cov_star_bounds, ~ runif(1, -.x, .x) / proposal_corr_0)
    # sigma_star_proposals <- cbind(proposal_star_variances, proposal_star_covs, deparse.level = 0) %>%
    #   apply(., 1, function(x) {diag(x[1:(length(x) - 1)] - x[length(x)]) + x[length(x)]},
    #         simplify = F)
    # sigma_star_proposals[[K]] <- sigma_chain[[i - 1]]
    
    # star_proposal_args <- current_args
    # star_proposal_args$mu <- y$mu
    log_star_weights <- c()
    for (k in 1:K) {
      # star_proposal_args[c("p", "mu", "sigma")] <- list(
      #   "p" = p_star_proposals[k,],
      #   "mu" = mu_star_proposals[k,],
      #   "sigma" = sigma_star_proposals[[k]]
      # )
      # h_k <- suppressMessages(do.call(calculate_mixture_posterior, star_proposal_args))
      h_k <- max(calculate_mh_gibbs_mean_posterior(
        data, p_chain[i,], z, mu_star_proposals[k,], true_sigma, n_components, priors$mu_0, priors$tau_2),
        -1e7)
      # q_x_y <- calculate_log_mip_proposal_prob(p_star_proposals[k,], mu_star_proposals[k,],
      #                                          sigma_star_proposals[[k]],
      #                                          rep(proposal_alpha, n_components),
      #                                          y$mu, proposal_param,
      #                                          max_proposal_sigma)
      q_x_y <- calculate_log_mip_proposal_prob(mu_star_proposals[k,], y$mu, proposal_param)
      l_x_y <- -q_x_y
      log_star_weights <- c(log_star_weights, h_k)
    }
    
    # log hastings ratio log(sum(weights(y))) - log(sum(weights(x)))
    log_hastings_ratio <- log(sum(exp(log_weights))) - log(sum(exp(log_star_weights)))
    
    acceptance_prob <- min(1, exp(log_hastings_ratio))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      # p_chain <- rbind(p_chain, y$p, deparse.level = 0)
      mu_chain <- rbind(mu_chain, y$mu, deparse.level = 0)
      # sigma_chain[[i]] <- y$sigma
      # current_args[c("p", "mu", "sigma")] <- y
    } else {
      # p_chain <- rbind(p_chain, p_chain[i - 1,], deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_chain[i - 1,], deparse.level = 0)
      # sigma_chain[[i]] <- sigma_chain[[i - 1]]
    }
    acceptance_probs <- c(acceptance_probs, acceptance_prob)
    acceptances <- c(acceptances, accept)
  }
  return(list(
    "p_chain" = p_chain,
    "mu_chain" = mu_chain,
    "acceptance_probs" = acceptance_probs,
    # "sigma_chain" = sigma_chain,
    "acceptances" = acceptances))
}

#' Calculate g(x) or g(y) for the Multiple Independent Proposals Sampler
calculate_log_mip_proposal_prob <- function(#p, mu, sigma,
                                            # dir_parameter,
                                            mu, n_mu_parameter, n_sigma_parameter) {
                                            # unif_var_parameter) {
  return(
    # log(gtools::ddirichlet(p, dir_parameter)) +
    # sum(log(dnorm(mu, n_mu_parameter, n_sigma_parameter))) +
    # #log(dWishart(sigma, wish_df_parameter, wish_scale_mat_parameter))
    # n_components * log(1 / (unif_var_parameter - 0.001)) +
    # log(1 / (min(diag(sigma)) + min(diag(sigma))))
    sum(log(dnorm(mu, n_mu_parameter, n_sigma_parameter)))
  )
}

