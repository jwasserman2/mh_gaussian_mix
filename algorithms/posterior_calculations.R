library(gtools)
library(dplyr)
library(purrr)
library(tidyr)

dInvWishart <- function(sigma, df, scale_mat) {
  p <- ncol(sigma)
  stopifnot(df > p)
  return(
    1 /
      (2^(df * p / 2) * pi^(choose(p, 2) / 2) * det(scale_mat)^(-df / 2) *
         prod(gamma((df + 1 - 1:p) / 2))) *
      det(sigma)^((df - p - 1) / 2) *
      exp(-0.5 * sum(eigen(sigma %*% scale_mat)$values))
  )
}

dMultivariateNormal <- function(x, mu, Sigma) {
  p <- length(x)
  return(
    (2 * pi)^(-p/2) * det(Sigma)^(-1/2) * exp((-1/2) * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  )
}

dWishart <- function(sigma, df, scale_mat) {
  p <- ncol(sigma)
  return(
    1 / 
      (2^(df * p / 2) * pi^(p * (p - 1) / 4) * det(scale_mat)^(df / 2) *
         prod(gamma((df + 1 - 1:p) / 2))) *
      det(sigma)^((df - p - 1) / 2) * exp(-0.5 * sum(eigen(sigma %*% solve(scale_mat))$values))
  )
}

#' Calculate the posterior probability for the parameters
#' @param p numeric vector, proposed mixture probabilities 
#' @param mu numeric vector, proposed component means
#' @param sigma matrix, proposed covariance matrix
#' @param n_components int, number of mixture components
#' @param alpha numeric vector, hyperparameter for the Dirichlet prior for the mixture probabilities
#' @param k_0 float, scale hyperparameter on the variance for component means prior
#' @param mu_0 numeric, means hyperparameter for component means prior
#' @param S_0 matrix, scale matrix hyperparameter for Inverse-Wishart prior on covariance matrix
#' @param v_0 int, degrees of freedom hyperparameter for Inverse-Wishart prior on covariance matrix
#' @param tau_2 int, variance hyperparameter for component means prior
calculate_mixture_posterior <- function(data,
                                        p,
                                        mu,
                                        sigma,
                                        n_components,
                                        alpha,
                                        k_0,
                                        mu_0,
                                        S_0,
                                        v_0,
                                        tau_2) {
  ## PRIORS
  # component means prior: N(mu_0, sigma_0)
  # log_mean_prior_prob <- log(dMultivariateNormal(mu, mu_0, k_0 * sigma))
  # log_mean_prior_prob <- log(prod(dnorm(mu, mu_0, tau_0)))
  message(paste0("Log prior probability for the component means is: ", log_mean_prior_prob))
  # component variances prior: Inv-Wishart(df_0, sigma_0)
  # log_sigma_prior_prob <- log(dInvWishart(sigma, v_0, S_0))
  # message(paste0("Log prior probability for the component variance structure is: ", log_sigma_prior_prob))

  # log(p(zi | p) * p(p))
  z <- sample(1:n_components, length(data), replace = TRUE, prob = p)
  # log_x_sum <- data.frame(
  #   "x" = data,
  #   "p_z" = current_args$p[z],
  #   "mu" = current_args$mu[z],
  #   "sd" = sqrt(diag(current_args$sigma))[z]) %>%
  #   dplyr::mutate(log_p_x = log(dnorm(x, mu, sd) * p_z)) %>%
  #   dplyr::summarize(l = sum(log_p_x)) %>%
  #   dplyr::pull(l)
  p_sum <- data.frame("z" = z) %>%
    dplyr::count(z) %>%
    dplyr::full_join(data.frame("z" = 1:n_components), by = "z") %>%
    dplyr::mutate_at("n", tidyr::replace_na, 0) %>%
    dplyr::arrange(z) %>%
    cbind(., "p" = p, "alpha" = alpha, "mu" = mu) %>%
    dplyr::summarize(l = sum((n + alpha - 1) * log(p) + dnorm(mu, mean = mu_0, sd = sqrt(tau_2)))) %>%
    dplyr::pull(l)
  message(paste0("P sum is :", p_sum))
  # p(xi | zi = j, mu_k, sigma_k)
  log_x_sum <- purrr::pmap_dbl(list(
    "x" = data,
    "mean" = mu[z],
    "sd" = sqrt(diag(sigma)[z])),
    dnorm) %>%
    log() %>%
    sum()
  message(paste0("Log X sum is: ", log_x_sum))

  # return(sum(log_mean_prior_prob, log_sigma_prior_prob, p_sum, log_x_sum))
  return(sum(p_sum, log_x_sum))
}

#' MH/Gibbs Mean Posterior Calculation
calculate_mh_gibbs_mean_posterior <- function(data, p, z, mu, sigma, n_components, mu_0, tau_2) {
  log_posterior_prob_sum <- 0
  for (k in 1:n_components) {
    log_posterior_prob_sum <- log_posterior_prob_sum +
      log(dnorm(mu[k], (mu_0[k] / tau_2 + sum(data[z == k]) / diag(sigma)[k]) /
                           (1 / tau_2 + sum(z == k) / diag(sigma)[k]),
          sqrt(1 / (1 / tau_2 + sum(z == k) / diag(sigma)[k]))))
  }

  return(log_posterior_prob_sum)
}

#' Gibbs sampler for latent assignments Z
#' @return numeric vector of Z draws
calculate_gibbs_z_updates <- function(data, p, mu, sigma, n_components) {
  z <- sample(1:n_components, length(data), replace = TRUE, prob = p)
  p_x_all <- data.frame("x" = data, "z" = z, "p_z" = p[z],  "mu_z" = mu[z],
                        "sd_z" = sqrt(diag(sigma))[z], "p_x_all" = rep(0, length(data)))
  p_x_all <- data.frame("p_x_z_1" = p[1] * dnorm(data, mu[1], sqrt(diag(sigma))[1]))
  for (i in 2:n_components) {
    p_x_all <- p_x_all %>%
      dplyr::mutate(!!sym(paste0("p_x_z_", i)) := p[i] * dnorm(data, mu[i], sqrt(diag(sigma))[i]))
  }

  z_draws <- t(p_x_all) %>%
    rbind(., "p_x_all" = colSums(.)) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate_at(vars(tidyselect::starts_with("p_x_z")), list("multinom" = ~ .x / p_x_all)) %>%
    dplyr::select(tidyselect::contains("multinom")) %>%
    purrr::pmap_dbl(., function(...) {
      return(sample(1:n_components, 1, replace = TRUE, prob = c(...)))
      })
  return(z_draws)
}

#' Gibbs sampler for mixture probabilities pi
#' @return (1, n_components) dimension matrix of mixture probabilities
calculate_gibbs_p_updates <- function(z, p, n_components, alpha) {
  p_update_args <- data.frame("z" = z) %>%
    dplyr::count(z) %>%
    dplyr::full_join(data.frame("z" = 1:n_components), by = "z") %>%
    dplyr::mutate_at("n", tidyr::replace_na, 0) %>%
    dplyr::arrange(z) %>%
    cbind(., "p" = p, "alpha" = alpha)
  return(gtools::rdirichlet(1, p_update_args$n + p_update_args$alpha))
}
  
  