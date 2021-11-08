source("constants.R")
source("data_sim.R")
#' Run Random Walk Metropolis Metropolis-Hastings Algorithm
#' Proposal dist q(x, y) for x, y in R_{d} at time t: x_{t-1} + epsilon
#' Epsilon ~ spherically symmetric with variance = sigma^2
#' @param n_dim int
#' @param true_density_args list, must include "gaussian_mixture", "mu", "sigma", and "mixture_probs"
rw_metropolis <- function(n_dim, gaussian_mix, mu, sigma, mixture_probs, n_iters) {
  chain <- runif(n_dim)
  message(paste0("Running Random Walk Metropolis chain with ", n_iters, " steps"))
  for (i in 1:n_iters) {
    if (i %% 1000 == 0) {
      message(paste0("Proposing step ", i))
    }
    proposal <- propose(chain[i])
    acceptance_prob <- min(1, calculate_true_density(proposal, gaussian_mix, mu, sigma, mixture_probs) /
                             calculate_true_density(chain[i], gaussian_mix, mu, sigma, mixture_probs))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      chain <- rbind(chain, proposal, deparse.level = 0)
    } else {
      chain <- rbind(chain, chain[i], deparse.level = 0)
    }
  }
  return(chain)
}

#' Make a spherically symmetric proposal (using multivariate standard normal)
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

mu_vec <- c(-5, 0, 5)
sigma_vec <- rep(2, 3)
mixture_probs <- c(0.45, 0.1, 0.45)
chain <- rw_metropolis(1, TRUE, mu_vec, sigma_vec, mixture_probs, 5000)
ggplot2::ggplot(data.frame("chain" = chain)) +
  ggplot2::geom_histogram(ggplot2::aes(x = chain), bins = 30) +
  ggplot2::geom_line(
    data = data.frame("x" = seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 0.1),
                      "true_density" = purrr::map_dbl(
                        seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 0.1),
                        calculate_gaussian_mixture_prob,
                        mu_vec,
                        sigma_vec,
                        mixture_probs
                      )),
    ggplot2::aes(x = x, y = true_density * 10000)) +
  ggplot2::scale_x_continuous(
    breaks = round(seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 1)))

theta_2_chain <- rw_metropolis(2, TRUE, mu_vec, sigma_vec, mixture_probs, 10000)
ggplot2::ggplot(data.frame("chain" = theta_2_chain[,2])) +
  ggplot2::geom_histogram(ggplot2::aes(x = chain), bins = 30) +
  ggplot2::geom_line(
    data = data.frame("x" = seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 0.1),
                      "true_density" = purrr::map_dbl(
                        seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 0.1),
                        calculate_gaussian_mixture_prob,
                        mu_vec,
                        sigma_vec,
                        mixture_probs
                      )),
    ggplot2::aes(x = x, y = true_density * 10000)) +
  ggplot2::scale_x_continuous(
    breaks = round(seq(min(mu_vec) - 4 * max(sigma_vec), max(mu_vec) + 4 * max(sigma_vec), 1)))
