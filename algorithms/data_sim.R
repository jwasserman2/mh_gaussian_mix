library(rlang)

#' Calculate the true density for a mixture distribution
calculate_true_density <- function(x, true_density_kwargs) {
  run_kwargs <- true_density_kwargs[names(true_density_kwargs) != "gaussian_mixture"]
  run_kwargs["x"] <- x
  if (true_density_kwargs[["gaussian_mixture"]]) {
    return(do.call(calculate_gaussian_mixture_prob, run_kwargs))
  } else {
    return(do.call(calulate_heavytail_mixture_prob, run_kwargs))
  }
}

#' Calculate Gaussian mixture density
calculate_gaussian_mixture_prob <- function(x, mu, sigma, mixture_probs) {
  vec_lengths <- purrr::map_dbl(list(mu, sigma, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    rlang::abort("Mu, sigma, and mixture probabilities must have corresponding shapes")
  }
  return(purrr::map_dbl(x, ~ dnorm(.x, mu, sqrt(sigma)) %*% mixture_probs))
}

#' Calculate T mixture density
calculate_heavytail_mixture_prob <- function(x, df, noncentral_param, mixture_probs) {
  vec_lengths <- purrr::map_dbl(list(df, noncentral_param, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    msg <- paste0("Degrees of freedom, noncentral parameters, and mixture probabilities ",
                  "must have corresponding shapes")
    rlang::abort(msg)
  }
  return(purrr::map_dbl(x, ~ dt(.x, df = df, ncp = noncentral_param) %*% mixture_probs))
}

#' Generate a random sample or empirical cdf from a Gaussian mixture 
generate_normal_mixture_data <- function(n, mu, sigma, mixture_probs, random = TRUE, seed = NA) {
  if (is.na(seed) & random) {
    seed <- round(runif(1), 7) * 10e7
    set.seed(seed)
    message(paste0("Seed set to ", seed))
  }

  vec_lengths <- purrr::map_dbl(list(mu, sigma, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    rlang::abort("Mixture parameters and probabilities must have corresponding shapes")
  }

  if (random) {
    dist_idx <- sample(1:length(mixture_probs), n, replace = TRUE, prob = mixture_probs)
    data <- purrr::map2_dbl(mu[dist_idx], sigma[dist_idx], rnorm, n = 1)
  } else {
    set.seed(2045)
    dist_idx <- sample(1:length(mixture_probs), N_FOR_CDF, replace = TRUE, prob = mixture_probs)
    cdf <- purrr::map2_dbl(mu[dist_idx], sigma[dist_idx], rnorm, n = 1)
    data <- sample(cdf, n)
  }
  
  return(data)
}

#' Generate a random sample or empirical cdf from a T mixture 
generate_t_mixture_data <- function(n, ncp, df, mixture_probs, random = TRUE, seed = NA) {
  if (is.na(seed) & random) {
    seed <- round(runif(1), 7) * 10e7
    set.seed(seed)
    message(paste0("Seed set to ", seed))
  }
  
  vec_lengths <- purrr::map_dbl(list(ncp, df, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    rlang::abort("Mixture parameters and probabilities must have corresponding shapes")
  }
  
  if (random) {
    dist_idx <- sample(1:length(mixture_probs), n, replace = TRUE, prob = mixture_probs)
    data <- purrr::map2_dbl(df[dist_idx], ncp[dist_idx], rt, n = 1)
  } else {
    set.seed(2045)
    dist_idx <- sample(1:length(mixture_probs), N_FOR_CDF, replace = TRUE, prob = mixture_probs)
    cdf <- purrr::map2_dbl(df[dist_idx], ncp[dist_idx], rt, n = 1)
    data <- sample(cdf, n)
  }
  
  return(data)
}
