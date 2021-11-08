#' Calculate the true density for a mixture distribution
calculate_true_density <- function(x, gaussian_mixture, mu, sigma, mixture_probs) {
  if (gaussian_mixture) {
    return(calculate_gaussian_mixture_prob(x, mu, sigma, mixture_probs))
  } else {
    return(calulate_heavytail_mixture_prob(x, mu, sigma, mixture_probs))
  }
}

#' Calculate Gaussian mixture density
calculate_gaussian_mixture_prob <- function(x, mu, sigma, mixture_probs) {
  vec_lengths <- purrr::map_dbl(list(mu, sigma, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    rlang::abort("Mu, sigma, and mixture probabilities must have corresponding shapes")
  }
  return((dnorm(x, mu, sigma) %*% mixture_probs)[1, 1])
}
