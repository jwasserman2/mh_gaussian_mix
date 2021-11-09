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
  return((dnorm(x, mu, sigma) %*% mixture_probs)[1, 1])
}

#' Calculate T mixture density
calculate_heavytail_mixture_prob <- function(x, df, noncentral_param, mixture_probs) {
  vec_lengths <- purrr::map_dbl(list(df, noncentral_param, mixture_probs), length)
  if (any(vec_lengths != vec_lengths[1])) {
    msg <- paste0("Degrees of freedom, noncentral parameters, and mixture probabilities ",
                  "must have corresponding shapes")
    rlang::abort(msg)
  }
  return((dt(x, df = df, ncp = noncentral_param) %*% mixture_probs)[1, 1])
}
