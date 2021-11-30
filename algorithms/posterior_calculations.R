library(dplyr)
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

calculate_mixture_posterior <- function(data,
                                        p,
                                        mu,
                                        sigma,
                                        n_components,
                                        alpha,
                                        k_0,
                                        mu_0,
                                        S_0,
                                        v_0) {

  ## PRIORS
  # mixture probs prior: dir(alpha)
  log_p_prior_prob <- sum((alpha - 1) * p)
  #message(paste0("Log prior probability for p is: ", log_p_prior_prob))
  # component means prior: N(mu_0, sigma_0)
  log_mean_prior_prob <- log(dMultivariateNormal(mu, mu_0, k_0 * sigma))
  #message(paste0("Log prior probability for the component means is: ", log_mean_prior_prob))
  # component variances prior: Inv-Wishart(df_0, sigma_0)
  log_sigma_prior_prob <- log(dInvWishart(sigma, v_0, S_0))
  #message(paste0("Log prior probability for the component variance structure is: ", log_sigma_prior_prob))
  
  ## LIKELIHOOD
  # p(zi | p)
  z <- sample(1:n_components, length(data), replace = TRUE, prob = p)
  log_z_likelihood <- data.frame("z" = z) %>%
    dplyr::count(z) %>%
    dplyr::full_join(data.frame("z" = 1:n_components)) %>%
    dplyr::mutate_at("n", tidyr::replace_na, 0) %>%
    dplyr::arrange(z) %>%
    cbind(., "p" = p) %>%
    dplyr::summarize(l = sum(n * p)) %>%
    dplyr::pull(l)
  #message(paste0("Log Z likelihood is :", log_z_likelihood))
  # p(xi | zi = j, mu_k, sigma_k)
  log_x_likelihood <- purrr::pmap_dbl(list(
    "x" = data,
    "mean" = mu[z],
    "sd" = sqrt(diag(sigma)[z])),
    dnorm) %>%
    log() %>%
    sum()
  #message(paste0("Log X likelihood is: ", log_x_likelihood))

  return(sum(log_p_prior_prob, log_mean_prior_prob, log_sigma_prior_prob, log_z_likelihood,
         log_x_likelihood))
}
