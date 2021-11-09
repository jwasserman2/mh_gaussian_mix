library(testthat)

source("./algorithms/data_sim.R")
mu_vec <- c(-5, 0, 5)
sigma_vec <- rep(2, 3)
df_vec <- rep(10, 3)
mixture_probs <- c(0.45, 0.1, 0.45)

testthat::test_that(
  "gaussian_data_sim_correct_value_one_param",
  testthat::expect_equal((dnorm(0, mu_vec, sigma_vec) %*% mixture_probs)[1, 1],
                         calculate_gaussian_mixture_prob(0, mu_vec, sigma_vec, mixture_probs))
)

testthat::test_that(
  "gaussian_data_sim_unequal_lengths",
  testthat::expect_error(calculate_gaussian_mixture_prob(0, c(mu_vec, 0), sigma_vec, mixture_probs),
                         regexp = "Mu,")
)

testthat::test_that(
  "gaussian_data_sim_correct_dim",
  testthat::expect_equal(1,
                         length(calculate_gaussian_mixture_prob(0, mu_vec, sigma_vec, mixture_probs)))
)

testthat::test_that(
  "t_data_sim_correct_value_one_param",
  testthat::expect_equal((dt(0, df_vec, mu_vec) %*% mixture_probs)[1, 1],
                         calculate_heavytail_mixture_prob(0, df_vec, mu_vec, mixture_probs))
)

testthat::test_that(
  "t_data_sim_unequal_lengths",
  testthat::expect_error(calculate_heavytail_mixture_prob(0, df_vec, mu_vec, c(mixture_probs, 0.1)),
                         regexp = "Degrees")
)