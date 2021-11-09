library(testthat)
library(mockery)

source("./algorithms/rw_metropolis.R")
mu_vec <- c(-5, 0, 5)
sigma_vec <- rep(2, 3)
mixture_probs <- c(0.45, 0.1, 0.45)
proposal_sigma <- 5
n_iters <- 10

testthat::test_that("rw_metropolis_randomly_initialized", {
  mock_runif <- mockery::mock(1)
  testthat::with_mock(runif = mock_runif, {
    rw_metropolis(1, mu_vec, sigma_vec, mixture_probs, proposal_sigma, n_iters, seed = 2045)
  })
  mockery::expect_called(mock_runif, 1)
  mockery::expect_args(mock_runif, 1, n = 1)
})

testthat::test_that("rw_metropolis_correct_always_accept_acceptance_prob", {
  mock_runif <- mockery::mock(0.5)
  mock_rnorm <- mockery::mock(-.5)
  mock_rbinom <- mockery::mock(1)
  testthat::with_mock(runif = mock_runif, rnorm = mock_rnorm, rbinom = mock_rbinom, {
    chain <- rw_metropolis(1, mu_vec, sigma_vec, mixture_probs, proposal_sigma, n_iters = 2, seed = 2045)
  })
  mockery::expect_args(mock_rbinom, 1, n = 1, size = 1, prob = 1)
  testthat::expect_equal(chain, rbind(0.5, 0))
})

testthat::test_that("rw_metropolis_correct_maybe_accept_acceptance_prob", {
  mock_runif <- mockery::mock(0.5)
  mock_rnorm <- mockery::mock(10)
  mock_rbinom <- mockery::mock(0)
  testthat::with_mock(runif = mock_runif, rnorm = mock_rnorm, rbinom = mock_rbinom, {
    chain <- rw_metropolis(1, mu_vec, sigma_vec, mixture_probs, proposal_sigma, n_iters = 2, seed = 2045)
  })
  mockery::expect_args(mock_rbinom, 1, n = 1, size = 1, prob = 0.002409815)
  testthat::expect_equal(chain, rbind(0.5, 0.5))
})
