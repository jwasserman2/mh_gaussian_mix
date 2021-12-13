SEED <- 2380
chains <- list()
for (seed in SEED:(SEED + 9)) {
  message(paste0("Running chain ", seed - SEED + 1))
  chains[[seed - SEED + 1]] <- run_metropolized_independent_sampler(data = data,
                                             priors = priors,
                                             n_components = N_COMPONENTS,
                                             proposal_param = 15,
                                             n_iters = 10000,
                                             true_sigma = diag(TRUE_VARIANCES),
                                             seed = seed)
}
json <- jsonlite::toJSON(chains)
file_suffix <- paste(c("mis_sampler", "normal", 10, 10000, 29791043),
                     collapse = "_")
json_file <- paste0("./jsons/", file_suffix)
jsonlite::write_json(json, json_file)


data <- generate_normal_d_mixture_data(
  500,
  matrix(c(-1, -1, 1, 1), nrow = 2, byrow = T),
  list(diag(2), diag(2)),
  mixture_probs = c(.25, 0.75),
  seed = 4898)
t_data_1 <- generate_t_d_mixture_data(
  500,
  matrix(c(-1, -1, 1, 1), nrow = 2, byrow = T),
  c(8,8),
  mixture_probs = c(.25, 0.75),
  random = TRUE,
  seed = 4898)
true_mu <- matrix(c(-10, -5, 0, 0, 5, 10), nrow = 3, byrow = T)
data_2 <- generate_normal_d_mixture_data(
  500,
  true_mu,
  list(2 * diag(2), 2 * diag(2), 2 * diag(2)),
  mixture_probs = c(0.4, 0.2, 0.4),
  seed = 4898)
t_data_2 <- generate_t_d_mixture_data(
  500,
  true_mu,
  rep(3, 3),
  mixture_probs = c(0.4, 0.2, 0.4),
  seed = 4898)
data_3_mu <- matrix(c(sample(c(1, -1), 10, replace = T) * sample(1:10, 10, replace = T)), ncol = 2)
data_3_sigma <- list(4 * diag(2), 4 * diag(2), 4 * diag(2), 4 * diag(2), 4 * diag(2))
data_3 <- generate_normal_d_mixture_data(
  500,
  data_3_mu,
  data_3_sigma,
  mixture_probs = rep(0.2, 5),
  seed = 4898)
t_data_3_df <- rep(6, 5)
t_data_3_sigma <- list(diag(2) * 6/4, diag(2) * 6/4, diag(2) * 6/4, diag(2) * 6/4, diag(2) * 6/4)
t_data_3 <- generate_t_d_mixture_data(
  500,
  data_3_mu,
  t_data_3_df,
  mixture_probs = rep(0.2, 5),
  seed = 4899)
N_COMPONENTS <- 5
priors <- list(
  "alpha" = rep(1 / N_COMPONENTS, N_COMPONENTS),
  "mu_0" = matrix(rep(0, N_COMPONENTS * 2), nrow = 1),
  "tau_2" = 5 * diag(N_COMPONENTS)
)
chains <- list()
SEED <- 19200
initial_params <- initial_anneal(data_3, data_3_sigma, -20, 20, 100)
t_initial_params <- initial_anneal(t_data_3, t_data_3_sigma, -20, 20, 100)
data_1_initial_params <- initial_anneal(data, list(diag(2), diag(2)), -20, 20)
t_data_1_initial_params <- initial_anneal(t_data_1, list(4/3 * diag(2), 4/3 * diag(2)), -20, 20, 200)
for (seed in SEED:(SEED + 9)) {
  chains[[seed - SEED + 1]] <- run_d_rw_metropolis(t_data_1, d = 2,
                                priors, 2, proposal_param = 0.01,
                                n_iters = 5000,
                                true_sigma = list(4/3 * diag(2), 4/3 * diag(2)),
                                initial_p = t_data_1_initial_params[
                                  1:2],
                                initial_mu = t_data_1_initial_params[
                                  3:6],
                                seed = seed)
}

json <- jsonlite::toJSON(chains)
file_suffix <- paste(c("mip", "simple", "normal", 4, 5000),
                     collapse = "_")
json_file <- paste0("./jsons/", file_suffix)
jsonlite::write_json(json, json_file)

for (seed in SEED:(SEED + 8)) {
  message(paste0("Running chain ", seed - SEED + 1))
  chains[[seed - SEED + 1]] <- run_multiple_d_independent_proposals(data,
                                                                    K = 6,
                                                                    d = 2,
                                                                    priors = priors,
                                                                    n_components = N_COMPONENTS,
                                                                    proposal_param = 0.05,
                                                                    n_iters = 5000,
                                                                    initial_p = data_1_initial_params[
                                                                      1:2],
                                                                    initial_mu = data_1_initial_params[
                                                                      3:6],
                                                                    true_sigma = list(diag(2),  diag(2)),
                                                                    seed = seed)
}

data_1_initial_params <- initial_anneal(data, list(diag(2), diag(2)), -20, 20)
data_2_initial_params <- initial_anneal(data_2, list(2 * diag(2), 2 * diag(2), 2 * diag(2)), -20, 20, 100)
t_initial_params <- initial_anneal(t_data_3, t_data_3_sigma, -20, 20, 100)
independent_d_chains <- run_multiple_d_independent_sampler(data_2,
                                                           d = 2,
                                                           priors,
                                                           n_components = 3,
                                                           proposal_param = 1,
                                                           n_iters = 2000,
                                                           data_2_initial_params[1:3],
                                                           data_2_initial_params[4:9],
                                                           true_sigma = list(2 * diag(2), 2 * diag(2), 2 * diag(2)),
                                                           seed = 999)

# combined chains density plots
heat_map_posterior_predictive_density <- function(data, p_chains, mu_chains, true_p, true_mu,
                                                   true_sigma_or_df, density_scale, true_dist) {
  n_steps <- dim(p_chains[[1]])[1]
  burn_in_iters <- ceiling(n_steps / 2)
  density_space <- expand.grid("x" = density_scale, "y" = density_scale)
  
  p_estimates <- purrr::map(
    p_chains,
    ~ .x[(burn_in_iters + 1):n_steps,]) %>%
    Reduce(rbind, .) %>%
    apply(., 2, mean)
  
  mu_estimates <- purrr::map(
    mu_chains,
    ~ .x[(burn_in_iters + 1):n_steps,]) %>%
    Reduce(rbind, .) %>%
    apply(., 2, mean) %>%
    matrix(., ncol = 2, byrow = T)
  
  density_func <- switch(
    true_dist,
    "normal" = dMultivariateNormal,
    "t" = mvtnorm::dmvt
  )
  
  estimated_densities <- purrr::pmap(list(p_estimates, mu_estimates[,1], mu_estimates[,2]),
                                function(X, Y, Z) {
                                  return(
                                    X * purrr::map2_dbl(density_space$x, density_space$y,
                                                        function(x, y) {
                                                          args <- list("x" = c(x, y))
                                                          if (true_dist == "normal") {
                                                            args[["mu"]] <- c(Y, Z)
                                                            args[["Sigma"]] <- true_sigma_or_df
                                                          } else {
                                                            args[["delta"]] <- c(Y, Z)
                                                            args[["df"]] <- true_sigma_or_df
                                                            args[["log"]] <- F
                                                          }
                                                          return(do.call(density_func, args))
                                                        })
                                  )
                                }) %>%
    Reduce(cbind, .)
  density_space$estimated <- estimated_densities[,1] + estimated_densities[,2] # +
    # estimated_densities[,3] + estimated_densities[,4] + estimated_densities[,5]
  density_space$estimated <- density_space$estimated / sum(density_space$estimated)

  true_densities <- purrr::pmap(list(true_p, true_mu[,1], true_mu[,2]),
              function(X, Y, Z) {
                return(
                  X * purrr::map2_dbl(density_space$x, density_space$y,
                                      function(x, y) {
                                        args <- list("x" = c(x, y))
                                        if (true_dist == "normal") {
                                          args[["mu"]] <- c(Y, Z)
                                          args[["Sigma"]] <- true_sigma_or_df
                                        } else {
                                          args[["delta"]] <- c(Y, Z)
                                          args[["df"]] <- true_sigma_or_df
                                          args[["log"]] <- F
                                        }
                                        return(do.call(density_func, args))
                                      })
                )
              }) %>%
    Reduce(cbind, .)
  density_space$true <- true_densities[,1] + true_densities[,2] # + true_densities[,3] +
    # true_densities[,4] + true_densities[,5]
  density_space$true <- density_space$true / sum(density_space$true)

  out <- ggplot(density_space, aes(x = x, y = y)) +
    geom_raster(aes(fill = true)) +
    geom_contour(aes(z = estimated), color = "#d48f04", size = 0.25) +
    theme(panel.background = element_blank()) +
    labs(x = "X", y = "Y",
         subtitle = paste0("Estimated P's: ", paste(round(p_estimates, 3), collapse = ", "),
                           ", estimated ", expression(mu), "'s: ",
                           paste(round(mu_estimates, 3), collapse = ", "))) +
    theme(axis.text = element_text(size = 12),
          title = element_text(size = 16)) +
    guides(fill = F)
  
  return(out)
}

heat_map_posterior_predictive_differences <- function(data, p_chains, mu_chains, true_p, true_mu,
                                                   true_sigma_or_df, density_scale, true_dist) {
  n_steps <- dim(p_chains[[1]])[1]
  burn_in_iters <- ceiling(n_steps / 2)
  density_space <- expand.grid("x" = density_scale, "y" = density_scale)
  
  p_estimates <- purrr::map(
    p_chains,
    ~ .x[(burn_in_iters + 1):n_steps,]) %>%
    Reduce(rbind, .) %>%
    apply(., 2, mean)
  
  mu_estimates <- purrr::map(
    mu_chains,
    ~ .x[(burn_in_iters + 1):n_steps,]) %>%
    Reduce(rbind, .) %>%
    apply(., 2, mean) %>%
    matrix(., ncol = 2, byrow = T)
  
  density_func <- switch(
    true_dist,
    "normal" = dMultivariateNormal,
    "t" = mvtnorm::dmvt
  )
  
  estimated_densities <- purrr::pmap(list(p_estimates, mu_estimates[,1], mu_estimates[,2]),
                                     function(X, Y, Z) {
                                       return(
                                         X * purrr::map2_dbl(density_space$x, density_space$y,
                                                             function(x, y) {
                                                               args <- list("x" = c(x, y))
                                                               if (true_dist == "normal") {
                                                                 args[["mu"]] <- c(Y, Z)
                                                                 args[["Sigma"]] <- true_sigma_or_df
                                                               } else {
                                                                 args[["delta"]] <- c(Y, Z)
                                                                 args[["df"]] <- true_sigma_or_df
                                                                 args[["log"]] <- F
                                                               }
                                                               return(do.call(density_func, args))
                                                             })
                                       )
                                     }) %>%
    Reduce(cbind, .)
  density_space$estimated <- estimated_densities[,1] + estimated_densities[,2] # +
    # estimated_densities[,3] + estimated_densities[,4] + estimated_densities[,5]
  density_space$estimated <- density_space$estimated / sum(density_space$estimated)

  true_densities <- purrr::pmap(list(true_p, true_mu[,1], true_mu[,2]),
                                function(X, Y, Z) {
                                  return(
                                    X * purrr::map2_dbl(density_space$x, density_space$y,
                                                        function(x, y) {
                                                          args <- list("x" = c(x, y))
                                                          if (true_dist == "normal") {
                                                            args[["mu"]] <- c(Y, Z)
                                                            args[["Sigma"]] <- true_sigma_or_df
                                                          } else {
                                                            args[["delta"]] <- c(Y, Z)
                                                            args[["df"]] <- true_sigma_or_df
                                                            args[["log"]] <- F
                                                          }
                                                          return(do.call(density_func, args))
                                                        })
                                  )
                                }) %>%
    Reduce(cbind, .)
  density_space$true <- true_densities[,1] + true_densities[,2] # + true_densities[,3] +
    # true_densities[,4] + true_densities[,5]
  density_space$true <- density_space$true / sum(density_space$true)
  density_space$diff <- density_space$estimated - density_space$true
  
  wasserstein_1 <- sum(abs(density_space$diff))

  out <- ggplot(density_space, aes(x = x, y = y)) +
    geom_raster(aes(fill = diff)) +
    geom_contour(aes(z = estimated), bins = 5, color = "#d48f04", size = 0.25) +
    scale_fill_gradient2() +
    theme(panel.background = element_blank()) +
    labs(x = "X", y = "Y",
         subtitle = paste0("Approximate 1-Wasserstein distance: ", round(wasserstein_1, 3))) +
    theme(axis.text = element_text(size = 12),
          title = element_text(size = 16)) +
    guides(fill = guide_colorbar(title = "Difference", label.theme = element_text(size = 10)))
  
  return(out)
}

heat_map_posterior_mean_estimates <- function(mean_chains, true_parameter, density_scale) {
  n_steps <- dim(mean_chains[[1]])[1]
  burn_in_iters <- ceiling(n_steps / 2)
  density_space <- expand.grid("x" = density_scale, "y" = density_scale)
  
  estimates <- purrr::map(mean_chains, ~ .x[(burn_in_iters + 1):n_steps,]) %>%
    Reduce(rbind, .) %>%
    as.data.frame()
  
  mean_elementwise_bias <- purrr::map(
    mean_chains,
    function(x) {
      (matrix(apply(x[(burn_in_iters + 1):n_steps,], 2, mean), ncol = 2, byrow = T) -
           true_parameter)^2
    }) %>%
    Reduce(rbind, .) %>%
    mean()

  n <- burn_in_iters
  m <- length(mean_chains)
  W_num <- purrr::map(
    mean_chains,
    ~ apply((.x[(burn_in_iters + 1):n_steps,] -
               apply(.x[(burn_in_iters + 1):n_steps,], 2, mean))^2, 2, sum)) %>%
    Reduce(rbind, .) %>%
    apply(., 2, sum)
  W <- W_num / m / (n - 1)
  B_div_n_sq_terms <- purrr::map(
    mean_chains,
    ~ apply(.x[(burn_in_iters + 1):n_steps,], 2, mean)
  ) %>%
    Reduce(rbind, .)
  B_div_n <- colSums(t((t(B_div_n_sq_terms) - sapply(estimates, mean))^2)) / (m - 1)
  sigma_hat_2 <- (n - 1) * W / n + B_div_n
  r_hat <- (m + 1) * sigma_hat_2 / m / W - (n - 1) / (m * n)
  
  true <- as.data.frame(true_parameter)
  
  message(paste0("The mean elementwise bias for the mean components is" , mean_elementwise_bias))
  message(paste0("R-hat for each element of the mean components is: ", paste(r_hat, collapse = ", ")))
  
  out <- ggplot(density_space, aes(x = x, y = y)) +
    geom_density_2d(data = estimates, aes(x = V1, y = V2), color = "#d48f04", size = 0.25) +
    geom_point(data = true, aes(x = V1, y = V2), color = "#eb4034", size = 2) +
    theme(panel.background = element_blank()) +
    labs(x = "X", y = "Y") +
    theme(axis.text = element_text(size = 12),
          title = element_text(size = 16))
  
  if (ncol(estimates) > 2) {
    for (comp in 2:(ncol(estimates) / 2)) {
      cols <- paste0("V", c(comp * 2 - 1, comp * 2))
      out <- out +
        geom_density_2d(data = estimates, aes(x = !!sym(cols[1]), y = !!sym(cols[2])),
                        color = "#d48f04", size = 0.25)      
    }
  }
  
  return(out)
}

new_p_calc <- function() {
  purrr::map(
    p_chains,
    function(x) {
      (apply(x[(burn_in_iters + 1):n_steps,], 2, mean) - c(0.75, 0.25))^2
    }) %>% Reduce(rbind, .) %>% mean
}

N_STEPS <- 10000
BURN_IN_ITERS <- N_STEPS / 2 # Gelman and Rubin, 1998

# data read
json <- jsonlite::read_json("jsons/rw_simple_normal_4_5000")
prac_chains <- jsonlite::parse_json(json[[1]])
p_chains <- purrr::map(
  prac_chains,
  function(x) {
    purrr::map(x$p_chain, ~ matrix(unlist(.x), nrow = 1)) %>%
      Reduce(rbind, .)
  }
)
mu_chains <- purrr::map(
  prac_chains,
  function(x) {
    purrr::map(x$mu_chain, ~ matrix(unlist(.x), nrow = 1)) %>%
      Reduce(rbind, .)
  }
)
p_chains <- purrr::map(
  chains,
  ~ .x$p_chain)
mu_estimates <- rbind(mu_chains[[1]][2500:5000,c(9:10, 5:6, 7:8, 1:2, 3:4)],
                     mu_chains[[2]][2500:5000,c(1:2, 9:10, 7:8, 3:4, 5:6)],
                     mu_chains[[3]][2500:5000,c(1:2, 9:10, 7:8, 5:6, 3:4)],
                     mu_chains[[4]][2500:5000,c(7:8, 5:6, 9:10, 3:4, 1:2)]) %>%
  apply(., 2, mean) %>%
  matrix(., ncol = 2, byrow = T)
p_estimates <- rbind(p_chains[[1]][2500:5000,c(5, 3, 4, 1, 2)],
                      p_chains[[2]][2500:5000,c(1, 5, 4, 2, 3)],
                      p_chains[[3]][2500:5000,c(1, 5, 4, 3, 2)],
                      p_chains[[4]][2500:5000,c(4, 3, 5, 2, 1)]) %>%
  apply(., 2, mean)

p_chains <- list(p_chains[[1]],
                 p_chains[[2]][, c(3, 2, 5, 4, 1)],
                 p_chains[[3]],
                 p_chains[[4]])
mu_chains <- list(mu_chains[[1]],
                  mu_chains[[2]][,c(5:6, 3:4, 9:10, 7:8, 1:2)],
                  mu_chains[[3]],
                  mu_chains[[4]])
rw_simple_normal_y_heat_map <- heat_map_posterior_predictive_density(
  data,
  p_chains = p_chains,
  mu_chains = mu_chains,
  true_p = c(0.25, 0.75),
  true_mu = matrix(c(-1, -1, 1, 1), nrow = 2, byrow = T),
  true_sigma_or_df = diag(2),
  density_scale = seq(-14, 14, 0.1),
  true_dist = "normal")
rw_simple_normal_diff_heat_map <- heat_map_posterior_predictive_differences(
  t_data_1,
  p_chains = p_chains,
  mu_chains = mu_chains,
  true_p = c(0.25, 0.75),
  true_mu = matrix(c(-1, -1, 1, 1), nrow = 2, byrow = T),
  true_sigma_or_df = diag(2),
  density_scale = seq(-14, 14, 0.1),
  true_dist = "normal")
rw_simple_normal_mu_heat_map <- heat_map_posterior_mean_estimates(
  mu_chains, data_3_mu, seq(-14, 14, 0.1)
)


initial_em <- function(data, n_components, true_sigma) {
  p <- gtools::rdirichlet(1, rep(4, n_components))
  mu <- cbind(t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))),
                      t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))),
                      t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))))
  log_lik_diff <- 1e7
  
  while(log_lik_diff > 1e-4) {
    log_lik <- sum(log(p[1] * apply(data, 1, dMultivariateNormal, mu[1:2], true_sigma[[1]]))) +
      sum(log(p[2] * apply(data, 1, dMultivariateNormal, mu[3:4], true_sigma[[2]]))) +
      sum(log(p[3] * apply(data, 1, dMultivariateNormal, mu[5:6], true_sigma[[3]])))
    p_x_all <- data.frame(data)
    for (i in 1:n_components) {
      p_x_all[paste0("p_x_z_", i)] <- p[i] * apply(data, 1, dMultivariateNormal,
                                                   mu = mu[(i * d - 1):(i * d)],
                                                   true_sigma[[i]])
    }
    
    p_x_all <- p_x_all %>%
      dplyr::select(tidyselect::starts_with("p_x")) %>%
      t() %>%
      rbind(., "p_x_all" = colSums(.)) %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate_at(vars(tidyselect::starts_with("p_x_z")), list("multinom" = ~ .x / p_x_all)) %>%
      cbind(data, .)
    
    mu <- c()
    p <- c()
    for (i in 1:n_components) {
      mu <- c(mu, sum(p_x_all[,paste0("p_x_z_", i, "_multinom")] * p_x_all[,1]) /
                sum(p_x_all[,paste0("p_x_z_", i, "_multinom")]),
              sum(p_x_all[,paste0("p_x_z_", i, "_multinom")] * p_x_all[,2]) /
                sum(p_x_all[,paste0("p_x_z_", i, "_multinom")]))
      p <- c(p, sum(p_x_all[,paste0("p_x_z_", i, "_multinom")]) / dim(data)[1])
    }
    
    log_lik_new <- sum(log(p[1] * apply(data, 1, dMultivariateNormal, mu[1:2], true_sigma[[1]]))) +
      sum(log(p[2] * apply(data, 1, dMultivariateNormal, mu[3:4], true_sigma[[2]]))) +
      sum(log(p[3] * apply(data, 1, dMultivariateNormal, mu[5:6], true_sigma[[3]])))
    log_lik_diff <- log_lik_new - log_lik
  }
  return(mu)
}

# initial_anneal <- function(data, n_epochs) {
#   lambda <- 0.01
#   gamma <- 0.15
#   tau <- gamma * (lambda / gamma)^(-1/(n_epochs - 1))
#   
#   p <- gtools::rdirichlet(1, rep(4, n_components))
#   mu <- cbind(t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))),
#               t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))),
#               t(MASS::mvrnorm(1, mu = rep(0, 2), Sigma = 25 * diag(2))))
#   params <- c(p, mu)
#   
#   fn <- function(p) {
#     return(
#       -sum(log(p[1] * apply(data, 1, dMultivariateNormal, p[4:5], true_sigma[[1]]) +
#               p[2] * apply(data, 1, dMultivariateNormal, p[6:7], true_sigma[[2]]) +
#               p[3] * apply(data, 1, dMultivariateNormal, p[8:9], true_sigma[[3]])))
#         + 1000 * abs(p[1] + p[2] + p[3] - 1)
#     )
#   }
#   
#   for (epoch in 1:n_epochs) {
#     for (par in 1:4) {
#       tau <- gamma * (lambda / gamma)^((epoch-1) / (n_epochs - 1))
#       if (par == 1) {
#         p <- gtools::rdirichlet(1, rep(4, n_components))
#       } else if (par == 2) {
#         mu[1:2] <- MASS::mvrnorm(1, mu = mu[1:2], Sigma = 25 * diag(2))
#       } else if (par == 3) {
#         mu[3:4] <- MASS::mvrnorm(1, mu = mu[3:4], Sigma = 25 * diag(2))
#       } else {
#         mu[5:6] <- MASS::mvrnorm(1, mu = mu[5:6], Sigma = 25 * diag(2))
#       }
#       new_params <- c(p, mu)
#       e_new_par <- exp(-(fn(new_params) - e_par) / tau)
#       s <- runif(1)
#       if (e_new_par > s) {
#         best_e_par <- e_par <- e_new_par
#         best_params <- params <- new_params
#       }
#     }
#     gamma <- gamma * 0.85
#     if (gamma <= lambda) {
#       break
#     }
#   }
#   return(best_params)
# }
initial_anneal <- function(data, true_sigma, min_mu, max_mu, maxit = 20) {
  fn <- function(p) {
    return(
      # -sum(log(p[1] * apply(data, 1, dMultivariateNormal, p[6:7], true_sigma[[1]]) +
      #          p[2] * apply(data, 1, dMultivariateNormal, p[7:8], true_sigma[[2]]) +
      #          p[3] * apply(data, 1, dMultivariateNormal, p[9:10], true_sigma[[3]]) +
      #          p[4] * apply(data, 1, dMultivariateNormal, p[11:12], true_sigma[[4]]) +
      #          p[5] * apply(data, 1, dMultivariateNormal, p[13:14], true_sigma[[5]])))
      # + 1000 * abs(p[1] + p[2] + p[3] + p[4] + p[5] - 1)
      # -sum(log(p[1] * apply(data, 1, dMultivariateNormal, p[4:5], true_sigma[[1]]) +
      #            p[2] * apply(data, 1, dMultivariateNormal, p[6:7], true_sigma[[2]]) +
      #            p[3] * apply(data, 1, dMultivariateNormal, p[8:9], true_sigma[[3]])))
      -sum(log(p[1] * apply(data, 1, dMultivariateNormal, p[3:4], true_sigma[[1]]) +
                 p[2] * apply(data, 1, dMultivariateNormal, p[5:6], true_sigma[[2]])))
      + 1000 * abs(p[1] + p[2] - 1)
    )
  }
  opt <- GenSA::GenSA(fn = fn,
                      lower = c(rep(0, 2), rep(min_mu, 4)),
                      upper = c(rep(1, 2), rep(max_mu, 4)),
                      control = list("maxit" = maxit))
  message(paste0("Simulated annealing yielded parameters with a penalized negative log likelihood",
                 " of ", opt$value))
  return(opt$par)
}

run_d_rw_metropolis <- function(data,
                                d,
                          priors,
                          n_components,
                          # proposal_alpha,
                          proposal_param,
                          # max_proposal_sigma,
                          # proposal_corr_0,
                          n_iters,
                          initial_p,
                          initial_mu,
                          true_sigma,
                          seed = NA) {
  if (n_iters < 2) {
    rlang::abort("Must run at least 2 iterations in the chain")
  }
  
  # incorrect_priors <- setdiff(names(priors),
  #                             c("alpha", "k_0", "mu_0", "S_0", "v_0", "tau_2"))
  # if (length(incorrect_priors) > 0) {
  #   msg <- paste0("Priors must have names \"alpha\", \"k_0\", \"mu_0\", \"S_0\", ",
  #                 "\"tau_2\", and \"v_0\"")
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
  
  message(paste0("Running Random Walk Metropolis chain with ", n_iters, " steps and seed: ", seed))
  set.seed(seed)
  
  # initial param values
  p_chain <- matrix(initial_p, nrow = 1)
  mu_chain <- matrix(ncol = n_components * d)
  mu_chain[1,] <- initial_mu
  z <- sample(1:n_components, dim(data)[1], replace = T, prob = p_chain[1,])
  # init_sigma_variances <- runif(n_components, min = 0.001, max = max_proposal_sigma)
  # init_sigma_cov <- runif(1, -min(abs(init_sigma_variances)), min(abs(init_sigma_variances))) / proposal_corr_0
  # sigma_chain <- list(diag(init_sigma_variances - init_sigma_cov) + init_sigma_cov)
  
  # final pre-chain prep
  # current_args <- priors
  # current_args <- append(current_args, list(
  #   "n_components" = n_components,
  #   "data" = data,
  #   "p" = as.vector(p_chain),
  #   # "p" = true_p,
  #   "mu" = mu_chain[1,],
  #   # "sigma" = sigma_chain[[1]]))
  #   "sigma" = true_sigma))
  
  acceptance_probs <- c()
  acceptances <- c()
  for (i in 2:n_iters) {
    if (i %% round(n_iters / 5) == 0) {
      message(paste0("Proposing step ", i))
    }
    # proposals
    # symmetric Dirichlet with concentration parameter alpha for mixture probs p
    # p_proposal <- gtools::rdirichlet(1, rep(proposal_alpha, n_components))
    # symmetric distributions for drawing variances and covariance independently
    # proposal_variances <- runif(n_components, min = 0.001, max = max_proposal_sigma)
    # proposal_cov <- runif(1, -min(abs(proposal_variances)), min(abs(proposal_variances))) / proposal_corr_0
    # sigma_proposal <- diag(proposal_variances - proposal_cov) + proposal_cov
    p_chain <- rbind(p_chain, calculate_gibbs_p_updates(z, p_chain[i - 1,], n_components, priors$alpha))
    z <- calculate_gibbs_z_d_updates(data, p_chain[i,], mu_chain[i - 1,], true_sigma, n_components, d)
    
    # proposal_args <- current_args
    # proposal_args[c("p", "mu", "sigma")] <- list(
    #   "p" = as.vector(p_proposal),
    #   # "p" = true_p,
    #   "mu" = mu_proposal,
    #   # "sigma" = sigma_proposal
    #   "sigma" = true_sigma
    # )
    # 
    # log hastings ratio log((h(x)/c) / (h(y)/c)) = log(h(x)) - log(h(y))
    # h_current <- suppressMessages(do.call(calculate_mixture_posterior, current_args))
    # h_current <- do.call(calculate_mixture_posterior, current_args)
    # h_proposal <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
    # h_proposal <- do.call(calculate_mixture_posterior, proposal_args)
    h_current <- max(calculate_mh_gibbs_d_mean_posterior(
      data, 2, p_chain[i,], z, mu_chain[i-1,], true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    # symmetric Gaussian with mean 0 and variance proposal_param^2 * I
    mu_proposal <- matrix(c(
      MASS::mvrnorm(1, mu = mu_chain[i-1,1:2], Sigma = proposal_param * diag(d)),
      MASS::mvrnorm(1, mu = mu_chain[i-1,3:4], Sigma = proposal_param * diag(d))),
      # MASS::mvrnorm(1, mu = mu_chain[i-1,5:6], Sigma = proposal_param * diag(d)),
      # MASS::mvrnorm(1, mu = mu_chain[i-1,7:8], Sigma = proposal_param * diag(d)),
      # MASS::mvrnorm(1, mu = mu_chain[i-1,9:10], Sigma = proposal_param * diag(d))),
                          nrow = 1)
    h_proposal <- max(calculate_mh_gibbs_d_mean_posterior(
      data, 2, p_chain[i,], z, mu_proposal, true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    log_hastings_ratio <- h_proposal - h_current

    acceptance_prob <- min(1, exp(log_hastings_ratio))
    accept <- as.logical(rbinom(1, 1, acceptance_prob))
    if (accept) {
      # p_chain <- rbind(p_chain, p_proposal, deparse.level = 0)
      mu_chain <- rbind(mu_chain, mu_proposal, deparse.level = 0)
      # sigma_chain[[i]] <- sigma_proposal
      # current_args <- proposal_args
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

run_multiple_d_independent_sampler <- function(data,
                                               d,
                                               priors,
                                               n_components,
                                               proposal_param,
                                               n_iters,
                                               initial_p,
                                               initial_mu,
                                               true_sigma,
                                               seed = NA) {
  if (n_iters < 2) {
    rlang::abort("Must run at least 2 iterations in the chain")
  }
  
  incorrect_priors <- setdiff(names(priors),
                              c("alpha", "mu_0", "tau_2"))
  if (length(incorrect_priors) > 0) {
    rlang::abort("Priors must have names \"alpha\", \"mu_0\", and \"tau_2\"")
  }
  
  if (is.na(seed)) {
    seed <- round(runif(1) * 1e7)
  }
  
  message(paste0("Running MIS chain with ", n_iters, " steps and seed: ", seed))
  set.seed(seed)
  
  # initial param values
  p_chain <- matrix(initial_p, nrow = 1)
  mu_chain <- matrix(ncol = n_components * d)
  mu_chain[1,] <- initial_mu
  z <- sample(1:n_components, dim(data)[1], replace = T, prob = p_chain[1,])
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
    p_chain <- rbind(p_chain, calculate_gibbs_p_updates(z, p_chain[i - 1,], n_components, priors$alpha))
    z <- calculate_gibbs_z_d_updates(data, p_chain[i,], mu_chain[i - 1,], true_sigma, n_components, d)
    
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
    # mu_proposal <- matrix(mvtnorm::rmvt(n_components, df = proposal_param), nrow = 1)
    # symmetric Gaussian with mean 0 and variance proposal_param^2 * I
    mu_proposal <- matrix(
      MASS::mvrnorm(n_components, mu = rep(0, d), Sigma = proposal_param^2 * diag(d)),
      nrow = 1)
    
    h_y <- max(calculate_mh_gibbs_d_mean_posterior(
      data, 2, p_chain[i,], z, mu_proposal, true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    # q_y <- as.vector(calculate_log_mis_proposal_prob(mu_proposal, proposal_param))
    q_y <- as.vector(calculate_log_mis_d_proposal_prob(
      mu_proposal, n_components, d, rep(0, d), proposal_param^2 * diag(d)))

    h_x <- max(calculate_mh_gibbs_d_mean_posterior(
      data, 2, p_chain[i,], z, mu_chain[i-1,], true_sigma, n_components, priors$mu_0, priors$tau_2),
      -1e7)
    # q_x <- as.vector(calculate_log_mis_d_proposal_prob(mu_chain[i-1,], proposal_param))
    q_x <- as.vector(calculate_log_mis_d_proposal_prob(
      mu_chain[i-1,],n_components, d, rep(0, d), proposal_param^2 * diag(d)))
    
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
    acceptance_probs <- c(acceptance_probs, acceptance_prob)
  }
  return(list(
    "p_chain" = p_chain,
    "mu_chain" = mu_chain,
    "acceptance_probs" = acceptance_probs,
    # "sigma_chain" = sigma_chain,
    "acceptances" = acceptances))
}

#' Calculate g(x) or g(y) for the Metropolized Independent Sampler
calculate_log_mis_d_proposal_prob <- function(#p, mu, sigma,
  mu, n_components, d, mvn_mu_parameter, mvr_sigma_parameter) { #dir_parameter, 
  #wish_df_parameter, wish_scale_mat_parameter) {
  # mu, t_df_parameter) {
  log_prob_sum <- 0
  for (k in 1:n_components) {
    # log_prob_sum <- log_prob_sum + mvtnorm::dmvt(mu[(k * d - 1):(k * d)], df = t_df_parameter)
    log_prob_sum <- log_prob_sum +
      log(dMultivariateNormal(mu[(k * d - 1):(k * d)], mvn_mu_parameter[(k * d - 1):(k * d)],
                              mvr_sigma_parameter[[k]]))
  }
  return(log_prob_sum)
}

run_multiple_d_independent_proposals <- function(data, d, K, priors, n_components,
                                               # proposal_alpha,
                                               proposal_param,
                                               # max_proposal_sigma,
                                               # proposal_corr_0,
                                               n_iters,
                                               initial_p,
                                               initial_mu,
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
  p_chain <- matrix(initial_p, nrow = 1)
  mu_chain <- matrix(ncol = n_components * d)
  mu_chain[1,] <- initial_mu
  z <- sample(1:n_components, dim(data)[1], replace = T, prob = p_chain[1,])
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
    z <- calculate_gibbs_z_d_updates(data, p_chain[i,], mu_chain[i - 1,], true_sigma, n_components, d)
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
    
    mu_proposals <- cbind(
      MASS::mvrnorm(K, mu = mu_chain[i-1,1:2], Sigma = proposal_param * diag(d)),
      MASS::mvrnorm(K, mu = mu_chain[i-1,3:4], Sigma = proposal_param * diag(d))
      # MASS::mvrnorm(K, mu = mu_chain[i-1,5:6], Sigma = proposal_param * diag(d)),
      # MASS::mvrnorm(K, mu = mu_chain[i-1,7:8], Sigma = proposal_param * diag(d)),
      # MASS::mvrnorm(K, mu = mu_chain[i-1,9:10], Sigma = proposal_param * diag(d))
    )
    
    log_weights <- c()
    for (k in 1:K) {
      # proposal_args <- current_args
      # proposal_args[c("p", "mu", "sigma")] <- list(
      #   "p" = p_proposals[k,],
      #   "mu" = mu_proposals[k,],
      #   "sigma" = sigma_proposals[[k]]
      # )
      # h_k <- suppressMessages(do.call(calculate_mixture_posterior, proposal_args))
      h_k <- max(calculate_mh_gibbs_d_mean_posterior(
        data, 2, p_chain[i,], z, mu_proposals[k,], true_sigma, n_components, priors$mu_0, priors$tau_2),
        -1e7)

      q_x_y <- calculate_log_mis_d_proposal_prob(mu_proposals[k,], n_components, d,
                                                 mu_chain[i - 1,],
                                                 list(proposal_param * diag(d),
                                                      proposal_param * diag(d)))#,
                                                      # proposal_param * diag(d),
                                                      # proposal_param * diag(d),
                                                      # proposal_param * diag(d)))
      l_x_y <- -q_x_y
      log_weights <- c(log_weights, h_k)
    }
    
    # y_idx <- as.vector(rmultinom(1, 1, (-1 / log_weights) / sum(-1 / log_weights)))
    if (all(exp(log_weights) == 0)) {
      m_probs <- rep(1 / K, K)
    } else if (sum(exp(log_weights)) != 0) {
      m_probs <- exp(log_weights) / sum(exp(log_weights))
    } else {
      m_probs <- (-1 / log_weights) / sum(-1 / log_weights)
    }
    y_idx <- as.vector(rmultinom(1, 1, m_probs))
    y <- list(# "p" = p_proposals[which.max(y_idx), ],
      "mu" = mu_proposals[which.max(y_idx), ])
    # "sigma" = sigma_proposals[[which.max(y_idx)]])
    
    # p_star_proposals <- rbind(gtools::rdirichlet(K - 1, rep(proposal_alpha, n_components)),
    #                           p_chain[i - 1,])
    # mu_star_proposals <- rbind(matrix(rnorm((K - 1) * n_components, mean = y$mu,
    #                                         sd = proposal_param), nrow = K - 1),
    #                            mu_chain[i - 1,])
    mu_star_proposals <- rbind(
      cbind(MASS::mvrnorm(K - 1, mu = y$mu[1:2], Sigma = proposal_param * diag(d)),
            MASS::mvrnorm(K - 1, mu = y$mu[3:4], Sigma = proposal_param * diag(d))),
            # MASS::mvrnorm(K - 1, mu = y$mu[5:6], Sigma = proposal_param * diag(d)),
            # MASS::mvrnorm(K - 1, mu = y$mu[7:8], Sigma = proposal_param * diag(d)),
            # MASS::mvrnorm(K - 1, mu = y$mu[9:10], Sigma = proposal_param * diag(d))),
      mu_chain[i-1,])
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
      h_k <- max(calculate_mh_gibbs_d_mean_posterior(
        data, 2, p_chain[i,], z, mu_star_proposals[k,], true_sigma, n_components, priors$mu_0, priors$tau_2),
        -1e7)
      # q_x_y <- calculate_log_mip_proposal_prob(p_star_proposals[k,], mu_star_proposals[k,],
      #                                          sigma_star_proposals[[k]],
      #                                          rep(proposal_alpha, n_components),
      #                                          y$mu, proposal_param,
      #                                          max_proposal_sigma)
      q_x_y <- calculate_log_mis_d_proposal_prob(mu_star_proposals[k,], n_components, d,
                                                 y$mu, list(proposal_param * diag(d),
                                                            proposal_param * diag(d)))#,
                                                            # proposal_param * diag(d),
                                                            # proposal_param * diag(d),
                                                            # proposal_param * diag(d)))
      l_x_y <- -q_x_y
      log_star_weights <- c(log_star_weights, h_k)
    }
    
    # log hastings ratio log(sum(weights(y))) - log(sum(weights(x)))
    if (log_weights)
    log_hastings_ratio <- log(sum(exp(log_weights / 10000))) - log(sum(exp(log_star_weights / 10000)))

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

calculate_mh_gibbs_d_mean_posterior <- function(data, d, p, z, mu, sigma, n_components, mu_0, tau_2) {
  log_posterior_prob_sum <- 0
  for (k in 1:n_components) {
    # k_mean_numerator <- (mu_0[(k * d - 1):(k * d)] / diag(tau_2) +
    #   colSums(data[z == k,]) / diag(sigma[[k]]))
    # k_mean_denominator <- 1 / diag(tau_2) + sum(z == k) / diag(sigma[[k]])
    # k_mean <- k_mean_numerator / k_mean_denominator
    # k_variance <- 1 / k_mean_denominator
    if (sum(z == k) > 1) {
      x_z <- data[z == k,]
      k_mean <- colSums(x_z) / dim(x_z)[1]
      k_variance <- 1 / max(1e-4, sum(z == k)) * sigma[[k]]
    } else if (sum(z == k) == 1) {
      k_mean <- data[z == k,]
      k_variance <- sigma[[k]]
    } else {
      k_mean <- rep(0, d)
      k_variance <- sigma[[k]]
    }

    log_posterior_prob_sum <- log_posterior_prob_sum +
      log(dMultivariateNormal(mu[(k * d - 1):(k * d)], k_mean, k_variance))
  }
  
  return(log_posterior_prob_sum)
}

#' Gibbs sampler for latent assignments Z
#' @return numeric vector of Z draws
calculate_gibbs_z_d_updates <- function(data, p, mu, sigma, n_components, d) {
  z <- sample(1:n_components, dim(data)[1], replace = TRUE, prob = p)
  p_x_all <- data.frame(cbind(data, "z" = z, "p_z" = p[z]))
  for (i in 1:n_components) {
    p_x_all[paste0("p_x_z_", i)] <- p[i] * apply(data, 1, dMultivariateNormal,
                                                 mu = mu[(i * d - 1):(i * d)],
                                                 sigma[[i]])
  }
  
  z_draws <- p_x_all %>%
    dplyr::select(tidyselect::starts_with("p_x")) %>%
    t() %>%
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
