# Generate data
N_COMPONENTS <- 2
N_CHAIN <- 10000
N_DATA <- 500
true_mu <- c(-3, 3)
true_sds <- rep(2, N_COMPONENTS)
true_sigma <- diag(true_sds^2)
true_p <- c(.25, .75)
data <- generate_normal_mixture_data(N_DATA, true_mu, true_sds, true_p)
t_data <- generate_t_mixture_data(N_DATA, true_mu, true_sds, true_p)

# priors
priors <- list(
  "alpha" = rep(1 / N_COMPONENTS, N_COMPONENTS),
  # "k_0" = 10,
  "mu_0" = c(-5, 5),
  # "S_0" = diag(N_COMPONENTS) * 5,
  # "v_0" = N_COMPONENTS + 1,
  "tau_2" = 5)
proposal_alpha <- 5
proposal_sd <- 1
proposal_df <- 100
max_proposal_sigma <- 10
proposal_corr_0 <- 100

posterior <- rw_metropolis(data = data,
                            priors = priors,
                            n_components = N_COMPONENTS,
                            proposal_sd = 0.5,
                            n_iters = N_CHAIN,
                            true_sigma = true_sigma,
                            seed = 2045)
posterior <- rw_metropolis(data, priors, N_COMPONENTS, proposal_alpha, proposal_sd, 
                           max_proposal_sigma, proposal_corr_0, N_CHAIN,
                           true_sigma, true_p, seed = 2045)
sum(posterior$acceptances) / length(posterior$acceptances)
setNames(data.frame(posterior$p_chain), paste0("p_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("p_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(posterior$p_chain), paste0("p_", 1:N_COMPONENTS)) %>%
  tidyr::pivot_longer(cols = paste0("p_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))
setNames(data.frame(posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))
purrr::map_dfr(posterior$sigma_chain, ~ data.frame(t(sqrt(diag(.x))))) %>%
  setNames(paste0("sigma_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("sigma_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
mis_posterior <- run_metropolized_independent_sampler(data, priors, N_COMPONENTS,
                                                      proposal_df = 10,
                                                      N_CHAIN, true_sigma, seed = 2045)
sum(mis_posterior$acceptances) / length(mis_posterior$acceptances)
setNames(data.frame(mis_posterior$p_chain), paste0("p_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("p_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(mis_posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(mis_posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_density(aes(x = value, color = name))
purrr::map_dfr(mis_posterior$sigma_chain, ~ data.frame(t(sqrt(diag(.x))))) %>%
  setNames(paste0("sigma_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("sigma_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))

# Multiple Independent Proposals
mip_posterior <- run_multiple_independent_proposals(data, 5, priors, N_COMPONENTS,
                                                    proposal_alpha, proposal_sd,
                                                    max_proposal_sigma, proposal_corr_0,
                                                    N_CHAIN, seed = 2045)
sum(mip_posterior$acceptances) / length(mip_posterior$acceptances)
setNames(data.frame(mip_posterior$p_chain), paste0("p_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("p_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(mip_posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
purrr::map_dfr(mip_posterior$sigma_chain, ~ data.frame(t(sqrt(diag(.x))))) %>%
  setNames(paste0("sigma_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("sigma_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
true_cdf <- generate_normal_mixture_data(N_DATA, true_mu, true_sigma, true_p, random = F)

json <- jsonlite::read_json("./jsons/rw_metropolis_normal_20_20000_hash_27453")
chains <- jsonlite::parse_json(json[[1]])
p_chains <- purrr::map(
  chains,
  function(x) {
    purrr::map(chains[[1]]$p_chain, ~ matrix(unlist(.x), nrow = 1)) %>%
      Reduce(rbind, .)
  }
)
mu_chains <- purrr::map(
  chains,
  function(x) {
    purrr::map(chains[[1]]$mu_chain, ~ matrix(unlist(.x), nrow = 1)) %>%
      Reduce(rbind, .)
  }
)
purrr::map(mu_chains, ~ mean(.x[10000:20000, 2]))

# NEED TO DERIVE THE TRUE POSTERIORS FOR EACH PARAMETER TO COMPARE FOR WASSERSTEIN DISTANCE