# Generate data
N_COMPONENTS <- 2
N_CHAIN <- 1000
N_DATA <- 500
true_mu <- c(-1, 1)
true_sds <- rep(0.5, N_COMPONENTS)
true_p <- c(.25, .75)
data <- generate_normal_mixture_data(N_DATA, true_mu, true_sds, true_p)
t_data <- generate_t_mixture_data(N_DATA, true_mu, true_sigma, true_p)

# priors
priors <- list(
  "alpha" = rep(1, N_COMPONENTS),
  "k_0" = 1,
  "mu_0" = rep(0, N_COMPONENTS),
  "S_0" = diag(N_COMPONENTS) * 5,
  "v_0" = N_COMPONENTS + 1)
proposal_alpha <- 5
proposal_sd <- 0.1
proposal_df <- 100
max_proposal_sigma <- 10
proposal_corr_0 <- 100
true_sigma <- diag(true_sds^2)

posterior <- suppressMessages(rw_metropolis(data, priors, N_COMPONENTS, proposal_alpha, proposal_sd, 
                           max_proposal_sigma, proposal_corr_0, N_CHAIN,
                           true_sigma, seed = 2045))
sum(posterior$acceptances) / length(posterior$acceptances)
setNames(data.frame(posterior$p_chain), paste0("p_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("p_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
setNames(data.frame(posterior$mu_chain), paste0("mu_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("mu_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
purrr::map_dfr(posterior$sigma_chain, ~ data.frame(t(sqrt(diag(.x))))) %>%
  setNames(paste0("sigma_", 1:N_COMPONENTS)) %>%
  cbind(., "step" = 1:N_CHAIN) %>%
  tidyr::pivot_longer(cols = paste0("sigma_", 1:N_COMPONENTS)) %>%
  ggplot() +
  geom_line(aes(x = step, y = value, group = name, color = name))
mis_posterior <- run_metropolized_independent_sampler(data, priors, N_COMPONENTS,
                                                      proposal_alpha, proposal_sd,
                                                      proposal_df, N_CHAIN, seed = 2045)
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
# NEED TO DERIVE THE TRUE POSTERIORS FOR EACH PARAMETER TO COMPARE FOR WASSERSTEIN DISTANCE