source("./algorithms/data_sim.R")
#' Plot an MCMC chain against a Gaussian mixture distribution (currently only for one parameter of interest)
plot_mcmc_vs_true_density <- function(chain, true_density_kwargs) {
  if (any(names(true_density_kwargs) == "df")) {
    true_density_func <- calculate_heavytail_mixture_prob
  } else {
    true_density_func <- calculate_gaussian_mixture_prob
  }
  
  n_iters <- dim(chain)[1]
  chain_df <- data.frame("chain" = chain)
  mu <- rlang::`%||%`(true_density_kwargs[["noncentral_param"]], true_density_kwargs[["mu"]])
  sigma <- rlang::`%||%`(true_density_kwargs[["df"]], true_density_kwargs[["sigma"]])
  density_breaks <- seq(min(mu) - N_SDS_TO_PLOT * max(sqrt(sigma)), max(mu) + N_SDS_TO_PLOT * max(sqrt(sigma)), 0.1)
  true_df <- data.frame("x" = density_breaks,
                        "true_density" = purrr::map_dbl(
                          density_breaks,
                          function(x) {
                            run_kwargs <- true_density_kwargs
                            run_kwargs[['x']] <- x
                            do.call(true_density_func, run_kwargs)
                          }
                        ))
  x_plot_breaks <- round(c(unique(c(mu, 0)), min(mu) - N_SDS_TO_PLOT * max(sqrt(sigma)),
                           min(mu) - N_SDS_TO_PLOT / 2 * round(max(sqrt(sigma)), 1), min(mu),
                           max(mu) + N_SDS_TO_PLOT / 2 * round(max(sqrt(sigma)), 1),
                           max(mu) + N_SDS_TO_PLOT * round(max(sqrt(sigma)), 1)))

  ggplot2::ggplot(chain_df) +
    ggplot2::geom_histogram(ggplot2::aes(x = chain), bins = 30) +
    ggplot2::geom_line(data = true_df, ggplot2::aes(x = x, y = true_density * n_iters)) +
    ggplot2::scale_x_continuous(breaks = x_plot_breaks) +
    ggplot2::scale_y_continuous(labels = ~ .x / n_iters) +
    ggplot2::labs(x = expression(theta), y = "Density") +
    ggplot2::theme(panel.background = element_blank(),
                   axis.text = element_text(size = 10),
                   axis.ticks.length = unit(0.05, "inches"))
}
