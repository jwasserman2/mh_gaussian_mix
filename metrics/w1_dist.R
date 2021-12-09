library(doParallel)
library(foreach)
library(parallel)

#' Calculate the Wasserstein distance between two empirical distributions
#' Based off the scipy.stats.wasserstein_distance calculation, which is based
#' on the calculation devised by Bellemare et. al in "The Cramer Distance as a
#' Solution to Biased Wasserstein Gradients" (2017).
calculate_wasserstein_distance <- function(dist1, dist2, parallel = TRUE) {  
  dist1_sorted <- dist1[order(dist1)]
  dist2_sorted <- dist2[order(dist2)]
  all_sorted <- sort(c(dist1, dist2))
  deltas <- all_sorted[2:length(all_sorted)] - all_sorted[1:(length(all_sorted) -1)]
  
  if (parallel) {
    message("Sorting distributions in parallel on FORK cluster with 2 cores...")
    cl <- parallel::makeCluster(2, type="FORK")
    doParallel::registerDoParallel(cl)
    
    cdf_indices <- foreach::foreach(dist = list(dist1_sorted, dist2_sorted)) %dopar% {
      purrr::map_int(all_sorted[-length(all_sorted)], ~ length(dist[dist <= .x]))
    }
    parallel::stopCluster()

    dist1_cdf_indices <- cdf_indices[[1]]
    dist2_cdf_indices <- cdf_indices[[2]]
  } else {
    message("Sorting distributions sequentially...")
    dist1_cdf_indices <- purrr::map_int(
      all_sorted[-length(all_sorted)],
      ~ length(dist1_sorted[dist1_sorted <= .x])
    )
    dist2_cdf_indices <- purrr::map_int(
      all_sorted[-length(all_sorted)],
      ~ length(dist2_sorted[dist2_sorted <= .x])
    )
  }

  dist1_cdf <- dist1_cdf_indices / length(dist1)
  dist2_cdf <- dist2_cdf_indices / length(dist2)
  w1_dist <- sum(abs(dist1_cdf - dist2_cdf) * deltas)
  return(w1_dist)
}
