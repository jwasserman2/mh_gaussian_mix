library(foreach)
library(parallel)
library(doParallel)
library(jsonlite)

args <- commandArgs(TRUE)
SAMPLER <- args[1]
N_CHAINS <- as.numeric(args[2])
N_STEPS <- as.numeric(args[3])
N_DATA <- as.numeric(args[4])
DATA_DIST <- args[5]
TRUE_MEANS <- as.numeric(strsplit(args[6], ",")[[1]])
TRUE_VARIANCES <- as.numeric(strsplit(args[7], ",")[[1]])
TRUE_PROBS <- as.numeric(strsplit(args[8], ",")[[1]])
PROPOSAL_PARAM <- as.numeric(args[9])
SEED <- as.numeric(args[10])

stopifnot(all(length(TRUE_MEANS) == length(TRUE_VARIANCES), length(TRUE_MEANS) == length(TRUE_PROBS)))
N_COMPONENTS <- length(TRUE_MEANS)
TRUE_SIGMA <- diag(TRUE_VARIANCES)

source("./algorithms/constants.R")
source("./algorithms/data_sim.R")
algorithm_file <- switch(
  SAMPLER,
  "rw_metropolis" = "rw_metropolis",
  "mip_sampler" = "multiple_independent_proposals",
  "mis_sampler" = "metropolized_independent_sampler"
)
source(paste0("./algorithms/", algorithm_file, ".R"))

data_gen_func <- switch(
  DATA_DIST,
  "normal" = generate_normal_mixture_data,
  "t" = generate_t_mixture_data
)
data <- data_gen_func(N_DATA, TRUE_MEANS, sqrt(TRUE_VARIANCES), TRUE_PROBS, seed = 45496)

sampler_func <- switch(
  SAMPLER,
  "rw_metropolis" = rw_metropolis,
  "mis_sampler" = run_metropolized_independent_sampler,
  "mip_sampler" = run_multiple_independent_proposals
)
priors <- list(
  "alpha" = rep(1 / N_COMPONENTS, N_COMPONENTS),
  "mu_0" = rep(0, N_COMPONENTS),
  "tau_2" = 5
)

sampler_args <- list(
  "data" = data,
  "priors" = priors,
  "n_components" = N_COMPONENTS,
  "proposal_param" = PROPOSAL_PARAM,
  "n_iters" = N_STEPS,
  "true_sigma" = TRUE_SIGMA
)
if (SAMPLER == "mip_sampler") {
  sampler_args[["K"]] <- 5
}

n_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_cores, type = "FORK")
doParallel::registerDoParallel(cl)

chains <- foreach::foreach(seed = seq(SEED, SEED + N_CHAINS, by = 1)) %dopar% {
  sampler_args[["seed"]] <- seed
  do.call(sampler_func, sampler_args)
}
parallel::stopCluster(cl)

json <- jsonlite::toJSON(chains)
file_suffix <- paste(c(SAMPLER, DATA_DIST, N_CHAINS, N_STEPS, "hash", round(runif(1), 5) * 1e5),
                     collapse = "_")
json_file <- paste0("./jsons/", file_suffix)
jsonlite::write_json(json, json_file)
