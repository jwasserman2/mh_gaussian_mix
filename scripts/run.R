library(doParallel)
library(foreach)
library(parallel)

args <- commandArgs(TRUE)
SAMPLER <- args[1]
print(SAMPLER)