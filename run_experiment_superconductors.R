#######run_experiment_superconductors#####

library(parallel)
library(plyr)
library(ranger)

source("functions/experiment_superconductors.R")

N.threads = max(1, parallel::detectCores() - 1)
Sys.setenv(R_RANGER_NUM_THREADS = N.threads)
options(Ncpus = N.threads)

arg_list = list(
  rep.vect = 1:10,
  alpha = 0.1,
  rho = 0.8,
  models = "all",
  model.selection = "stab",#c("PCA,"stab")
  seed = 1,
  frac.test = 0.2,
  trasf.log.y = F,
  trasf.add.group = T,
  coverage_factor = 0.85,
  CV_hyperpar = T,#T,"load" (F not implemented)
  save_best.par = T,
  N.threads = N.threads
)

do.call(experiment_super, args = arg_list)




