#######run_experiment_superconductors#####

library(parallel)
library(plyr)
library(ranger)

source("functions/exmperiment_superconductors.R")

N.threads = max(1, parallel::detectCores() - 1)
Sys.setenv(R_RANGER_NUM_THREADS = N.threads)
options(Ncpus = N.threads)

settings = list()
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
  CV_hyperpar = T,#T,F,"load"
  save_best.par = T,
  N.threads = N.threads
)

experiment_super(rep.vect = 1:10, N.threads = N.threads) #list(arg_list, arg_list,....)
#setting: 1
settings[[1]]$dt.names = "community"
settings[[1]]$title = "example"
settings[[1]]$method = "all"
settings[[1]]$tail.error.control = "single"
settings[[1]]$tune_quantile_CV = F
#setting:2....

for (sett in settings)
{
  do.call(experiment_super, args = sett)
}