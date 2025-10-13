#######run_experiment_CQR_comparison###########

#library(conformalInference)
library(tictoc)
library(tidyverse)
library(parallel)
library(fastDummies)
library(data.table)
library(ranger)

source("functions/new_models.R")
source("functions/experiment_CQR_comparison.R")

N.threads = max(1, parallel::detectCores() - 1)

settings = list()
arg_list = list(
  dt.names = NULL,
  REP = 1:10,
  seed = 1,
  N.threads = N.threads,
  title = NULL,
  alpha = 0.1,
  rho = 0.8,
  tau_by = 0.005,
  max.n.tau = 50,
  funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F),
  funs_QR = quant.rf_ranger.funs(nthreads = N.threads, importance = F, keep.inbag = F),
  method = "all",
  tail.error.control = "both",
  score.type = "both",
  tune_quantile_CV = T,
  coverage_factor = 0.85,
  calculate.QR = T,
  only.numeric_X = F,
  std.var = T,
  frac.y = T,
  log.y = F,
  factor_to_dummy = T
)

settings = list(arg_list) #list(arg_list, arg_list,....)
#setting: 1
settings[[1]]$dt.names = "community"
settings[[1]]$title = "example"
settings[[1]]$method = "all"
settings[[1]]$tail.error.control = "single"
settings[[1]]$tune_quantile_CV = F
#setting:2....

for (sett in settings)
{
  do.call(experiment_CQR, args = sett)
}
