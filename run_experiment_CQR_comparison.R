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
  dt.names = "all",
  REP = 1:20,
  seed = 1,
  N.threads = N.threads,
  title = NULL,
  alpha = 0.1,
  test_frac = 0.2,
  rho = 0.8,
  tau_by = 0.005,
  max.n.tau = 50,
  funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F),
  method = "all",
  tail.error.control = "both",
  score.type = "both",
  tune_quantile_CV = T,
  coverage_factor = 1
)

#MORE SETTINGS
arg_list$dt.names = c(
  "bike", #Bike Sharing
  "bio", #Physicochemical Properties of Protein Tertiary Structure
  "blog", #BlogFeedback
  "community", #Communities and Crime
  "concrete", #Concrete Compressive Strength
  "facebook_1", #Facebook Comment Volume (Variant 1)
  "facebook_2", #Facebook Comment Volume (Variant 2)
  "homes", #House Sales in King County, USA
  "meps_19", #MEPS 2019
  "meps_20", #MEPS 2020
  "meps_21", #MEPS 2021
  "star", #STAR
  "super"#Superconductivity Data
)

settings = list(arg_list, arg_list) #list(arg_list, arg_list,....)
###setting: 1 (example)
settings[[1]]$tune_quantile_CV = F
###setting:2 (comparison of all CQR methods on all dataset)
settings[[2]]$coverage_factor = 1
for (sett in settings)
{
  do.call(experiment_CQR, args = sett)
}

#ONE SETTING
#settings = arg_list
#settings$dt.names = c("meps_19","meps_20")
#settings$coverage_factor = 1
#do.call(experiment_CQR, args = settings)




