#######run_experiment_DCP###########

#library(conformalInference)
library(tictoc)
library(tidyverse)
library(parallel)
library(fastDummies)
library(data.table)
library(ranger)

source("functions/new_models.R")
source("functions/experiment_DCP.R")

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
  taus = NULL,
  funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F),
  exact_coverage = T,
  tail.error.control = "both",
  score.type = "both"
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
  "super" #Superconductivity Data
)

# Se vuoi eseguire solo DCP e DCP_optim (senza le varianti _sep), imposta:
# arg_list$tail.error.control = "single"

#JUST ONE SETTING------
settings = arg_list
#settings$dt.names = c("meps_19", "meps_20")
do.call(experiment_DCP, args = settings)


#MULTIPLE SETTINGS----------
#settings = list(arg_list, arg_list) #list(arg_list, arg_list,....)
###setting: 1 (standard datasets)
#settings[[1]]$dt.names = c("bike", "community", "concrete", "homes", "star", "meps_19", "meps_20", "meps_21", "super")
###setting: 2 (large datasets with nodesize = 50)
#settings[[2]]$dt.names = c("bio", "blog", "facebook_1", "facebook_2")
#settings[[2]]$title = "Test"
#settings[[2]]$funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F, nodesize = 50)
#
#for (sett in settings)
#{
#  do.call(experiment_DCP, args = sett)
#}