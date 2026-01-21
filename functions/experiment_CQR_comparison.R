#########EXPERIMENT_CQR_COMPARISON_functions###########
#experiment-------
experiment_CQR = function(dt.names = NULL,
                      REP = 1:10,
                      seed = 1,
                      test_frac = 0.2,
                      N.threads = max(1, parallel::detectCores() - 1),
                      title = NULL,
                      alpha = 0.1,
                      rho = 0.8,
                      tau_by = 0.005,
                      max.n.tau = 50,
                      funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F),
                      method = "all",
                      tail.error.control = "both",
                      score.type = "both",
                      tune_quantile_CV = T,
                      coverage_factor = 0.85)
{
  #require(conformalInference)
  require(tictoc)
  require(tidyverse)
  require(parallel)
  require(fastDummies)
  require(data.table)
  
  source("functions/import_datasets.R")
  source("functions/new_models.R")
  source("functions/general_helper.R")
  source("functions/CQR_helper.R")
  source("functions/CQR.R")

  possible.dt = c(
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
  
  if (any(dt.names=="all"))
    dt.names = possible.dt
  
  Sys.setenv(R_RANGER_NUM_THREADS = N.threads)
  options(Ncpus = N.threads)
  options(readr.show_col_types = FALSE)
  
  #create folders where the results are saved
  if (!is.null(title))
    title = paste0("_", title)
  if (tune_quantile_CV) {
    path_base = paste0("comparison_CQR/q.cv_", coverage_factor, title, "/")
  } else {
    path_base = paste0("comparison_CQR/q.no.cv", title, "/")
  }
  if (!dir.exists(path_base)) {
    dir.create(path_base, recursive = TRUE)
  }
  
  for (dt.name in dt.names)
  {
    print(dt.name)
    dt = import_dataset(dt.name)
    X = as.data.table(dt$X)
    y = dt$Y
    rm(dt)
    
    for (rep in REP)
    {
      print(rep)
      #split
      set.seed(seed + which(dt.name == possible.dt) + 100 * rep)
      i.tr = sample(1:nrow(X),
                    size = floor(nrow(X) * (1 - test_frac)),
                    replace = F)
      x.train = X[i.tr, ]
      y.train = y[i.tr]
      x.test = X[-i.tr, ]
      y.test = y[-i.tr]
      
      #proper training set
      I1 = sample(1:nrow(x.train),
                  size = floor(nrow(x.train) * rho),
                  replace = F)
      
      #standardization
      std.stat = scale(x.train[I1,], center = T, scale = T)
      
      x.test = scale(
        x.test,
        center = attr(std.stat, "scaled:center"),
        scale = attr(std.stat, "scaled:scale")
      )|>
        as.data.table()
      
      x.train = scale(
        x.train,
        center = attr(std.stat, "scaled:center"),
        scale = attr(std.stat, "scaled:scale")
      )|>
        as.data.table()
      #trasform y
      y.mu_abs = mean(y.train[I1])|>
        abs()
      y.train = y.train / y.mu_abs
      y.test = y.test / y.mu_abs
      
      rm(std.stat, y.mu_abs)
      
      if(dt.name=="super")
      {
        y.min = 0
      } else
      {
        y.min = NULL
      }
      
      set.seed(seed)
      tic()
      out = conformal.quant.split(
        x = x.train,
        y = y.train,
        x0 =  x.test,
        train.quant.fun = funs$train.fun,
        predict.quant.fun = funs$predict.fun,
        predict.quant.bootstrap.fun = funs$predict.bootstrap.fun,
        alpha = alpha,
        tau_by = tau_by,
        max.n.tau = max.n.tau,
        rho = rho,
        method = method,
        tail.error.control = tail.error.control,
        score.type = score.type,
        tune_quantile_CV = tune_quantile_CV,
        coverage_factor = coverage_factor,
        exact_coverage = T,
        verbose = F,
        split = I1
      )
      toc()
      
      new_res = rbindlist(lapply(names(out$Intervals), function(ty) {
        data.table(
          dataset = dt.name,
          rep = rep,
          type = ty,
          y = y.test,
          create_tab_out.cp(
            out.cp = out$Intervals[[ty]],
            y.true = y.test,
            alpha = alpha,
            y_min = y.min
          )
        )
      }))
      
      #write rep result
      path_res = paste0(path_base, dt.name, "/")
      if (!dir.exists(path_res)) {
        dir.create(path_res, recursive = TRUE)
      }
      fwrite(new_res,
             file = paste0(path_res, dt.name, "._.", "rep_", rep, ".csv"))
      
      rm(new_res, out)
    }
  }
  #write avg result
  write_AVG.files(path_base)
}

