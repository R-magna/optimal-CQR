#########EXPERIMENT_CQR_COMPARISON_functions###########
#write_AVG.files-----
#write average result CSV file
write_AVG.files = function(paths = NULL)
{
  require(dplyr)
  require(data.table)
  if (is.null(paths))
    paths = list.dirs("comparison_CQR/", recursive = F)
  paths = sub("/$", "", paths)
  
  for (path in paths)
  {
    type = strsplit(path, "/")[[1]][2]
    path_avg = paste0(path, "/AVG_result/")
    if (!dir.exists(path_avg))
      dir.create(path_avg, recursive = TRUE)
    for (path_dataset in setdiff(list.dirs(path, recursive = F), paste0(path, "/AVG_result")))
    {
      dt.name = strsplit(x = path_dataset, split = "/")[[1]] |>
        tail(1)
      path.file.avg = paste0(path_avg, dt.name, "._.avg", ".csv")
      if (!file.exists(path.file.avg))
      {
        tab_avg = c()
        for (file.name in list.files(path_dataset,
                                     pattern = "\\.csv$",
                                     recursive = F))
        {
          rep = ((strsplit(file.name, "rep_")[[1]][2]) |>
                   strsplit(split = "\\."))[[1]][1]
          dt = fread(paste0(path_dataset, "/", file.name))
          tab_avg = rbind(tab_avg,
                          cbind(
                            rep = rep,
                            select(
                              dt,
                              dataset,
                              type,
                              is.out.dx,
                              is.out.sx,
                              coverage,
                              width,
                              w.score
                            ) |>
                              group_by(dataset, type) |>
                              summarise_all(mean)
                          ))
        }
        fwrite(tab_avg, file = path.file.avg)
      }
    }
    #create finals tab
    path_final = paste0(path_avg, "final_tab/")
    if (!dir.exists(path_final))
      dir.create(path_final, recursive = TRUE)
    
    tab_avg = c()
    for (file.name in list.files(path_avg, pattern = "\\.csv$", recursive = F))
    {
      dt.name = strsplit(x = file.name, split = "\\._.")[[1]][1]
      dt = fread(paste0(path_avg, file.name))
      if (!"dataset" %in% colnames(dt))
        dt =  cbind(dataset = dt.name, dt)
      tab_avg = rbind(tab_avg, dt)
    }
    
    fwrite(tab_avg,
           file = paste0(path_final, type, "._.", "1.all_rep_result", ".csv"))
    tab_avg = tab_avg |>
      group_by(dataset, type) |>
      summarise(
        err.dx = mean(is.out.dx),
        err.sx = mean(is.out.sx),
        coverage = mean(coverage),
        width = mean(width),
        w.score = mean(w.score)
      )
    fwrite(tab_avg,
           file = paste0(path_final, type, "._.", "2.unique_dataset_result", ".csv"))
    tab_avg = tab_avg |>
      group_by(type) |>
      summarise(
        err.dx = mean(err.dx),
        err.sx = mean(err.sx),
        coverage = mean(coverage),
        width = mean(width),
        w.score = mean(w.score)
      )
    fwrite(tab_avg,
           file = paste0(path_final, type, "._.", "3.unique_method_result", ".csv"))
  }
}
#experiment-------
experiment_CQR = function(dt.names = NULL,
                      REP = 1:10,
                      seed = 1,
                      N.threads = max(1, parallel::detectCores() - 1),
                      title = NULL,
                      alpha = 0.1,
                      rho = 0.8,
                      tau_by = 0.005,
                      max.n.tau = 50,
                      funs = quant.rf_ranger.funs(nthreads = N.threads, importance = F),
                      funs_QR = quant.rf_ranger.funs(nthreads = N.threads,
                                                     importance = F,
                                                     keep.inbag = F),
                      method = "all",
                      tail.error.control = "both",
                      score.type = "both",
                      tune_quantile_CV = T,
                      coverage_factor = 0.85,
                      calculate.QR = F,
                      only.numeric_X = F,
                      std.var = T,
                      frac.y = T,
                      log.y = F,
                      factor_to_dummy = T)
{
  #require(conformalInference)
  require(tictoc)
  require(tidyverse)
  require(parallel)
  require(fastDummies)
  require(data.table)
  
  source("import_datasets.R")
  source("new_models.R")
  source("general_helper.R")
  source("CQR_helper.R")
  source("CQR.R")

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
    "star" #STAR
  )
  
  if (is.null(dt.names))
    dt.names = possible.dt
  
  Sys.setenv(R_RANGER_NUM_THREADS = N.threads)
  options(Ncpus = N.threads)
  options(readr.show_col_types = FALSE)
  
  #create folders where the results are saved
  if (!is.null(title))
    title = paste0("_", title)
  if (tune_quantile_CV) {
    path_base = paste0("comparison_CQR/cv_", coverage_factor, title, "/")
  } else {
    path_base = paste0("comparison_CQR/no_cv", title, "/")
  }
  if (!dir.exists(path_base)) {
    dir.create(path_base, recursive = TRUE)
  }
  
  for (dt.name in dt.names)
  {
    print(dt.name)
    dt = import_dataset(dt.name, factor_to_dummy = factor_to_dummy)
    X = as.data.table(dt$X)
    y = dt$Y
    rm(dt)
    if (only.numeric_X) {
      X = select_if(X, is.numeric)
    }
    if (log.y) {
      y = log(y + 1)
    }
    for (rep in REP)
    {
      print(rep)
      set.seed(seed + which(dt.name == possible.dt) + 100 * rep)
      i.tr = sample(1:nrow(X), size = nrow(X) * .8, replace = F)
      x.train = X[i.tr, ]
      y.train = y[i.tr]
      x.test = X[-i.tr, ]
      y.test = y[-i.tr]
      
      if (std.var)
      {
        x.train = scale(x.train, center = T, scale = T)
        x.test = scale(
          x.test,
          center = attr(x.train, "scaled:center"),
          scale = attr(x.train, "scaled:scale")
        ) |>
          as.data.table()
        x.train = as.data.table(x.train)
      }
      if (frac.y)
      {
        y.mu = mean(y.train)
        y.train = y.train / abs(y.mu)
        y.test = y.test / abs(y.mu)
      }
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
        verbose = F
      )
      toc()
      
      if (calculate.QR)
      {
        out.qr = funs_QR$train.fun(x = x.train,
                                   y = y.train,
                                   q = c(alpha / 2, 1 - alpha / 2))
        tic()
        out$Intervals$QR = list(
          lo = funs_QR$predict.fun(out.qr, newx = x.test, q = alpha / 2),
          up = funs_QR$predict.fun(out.qr, newx = x.test, q = 1 - alpha / 2)
        )
        toc()
      }
      new_res = rbindlist(lapply(names(out$Intervals), function(ty) {
        data.table(
          dataset = dt.name,
          rep = rep,
          type = ty,
          y = y.test,
          create_tab_out.cp(
            out.cp = out$Intervals[[ty]],
            y.true = y.test,
            alpha = alpha
          )
        )
      }))
      
      ##write rep result
      path_res = paste0(path_base, dt.name, "/")
      if (!dir.exists(path_res)) {
        dir.create(path_res, recursive = TRUE)
      }
      fwrite(new_res,
             file = paste0(path_res, dt.name, "._.", "rep_", rep, ".csv"))
      
      rm(new_res)
    }
  }
  #write avg result
  write_AVG.files(path_base)
}