#######exmperiment_superconductors#########

#experiment_super-----
experiment_super = function(rep.vect = 1:10,
                            alpha = 0.1,
                            rho = 0.8,
                            models = "all",
                            model.selection = "stab", #c("PCA,"stab")
                            seed = 1,
                            frac.test = 0.2,
                            trasf.log.y = F,
                            trasf.add.group = T,
                            coverage_factor = 0.85,
                            CV_hyperpar = "load", #T,F
                            save_best.par = T,
                            N.threads = 1)
{
  #library, source and set----
  require(readr)
  require(tidyverse)
  require(devtools)
  require(conformalInference)#devtools::install_github(repo="ryantibs/conformal", subdir="conformalInference")
  require(LCP)#devtools::install_github("LeyingGuan/LCP")
  require(tictoc)
  require(parallel)
  require(fastDummies)
  require(data.table)
  #modelling and BO
  require(ranger)
  require(quantregForest)
  require(quantreg)
  require(xgboost)
  require(ranger)
  require(mlr)
  require(mlrMBO)
  require(glmnet)
  require(kknn)
  require(e1071)
  dist = stats::dist#LCP::dist
  
  source("general_helper.R")
  source("new_models.R")
  source("CQR_helper.R")
  source("CQR.R")
  source("DCP_function.R")
  source("LCP_function.R")
  source("stability_selection_lasso.R")
  source("BO_tuning.R")
  
  if (any(models == "all"))
    models = c("XGB", "KNN", "RDG")#,"SVM")
  models = append(models, "RF") |>
    unique()
  models.var.selection = c("KNN", "EN", "SVM", "RDG", "RF.sel", "XGB.sel")
  
  if (!dir.exists("result.superconductor")) {
    dir.create("result.superconductor", recursive = TRUE)
  }
  
  #import data----
  molecules = read_csv("superconductivty+data/unique_m.csv")
  #data = fetch_ucirepo(name = "Superconductivty Data")$data$original #a volte non si connette al server UCI
  data = read_csv("superconductivty+data/train.csv")
  #all_molecules = molecules$material
  
  if (trasf.add.group)
  {
    data = mutate(
      data,
      G_Iron = as.numeric(molecules$Fe > 0),
      G_Cuprate = as.numeric(molecules$Cu > 0 &
                               molecules$O > 0)
    ) |>
      #G_MgB = as.numeric(molecules$Mg>0 & molecules$B>0))|>
      relocate(critical_temp, .after = last_col())
  }
  rm(molecules)
  name.dir = paste0(
    "cov",
    1 - alpha,
    "_cv.factor",
    coverage_factor,
    ifelse(trasf.log.y, "_log", ""),
    ifelse(trasf.add.group, "_Group", "")
  )
  if (!dir.exists(paste0("result.superconductor/", name.dir))) {
    dir.create(paste0("result.superconductor/", name.dir), recursive = TRUE)
  }
  #rep----
  for (rep in rep.vect)
  {
    print(rep)
    set.seed(seed + rep)
    #data split----
    i.test = sample(1:nrow(data), size = floor(0.2 * nrow(data)))
    
    test.y = data[i.test, ]$critical_temp
    test.x = select(data[i.test, ], -critical_temp)
    train.x = select(data[-i.test, ], -critical_temp)
    train.y = data[-i.test, ]$critical_temp
    train.x.origin = train.x
    
    n = nrow(train.x)
    p = ncol(train.x)
    #proper training set
    i1 = sample(seq_len(n), floor(n * rho))
    
    dt = scale(train.x[i1, ], center = T, scale = T)
    
    ts.std = scale(
      test.x,
      center = attr(dt, "scaled:center"),
      scale = attr(dt, "scaled:scale")
    ) |>
      as_tibble()
    
    tr.std = train.x
    tr.std[-i1, ] = scale(
      tr.std[-i1, ],
      center = attr(dt, "scaled:center"),
      scale = attr(dt, "scaled:scale")
    )
    tr.std[i1, ] = dt
    
    y_mu = mean(train.y[i1])
    y_sd = sqrt(var(train.y[i1]))
    tr.y.std = (train.y - y_mu) / y_sd
    ts.y.std = (test.y - y_mu) / y_sd
    
    rm(dt)
    #trasfy log y
    if (trasf.log.y) {
      train.y = log(train.y)
      test.y = log(test.y)
    }
    
    #stability selection----
    if ("stab" %in% model.selection)
    {
      stab = f_stability.selection_lasso(
        x.train =  as_tibble(tr.std[i1, ]),
        y.train = train.y[i1],
        B = 50,
        n.lambda = 200,
        lambda.min.ratio = 1e-3
      )
      
      pi_j.max_L = select(stab, -lambda, -err, -num.selected) |>
        apply(2, max) |>
        sort(decreasing = T)
      
      #definition of possible thresholds
      tau = c(max(min(pi_j.max_L), 0.3), 0.6, 0.75, 0.9) |> #max(pi_j.max_L))|>
        sort(decreasing = T)
      p.uniq = pi_j.max_L |>
        unique()
      tau.new = tau[1]
      #the thresholds are defined, avoiding the possibility that they select
      #the same number of variables
      for (i in 2:length(tau))
      {
        t = tau[i]
        if (any(p.uniq >= t & p.uniq < tau.new[i - 1])) {
          tau.new = append(tau.new, t)
        } else {
          if (i < length(tau))
          {
            t = median(p.uniq[p.uniq < tau.new[i - 1] & p.uniq > tau[i + 1]])
            if (!is.na(t))
              tau.new = append(tau.new, t)
          }
        }
      }
      tau = tau.new
      #variables selected for each tau
      var_tau.stability = sapply(tau, function(t) {
        names(pi_j.max_L[pi_j.max_L >= t])
      })
      names(var_tau.stability) = paste0("t_", tau)
      
      rm(tau.new, t, p.uniq, i, pi_j.max_L, tau)
    }
    #PCA-----
    if ("PCA" %in% model.selection)
    {
      pca_train = prcomp(as_tibble(tr.std[i1, ]),
                         center = F,
                         scale. = F)
      
      tr.pca = tr.std
      tr.pca[-i1, ] = predict(pca_train, newdata = tr.std[-i1, ])
      tr.pca[i1, ] = pca_train$x
      ts.pca = predict(pca_train, newdata = ts.std)
      
      cum_var = cumsum(pca_train$sdev^2) / sum(pca_train$sdev^2) * 100
      
      tau = c(0.85, 0.9, 0.95, 0.975) * 100
      
      var_tau.pca = sapply(tau, function(t) {
        colnames(pca_train$x)[cum_var <= t]
      })
      names(var_tau.pca) = paste0("t_", tau)
      rm(tau, cum_var)
    }
    
    #best_par----
    #load already tuned
    if (is.null(CV_hyperpar))
      CV_hyperpar = F
    
    if (CV_hyperpar == "load")
    {
      best_par = readRDS(
        paste0(
          "result.superconductor/",
          name.dir,
          "/best_par/best_par_rep_",
          rep,
          ".rds"
        )
      )
    } else if (CV_hyperpar == T)
    {
      K.cv = 10
      best_par = list()
      param.set_list = list()
      
      #HP def----
      #HP RF
      p.rf = ceiling(sqrt(p))
      param.set_list[["RF"]] = list(
        mod.name = "regr.ranger",
        par.set = makeParamSet(
          makeIntegerParam(
            "mtry",
            lower = max(2, p.rf - 10),
            upper = min(p.rf + 10, p)
          ),
          makeIntegerParam("min.node.size", lower = 1, upper = 25)
          #makeNumericParam("sample.fraction", lower = 0.5, upper = 1)#,
          #    makeIntegerParam("max.depth", lower = 0, upper = 30)
        ),
        par.fix = list(
          "num.trees" = 500,
          "num.threads" = N.threads,
          "importance" = "none",
          "keep.inbag" = F,
          "sample.fraction" = 1
        )
      )
      param.set_list[["RF.sel"]] =  param.set_list[["RF"]]
      param.set_list[["RF.sel"]]$par.fix$num.trees = 500
      #param.set_list[["RF.base"]] = list(
      #  mod.name = "regr.randomForest",
      #  par.set = makeParamSet(
      #    makeIntegerParam("mtry", lower = max(2, p.rf - 5), upper = min(p.rf + 8, p)),
      #    makeIntegerParam("nodesize", lower = 1, upper = 25),
      #    makeNumericParam("sample.fraction", lower = 0.5, upper = 1)
      #  ),
      #  par.fix = list("ntree" = 500,"importance" = F, "keep.inbag" = F)
      #)
      rm(p.rf)
      #hp XGB
      param.set_list[["XGB"]] = list(
        mod.name = "regr.xgboost",
        par.set = makeParamSet(
          makeNumericParam("eta", lower = 0.01, upper = 0.5),
          makeIntegerParam("max_depth", lower = 3, upper = 20),
          makeIntegerParam("nrounds", lower = 50, upper = 500),
          makeIntegerParam("min_child_weight", lower = 1, upper = 10),
          makeNumericParam("colsample_bytree", lower = 0.4, upper = 1.0),
          makeNumericParam("gamma", lower = 0, upper = 5)
        ),
        par.fix = list("nthread" = N.threads)
      )
      param.set_list[["XGB.sel"]] =  param.set_list[["XGB"]]
      #HP glmnet
      param.set_list[["EN"]] = list(
        mod.name = "regr.cvglmnet",
        par.set = makeParamSet(
          makeNumericParam("alpha", lower = 0.1, upper = 1),
          makeDiscreteParam("s", values = c("lambda.min", "lambda.1se"))
        ),
        par.fix = NULL
      )
      #HP ridge
      param.set_list[["RDG"]] = list(
        mod.name = "regr.glmnet",
        par.set = makeParamSet(
          makeNumericParam(
            "lambda",
            lower = log(1e-4),
            upper = log(1e+2),
            trafo = exp
          )
        ),
        par.fix = list("alpha" = 0)
      )
      
      #param.set_list[["EN"]] = list(
      #  mod.name = "regr.glmnet",
      #  par.set = makeParamSet(
      #    makeNumericParam("alpha", lower = 0.2, upper = 1),
      #    makeNumericParam("lambda", lower = 10^-4, upper = 10)
      #  ),
      #  par.fix = list("nthread" = N.threads)
      #)
      
      #HP KNN
      param.set_list[["KNN"]] = list(
        mod.name = "regr.kknn",
        par.set = makeParamSet(
          makeIntegerParam("k", lower = 1, upper = 50),
          makeIntegerParam("distance", lower = 1, upper = 5),
          makeDiscreteParam(
            "kernel",
            values =  c("rectangular", "gaussian", "cos", "epanechnikov")
          )
        ),
        par.fix = list(scale = F)
      )
      
      #HP SVM
      param.set_list[["SVM"]] = list(
        mod.name = "regr.svm",
        par.set = makeParamSet(
          makeDiscreteParam("kernel", values = c("radial", "polynomial", "linear")),
          makeNumericParam(
            "cost",
            -2,
            2,
            trafo = function(x)
              10^x
          ),
          makeNumericParam(
            "epsilon",
            lower = -3,
            upper = -0.7,
            trafo = function(x)
              10^x
          ),
          makeNumericParam(
            "gamma",
            -12,
            2,
            trafo = function(x)
              2^x,
            requires = quote(kernel %in% c("radial", "polynomial", "sigmoid"))
          ),
          makeNumericParam(
            "coef0",
            lower = -0.5,
            upper = 0.5,
            requires = quote(kernel %in% c("polynomial", "sigmoid"))
          )
        ),
        par.fix = list(scale = F, degree = 2)
      )
      
      #CV----
      param.set_list = param.set_list[names(param.set_list) %in% models]
      
      mod.x.tr = train.x[i1, ]
      mod.y.tr = train.y[i1]
      for (new.mod in intersect(names(param.set_list), c("RF", "XGB"))) {
        res = f_mlrMBO.bayesian.hyperpar.tuning(
          data_mod = mod.x.tr,
          y = mod.y.tr,
          mod = param.set_list[[new.mod]][["mod.name"]],
          param.SET = param.set_list[[new.mod]][["par.set"]],
          par.fix = param.set_list[[new.mod]][["par.fix"]],
          k.fold.cv = K.cv
        )
        best_par[[new.mod]] = append(res$x, param.set_list[[new.mod]][["par.fix"]])
      }
      
      par_selez.var = list()
      err_tab = tibble(
        mod = character(),
        # tipi dichiarati
        type = character(),
        tau = character(),
        err = numeric()
      )
      
      if (any(models.var.selection %in% names(param.set_list)))
      {
        for (new.mod in intersect(names(param.set_list), models.var.selection))
        {
          if (new.mod == "SVM")
          {
            mod.y.tr = tr.y.std[i1]
          }
          else
          {
            mod.y.tr = train.y[i1]
          }
          for (ty in model.selection)
          {
            if (ty == "PCA") {
              tau = var_tau.pca
              mod.x.tr = pca_train$x |>
                as_tibble()
            } else {
              tau = var_tau.stability
              mod.x.tr = tr.std[i1, ]
            }
            for (t in names(tau))
            {
              if (new.mod == "RF.sel") {
                p.p = ncol(mod.x.tr[, tau[[t]]])
                p.rf = ceiling(sqrt(p.p))
                param.set_list[["RF.sel"]][["par.set"]]$pars$mtry$lower = max(2, p.rf - 5)
                param.set_list[["RF.sel"]][["par.set"]]$pars$mtry$upper = min(p.rf + 8, p.p)
              }
              res = f_mlrMBO.bayesian.hyperpar.tuning(
                data_mod = mod.x.tr[, tau[[t]]],
                y = mod.y.tr,
                mod = param.set_list[[new.mod]][["mod.name"]],
                param.SET = param.set_list[[new.mod]][["par.set"]],
                par.fix = param.set_list[[new.mod]][["par.fix"]],
                k.fold.cv = K.cv,
                time.budget = 3600 / 2
              )
              
              par_selez.var[[new.mod]][[ty]][[t]] = append(res$x, param.set_list[[new.mod]][["par.fix"]])
              
              err_tab = add_row(
                err_tab,
                mod = new.mod,
                type = ty,
                tau = t,
                err = res$y
              )
            }
          }
        }
        
        err_tab.star = err_tab |>
          group_by(mod) |>
          arrange(desc(err),
                  type,
                  if_else(type == "PCA", tau, NA_character_),
                  desc(if_else(type == "stab", tau, NA_character_)),) |>
          slice_min(err, n = 1, with_ties = FALSE) |>
          ungroup()
        for (i in 1:nrow(err_tab.star))
        {
          j = err_tab.star[i, ] |>
            select(-err) |>
            unlist()
          new.mod = j[["mod"]]
          ty = j[["type"]]
          t = j[["tau"]]
          
          best_par[[new.mod]] = c(par_selez.var[[new.mod]][[ty]][[t]], list(type = ty, tau = t))
        }
      }
      rm(res, mod.x.tr, mod.y.tr, new.mod, t, tau, ty, j, i)
      if (save_best.par)
        saveRDS(best_par,
                file = paste0("result.superconductor/", name.dir, "/best_par.rds"))
    } else if (CV_hyperpar == F)
    {
      #best_par = par.default
    }
    res.CP = list()
    #benchmark-----
    rmse = 9.5
    mod.tr = select(train.x, -starts_with("G_"))
    mod.ts =  select(test.x, -starts_with("G_"))
    mod = xgboost(
      data = xgb.DMatrix(data = as.matrix(mod.tr), label = train.y),
      nrounds = 750,
      eta = 0.02,
      max_depth = 16,
      min_child_weight = 1,
      colsample_bytree = 0.5,
      subsample = 0.5
    )
    dt = xgb.DMatrix(data = as.matrix(mod.ts))
    pred = predict(mod, dt)
    pred = tibble(lo = pred - rmse, up = pred + rmse)
    create_tab_out.cp(pred,
                      y.true = test.y,
                      y_min = 0,
                      alpha = alpha) |>
      colMeans()
    res.CP$benchmark = create_tab_out.cp(pred,
                                         y.true = test.y,
                                         y_min = 0,
                                         alpha = alpha)
    
    rm(dt, pred, rmse, mod)
    #CP-----
    funs.rf = rf_ranger.funs(
      nthreads = N.threads,
      keep.inbag = F,
      importance = F,
      ntree = best_par$RF[["num.trees"]],
      replace = ifelse(best_par$RF[["sample.fraction"]] <
                         1, F, T),
      nodesize = best_par$RF[["min.node.size"]],
      varfrac = best_par$RF[["mtry"]] / p
    )
    
    if ("XGB" %in% models)
    {
      funs.xgb = xgb.funs(
        nrounds = best_par$XGB[["nrounds"]],
        eta = best_par$XGB[["eta"]],
        max_depth =  best_par$XGB[["max_depth"]],
        min_child_weight = best_par$XGB[["min_child_weight"]],
        varfrac = best_par$XGB[["colsample_bytree"]],
        nthreads = N.threads
      )
    }
    if ("KNN" %in% models)
    {
      funs.knn = kknn.funs(
        k = best_par$KNN[["k"]],
        distance = best_par$KNN[["distance"]],
        kernel =  best_par$KNN[["kernel"]],
        scale = best_par$KNN[["scale"]]
      )
    }
    if ("RF.sel" %in% models)
    {
      if (best_par$RF.sel$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RF.sel$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RF.sel$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      p.mod = ncol(mod.x)
      funs.rf.sel = rf_ranger.funs(
        nthreads = N.threads,
        keep.inbag = F,
        importance = F,
        ntree = best_par$RF.sel[["num.trees"]],
        replace = ifelse(best_par$RF.sel[["sample.fraction"]] <
                           1, F, T),
        nodesize = best_par$RF.sel[["min.node.size"]],
        varfrac = best_par$RF.sel[["mtry"]] /
          p.mod
      )
    }
    if ("RDG" %in% models)
    {
      funs.rdg = ridge.funs(lambda = best_par[["RDG"]][["lambda"]], standardize = F)
    }
    if ("SVM" %in% models)
    {
      funs.svm = svm.funs(
        kernel = best_par[["SVM"]][["kernel"]],
        cost =  best_par[["SVM"]][["cost"]],
        epsilon = best_par[["SVM"]][["epsilon"]],
        gamma = best_par[["SVM"]][["gamma"]],
        degree = best_par[["SVM"]][["degree"]],
        coef0 = best_par[["SVM"]][["coef0"]],
        scale = F
      )
    }
    #best_par$EN
    #funs.en = elastic.funs(cv=T)
    #SCP-----
    #slit CP: RF
    out.cp = conformal.pred.split(
      x = train.x,
      y = train.y,
      x0 = as.matrix(test.x),
      alpha = alpha,
      rho = rho,
      train.fun = funs.rf$train,
      predict.fun = funs.rf$predict,
      seed = seed,
      verbose = T,
      split = i1
    )
    res.CP$SCP$RF = create_tab_out.cp(out.cp,
                                      y.true = test.y,
                                      y_min = 0,
                                      alpha = alpha)
    
    #Split CP: XGB
    if ("XGB" %in% models)
    {
      out.cp = conformal.pred.split(
        train.x,
        train.y,
        x0 = as.matrix(test.x),
        alpha = alpha,
        rho = rho,
        train.fun = funs.xgb$train,
        predict.fun = funs.xgb$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      res.CP$SCP$XGB = create_tab_out.cp(out.cp,
                                         y.true = test.y,
                                         y_min = 0,
                                         alpha = alpha)
    }
    #Split CP: KNN
    if ("KNN" %in% models)
    {
      if (best_par$KNN$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$KNN$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$KNN$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = conformal.pred.split(
        mod.x,
        train.y,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.knn$train,
        predict.fun = funs.knn$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      res.CP$SCP[["KNN"]] = create_tab_out.cp(out.cp,
                                              y.true = test.y,
                                              y_min = 0,
                                              alpha = alpha)
    }
    #Split CP: RDG
    if ("RDG" %in% models)
    {
      if (best_par$RDG$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RDG$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RDG$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = conformal.pred.split(
        mod.x,
        train.y,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.rdg$train,
        predict.fun = funs.rdg$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      res.CP$SCP[["RDG"]] = create_tab_out.cp(out.cp,
                                              y.true = test.y,
                                              y_min = 0,
                                              alpha = alpha)
    }
    #Split CP: SVM
    if ("SVM" %in% models)
    {
      if (best_par$SVM$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$SVM$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$SVM$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      mod.y.tr = tr.y.std
      
      out.cp = conformal.pred.split(
        mod.x,
        mod.y.tr,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.svm$train,
        predict.fun = funs.svm$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      out.cp$lo = out.cp$lo * y_sd + y_mu
      out.cp$up = out.cp$up * y_sd + y_mu
      res.CP$SCP[["SVM"]] = create_tab_out.cp(out.cp,
                                              y.true = test.y,
                                              y_min = 0,
                                              alpha = alpha)
    }
    #Split CP: elastic net
    
    #out.cp = conformal.pred.split(train.x, train.y, x0 = as.matrix(test.x),
    #                              alpha=alpha ,
    #                              rho = rho,
    #                              train.fun=funs$train,
    #                              predict.fun=funs$predict,
    #                              seed = seed, verbose = T, split = i1)
    #res.CP$SCP$EN$out=out.cp
    #res.CP$SCP$EN = create_tab_out.cp(out.cp,y.true = test.y, y_min = 0, alpha = alpha)
    
    #LACP----
    #split Locally Adaprive cp: xgb.rf
    funs.res = funs.rf
    if ("XGB" %in% models)
    {
      out.cp = conformal.pred.split(
        train.x,
        train.y,
        x0 = as.matrix(test.x),
        alpha = alpha,
        rho = rho,
        train.fun = funs.xgb$train,
        predict.fun = funs.xgb$predict,
        mad.train.fun = funs.res$train,
        mad.predict.fun = funs.res$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      #res.CP$LACP$XGB.rf$out=out.cp
      res.CP$LACP$XGB.rf = create_tab_out.cp(out.cp,
                                             y.true = test.y,
                                             y_min = 0,
                                             alpha = alpha)
    }
    #probabilmente è modello troppo complesso e lascia poca informazione
    
    #split Locally Adaprive cp: rf.rf
    funs.res = funs.rf
    out.cp = conformal.pred.split(
      train.x,
      train.y,
      x0 = as.matrix(test.x),
      alpha = alpha,
      rho = rho,
      train.fun = funs.rf$train,
      predict.fun = funs.rf$predict,
      mad.train.fun = funs.res$train,
      mad.predict.fun = funs.res$predict,
      seed = seed,
      verbose = T,
      split = i1
    )
    #res.CP$LACP$RF.rf$out=out.cp
    res.CP$LACP$RF.rf = create_tab_out.cp(out.cp,
                                          y.true = test.y,
                                          y_min = 0,
                                          alpha = alpha)
    
    #kn.rf.sel
    if ("KNN" %in% models)
    {
      if (best_par$KNN$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$KNN$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$KNN$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = conformal.pred.split(
        mod.x,
        train.y,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.knn$train,
        predict.fun = funs.knn$predict,
        mad.train.fun = funs.res$train,
        mad.predict.fun = funs.res$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      #res.CP$LACP$RF.rf$out=out.cp
      res.CP$LACP$KNN.rf = create_tab_out.cp(out.cp,
                                             y.true = test.y,
                                             alpha = alpha,
                                             y_min = 0)
      
      if ("RF.sel" %in% models)
      {
        out.cp = conformal.pred.split(
          mod.x,
          train.y,
          x0 = as.matrix(mod.x.ts),
          alpha = alpha,
          rho = rho,
          train.fun = funs.knn$train,
          predict.fun = funs.knn$predict,
          mad.train.fun = funs.rf.sel$train,
          mad.predict.fun = funs.rf.sel$predict,
          seed = seed,
          verbose = T,
          split = i1
        )
        #res.CP$LACP$RF.rf$out=out.cp
        res.CP$LACP$KNN.rf.sel = create_tab_out.cp(out.cp,
                                                   y.true = test.y,
                                                   alpha = alpha,
                                                   y_min = 0)
      }
      
    }
    if ("RDG" %in% models)
    {
      if (best_par$RDG$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RDG$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RDG$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = conformal.pred.split(
        mod.x,
        train.y,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.rdg$train,
        predict.fun = funs.rdg$predict,
        mad.train.fun = funs.res$train,
        mad.predict.fun = funs.res$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      #res.CP$LACP$RF.rf$out=out.cp
      res.CP$LACP$RDG.rf = create_tab_out.cp(out.cp,
                                             y.true = test.y,
                                             alpha = alpha,
                                             y_min = 0)
      
      if ("RF.sel" %in% models)
      {
        out.cp = conformal.pred.split(
          mod.x,
          train.y,
          x0 = as.matrix(mod.x.ts),
          alpha = alpha,
          rho = rho,
          train.fun = funs.rdg$train,
          predict.fun = funs.rdg$predict,
          mad.train.fun = funs.rf.sel$train,
          mad.predict.fun = funs.rf.sel$predict,
          seed = seed,
          verbose = T,
          split = i1
        )
        #res.CP$LACP$RF.rf$out=out.cp
        res.CP$LACP$RDG.rf.sel = create_tab_out.cp(out.cp,
                                                   y.true = test.y,
                                                   alpha = alpha,
                                                   y_min = 0)
      }
      
    }
    if ("SVM" %in% models)
    {
      if (best_par$SVM$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$SVM$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$SVM$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      mod.y.tr = tr.y.std
      out.cp$lo = out.cp$lo * y_sd + y_mu
      out.cp$up = out.cp$up * y_sd + y_mu
      out.cp = conformal.pred.split(
        mod.x,
        mod.y.tr,
        x0 = as.matrix(mod.x.ts),
        alpha = alpha,
        rho = rho,
        train.fun = funs.svm$train,
        predict.fun = funs.svm$predict,
        mad.train.fun = funs.res$train,
        mad.predict.fun = funs.res$predict,
        seed = seed,
        verbose = T,
        split = i1
      )
      #res.CP$LACP$RF.rf$out=out.cp
      out.cp$lo = out.cp$lo * y_sd + y_mu
      out.cp$up = out.cp$up * y_sd + y_mu
      res.CP$LACP$SVM.rf = create_tab_out.cp(out.cp,
                                             y.true = test.y,
                                             alpha = alpha,
                                             y_min = 0)
      
      if ("RF.sel" %in% models)
      {
        out.cp = conformal.pred.split(
          mod.x,
          train.y,
          x0 = as.matrix(mod.x.ts),
          alpha = alpha,
          rho = rho,
          train.fun = funs.svm$train,
          predict.fun = funs.svm$predict,
          mad.train.fun = funs.rf.sel$train,
          mad.predict.fun = funs.rf.sel$predict,
          seed = seed,
          verbose = T,
          split = i1
        )
        #res.CP$LACP$RF.rf$out=out.cp
        res.CP$LACP$SVM.rf.sel = create_tab_out.cp(out.cp,
                                                   y.true = test.y,
                                                   alpha = alpha,
                                                   y_min = 0)
      }
      
    }
    
    rm(funs.res)
    
    #QR----
    funs.rf.qr = quant.rf_ranger.funs(
      nthreads = N.threads,
      keep.inbag = F,
      importance = F,
      ntree = best_par$RF[["num.trees"]],
      replace = ifelse(best_par$RF[["sample.fraction"]] <
                         1, F, T),
      nodesize = best_par$RF[["min.node.size"]],
      varfrac = best_par$RF[["mtry"]] /
        p
    )
    
    out.cp = funs.rf.qr$train.fun(x = train.x,
                                  y = train.y,
                                  q = c(alpha / 2, 1 - alpha / 2))
    
    res.CP$QR = create_tab_out.cp(
      out.cp = list(
        lo = funs.rf.qr$predict.fun(out.cp, newx = test.x, q = alpha / 2),
        up = funs.rf.qr$predict.fun(out.cp, newx = test.x, q = 1 - alpha /
                                      2)
      ),
      y.true = test.y,
      y_min = 0,
      alpha = alpha
    )
    #CQR-------
    funs.rf.cqr = quant.rf_ranger.funs(
      nthreads = N.threads,
      keep.inbag = T,
      importance = F,
      ntree = best_par$RF[["num.trees"]],
      replace = ifelse(best_par$RF[["sample.fraction"]] <
                         1, F, T),
      nodesize = best_par$RF[["min.node.size"]],
      varfrac = best_par$RF[["mtry"]] /
        p
    )
    out.cp = conformal.quant.split(
      x = train.x,
      y = train.y,
      x0 = as.matrix(test.x),
      train.quant.fun = funs.rf.cqr$train.fun,
      predict.quant.fun = funs.rf.cqr$predict.fun,
      alpha = alpha,
      rho = rho,
      method = "all",
      tail.error.control = "both",
      score.type = "both",
      q.probs = c(alpha / 2, 1 - alpha / 2),
      tau_by = 0.005,
      max.n.tau = 50,
      exact_coverage = T,
      tune_quantile_CV = T,
      coverage_factor = coverage_factor,
      predict.quant.bootstrap.fun = funs.rf.cqr$predict.bootstrap.fun,
      split = i1,
      seed = seed,
      verbose = F
    )
    
    new_res = lapply(names(out.cp$Intervals), function(ty) {
      create_tab_out.cp(
        out.cp = out.cp$Intervals[[ty]],
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
    })
    names(new_res) = names(out.cp$Intervals)
    
    res.CP[["CQR_cv"]] = new_res
    
    rm(new_res)
    #DCP----
    out.cp = dcp.qr(
      Y0 = train.y[i1],
      X0 = train.x[i1, ],
      Y1 = train.y[-i1],
      X1 = train.x[-i1, ],
      Y.test = test.y,
      X.test = test.x,
      taus = seq(0.001, 0.999, length = 150),
      alpha.sig = alpha,
      score.type = "both"
    )
    
    
    new_res = lapply(names(out.cp), function(ty) {
      create_tab_out.cp(
        out.cp = out.cp[[ty]],
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
    })
    names(new_res) = names(out.cp)
    
    res.CP[["DCP"]] = new_res
    
    out.cp = dcp.dr(
      Y0 = train.y[i1],
      X0 = train.x[i1, ],
      Y1 = train.y[-i1],
      X1 = train.x[-i1, ],
      Y.test = test.y,
      X.test = test.x,
      taus = seq(0.001, 0.999, length = 150),
      alpha.sig = alpha,
      ys = quantile(unique(train.y), seq(0.001, 0.999, length = 1500))
    )
    res.CP[["DCP"]][["DCP.dr"]] = create_tab_out.cp(
      out.cp = out.cp,
      y.true = test.y,
      alpha = alpha,
      y_min = 0
    )
    
    #LCP----
    #rf all data
    out.cp = lcp_lcr(
      x = tr.std,
      y = train.y,
      x0 = ts.std,
      alpha = alpha,
      rho = rho,
      train.fun = funs.rf$train.fun,
      predict.fun = funs.rf$predict.fun,
      i1 = i1,
      seed = seed,
      verbose = T
    )
    res.CP[["LCP"]][["RF"]] = create_tab_out.cp(
      out.cp = list(lo = out.cp$lower, up = out.cp$upper),
      y.true = test.y,
      alpha = alpha,
      y_min = 0
    )
    #xgb all data
    if ("XGB" %in% models)
    {
      out.cp = lcp_lcr(
        x = tr.std,
        y = train.y,
        x0 = ts.std,
        alpha = alpha,
        rho = rho,
        train.fun = funs.xgb$train.fun,
        predict.fun = funs.xgb$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["XGB"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #knn
    if ("KNN" %in% models)
    {
      if (best_par$KNN$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$KNN$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$KNN$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = lcp_lcr(
        x = tr.std,
        y = train.y,
        x0 = ts.std,
        alpha = alpha,
        rho = rho,
        train.fun = funs.knn$train.fun,
        predict.fun = funs.knn$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["KNN"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #rdg
    if ("RDG" %in% models)
    {
      if (best_par$RDG$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RDG$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RDG$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = lcp_lcr(
        x = tr.std,
        y = train.y,
        x0 = ts.std,
        alpha = alpha,
        rho = rho,
        train.fun = funs.rdg$train.fun,
        predict.fun = funs.rdg$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["RDG"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #svm
    if ("SVM" %in% models)
    {
      if (best_par$SVM$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$SVM$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$SVM$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      mod.y.tr = tr.y.std
      out.cp$lo = out.cp$lo * y_sd + y_mu
      out.cp$up = out.cp$up * y_sd + y_mu
      out.cp = lcp_lcr(
        x = tr.std,
        y = mod.y.tr,
        x0 = ts.std,
        alpha = alpha,
        rho = rho,
        train.fun = funs.svm$train.fun,
        predict.fun = funs.svm$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      out.cp$lower = out.cp$lower * y_sd + y_mu
      out.cp$upper = out.cp$upper * y_sd + y_mu
      res.CP[["LCP"]][["SVM"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #RF.sel
    if ("RF.sel" %in% models)
    {
      if (best_par$RF.sel$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RF.sel$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RF.sel$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      p.mod = ncol(mod.x)
      
      out.cp = lcp_lcr(
        x = mod.x,
        y = train.y,
        x0 = mod.x.ts,
        alpha = alpha,
        rho = rho,
        train.fun = funs.rf.sel$train.fun,
        predict.fun = funs.rf.sel$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["RF.sel"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #xgb.sel
    if ("XGB.sel" %in% models)
    {
      funs.xgb.sel = xgb.funs(
        nrounds = best_par$XGB.sel[["nrounds"]],
        eta = best_par$XGB.sel[["eta"]],
        max_depth =  best_par$XGB.sel[["max_depth"]],
        min_child_weight = best_par$XGB.sel[["min_child_weight"]],
        varfrac = best_par$XGB.sel[["colsample_bytree"]],
        nthreads = N.threads
      )
      if (best_par$XGB.sel$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$XGB.sel$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$XGB.sel$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      out.cp = lcp_lcr(
        x = mod.x,
        y = train.y,
        x0 = mod.x.ts,
        alpha = alpha,
        rho = rho,
        train.fun = funs.xgb.sel$train.fun,
        predict.fun = funs.xgb.sel$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["XGB.sel"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    #LCQR rf
    out.cp = lcp_lcqr(
      x = tr.std,
      y = train.y,
      x0 = ts.std,
      alpha = alpha,
      rho = rho,
      train.quant.fun = funs.rf.cqr$train.fun,
      predict.quant.fun = funs.rf.cqr$predict.fun,
      i1 = i1,
      seed = seed,
      verbose = T
    )
    res.CP[["LCP"]][["LCQR_rf"]] = create_tab_out.cp(
      out.cp = list(lo = out.cp$lower, up = out.cp$upper),
      y.true = test.y,
      alpha = alpha,
      y_min = 0
    )
    
    if ("RF.sel" %in% models)
    {
      if (best_par$RF.sel$type == "PCA")
      {
        var.sel = var_tau.pca[[best_par$RF.sel$tau]]
        mod.x = tr.pca[, var.sel]
        mod.x.ts = ts.pca[, var.sel]
      } else {
        var.sel = var_tau.stability[[best_par$RF.sel$tau]]
        mod.x = tr.std[, var.sel]
        mod.x.ts = ts.std[, var.sel]
      }
      p.mod = ncol(mod.x)
      funs.rf.cqr.sel = quant.rf_ranger.funs(
        nthreads = N.threads,
        keep.inbag = T,
        importance = F,
        ntree = best_par$RF.sel[["num.trees"]],
        replace = ifelse(best_par$RF.sel[["sample.fraction"]] <
                           1, F, T),
        nodesize = best_par$RF.sel[["min.node.size"]],
        varfrac = best_par$RF.sel[["mtry"]] /
          p.mod
      )
      
      
      out.cp = lcp_lcqr(
        x = mod.x,
        y = train.y,
        x0 = mod.x.ts,
        alpha = alpha,
        rho = rho,
        train.quant.fun = funs.rf.cqr.sel$train.fun,
        predict.quant.fun = funs.rf.cqr.sel$predict.fun,
        i1 = i1,
        seed = seed,
        verbose = T
      )
      res.CP[["LCP"]][["LCQR_rf.sel"]] = create_tab_out.cp(
        out.cp = list(lo = out.cp$lower, up = out.cp$upper),
        y.true = test.y,
        alpha = alpha,
        y_min = 0
      )
      
    }
    
    #CREATE TAB_RES----
    f.calc.stat.res = function(tab.out.cp) {
      vect = select(tab.out.cp, coverage, width, w.score) |>
        mutate(coverage = coverage * 100) |>
        filter(is.finite(width)) |>
        colMeans(na.rm = T)
      I.inf = mean(is.infinite(tab.out.cp[["width"]])) * 100
      if (I.inf > 0)
        vect[["coverage"]] = vect[["coverage"]] + I.inf
      names(I.inf) = "I.inf"
      vect = append(vect, I.inf)
      return(vect)
    }
    
    tab_res = sapply(names(res.CP), function(ty) {
      if (!ty %in% c("benchmark", "QR")) {
        sapply(names(res.CP[[ty]]), function(j) {
          f.calc.stat.res(res.CP[[ty]][[j]])
        })
      } else {
        f.calc.stat.res(res.CP[[ty]])
      }
    })
    
    tab_res = do.call(cbind, lapply(names(tab_res), function(nm) {
      m <- tab_res[[nm]]
      m = as.matrix(m, nrow = 3)
      if (!nm %in% c("benchmark", "QR", "CQR", "CQR_cv", "DCP")) {
        colnames(m) <- paste(nm, colnames(m), sep = "_")
      } else if (nm %in% c("benchmark", "QR")) {
        colnames(m) <- paste(nm)
      }
      m
    })) |>
      t()
    #WRITE----
    name.file = paste0(name.dir, "_rep", rep, ".csv")
    tab_res = tab_res |>
      as_tibble() |>
      mutate(type = rownames(tab_res)) |>
      relocate(type) |>
      as.data.frame()
    rownames(tab_res) = NULL
    
    write.csv(tab_res,
              file =  file.path(paste0("result.superconductor/", name.dir), name.file),
              row.names = F)
    if (!dir.exists(paste0("result.superconductor/", name.dir, "/single_prediction"))) {
      dir.create(
        paste0(
          "result.superconductor/",
          name.dir,
          "/single_prediction"
        ),
        recursive = TRUE
      )
    }
    write.csv(
      mutate(res.CP[["benchmark"]], y = test.y, i.pos = i.test),
      file =  file.path(
        paste0("result.superconductor/", name.dir, "/single_prediction"),
        paste0("benchmark_", rep, ".csv")
      ),
      row.names = F
    )
    write.csv(
      mutate(res.CP[["CQR_cv"]][["CQR_optim"]], y = test.y, i.pos = i.test),
      file =  file.path(
        paste0("result.superconductor/", name.dir, "/single_prediction"),
        paste0("CQR_optim_", rep, ".csv")
      ),
      row.names = F
    )
  }
  
  tab.avg = NULL
  for (name.file in list.files(path = paste0("result.superconductor/", name.dir),
                               pattern = "\\.csv$"))
  {
    if (!grepl(x = name.file, "_AVG"))
    {
      tab.avg = rbind(tab.avg,
                      read_csv(
                        file = paste("result.superconductor", name.dir, name.file, sep = "/")
                      ) |>
                        mutate(rep = strsplit(name.file, "rep|\\.csv")[[1]][2] |>
                                 as.numeric()))
    }
  }
  
  write.csv(tab.avg,
            file =  file.path(
              paste0("result.superconductor/", name.dir),
              paste0(name.dir, "_AVG.csv")
            ),
            row.names = F)
}
