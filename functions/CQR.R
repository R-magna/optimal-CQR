########CQR########
#Function to apply CQR split with different variants.
#The function was written maintaining a structure similar to that of the
#“conformal.pred.split” function (for Slit CP) of "conformalInference" package.
#https://github.com/ryantibs/conformal/tree/master/conformalInference

#possible CQR variations:
# - method = ("CQR", "CQR-m", "CQR-r", "UACQR-s", "UACQR-p", "all")
# - tail.error.control = ("single", "separate", "both")
# - score.type = ("basic", "optimal", "both")

# - tune_quantile_CV = (T, F) (if is T, the tuning is based on coverage_factor)
# - exact_coverage = (T, F) (if is T, and if it is necessary smoothed CP is applied)

#conformal.quant.split-----

conformal.quant.split = function(x,
                                 y,
                                 x0,
                                 train.quant.fun,
                                 predict.quant.fun,
                                 alpha = 0.1,
                                 rho = 0.5,
                                 method = "CQR",
                                 tail.error.control = "single",
                                 score.type = "basic",
                                 q.probs = c(alpha / 2, 1 - alpha / 2),
                                 tau_by = ifelse(score.type == "basic", NULL, 0.005),
                                 max.n.tau = ifelse(score.type == "basic", NULL, 50),
                                 exact_coverage = T,
                                 tune_quantile_CV = T,
                                 coverage_factor = ifelse(tune_quantile_CV, 0.85, NULL),
                                 predict.quant.bootstrap.fun = NULL,
                                 split = NULL,
                                 seed = NULL,
                                 verbose = FALSE)
{
  require(data.table)
  # Set up data----
  x = as.data.table(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = as.data.table(x0, ncol = p)
  nx0 = nrow(x0)
  
  possible.method = c("CQR", "CQR-m", "CQR-r", "UACQR-s", "UACQR-p", "all")
  possible.tail = c("single", "separate", "both")
  possible.score = c("basic", "optimal", "both")
  # Check input arguments ----
  #conformalInference:::check.args(x=x,y=y,x0=x0,alpha=alpha,train.fun=train.quant.fun,
  #                                predict.fun=predict.quant.fun)
  # Check rho
  if (is.null(rho) ||
      length(rho) != 1 || !is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("rho must be a number in between 0 and 1")
  
  # Check tau_by
  if (!is.null(tau_by) &
      (length(tau_by) != 1 ||
       !is.numeric(tau_by) || tau_by <= 0 || tau_by >= 1))
    stop("tau_by must be a number in between 0 and 1")
  
  # Check q.probs
  if (is.null(q.probs) ||
      !length(unique(q.probs)) == 2 ||
      !is.numeric(q.probs) || any(q.probs < 0) || any(q.probs > 1))
    stop(paste(
      deparse(substitute(q.probs)),
      "must be a vector of two different numbers between 0 and 1"
    ))
  
  #Check tail.error.control
  if (is.null(tail.error.control) ||
      !is.character(tail.error.control) ||
      !all(tolower(tail.error.control) %in% tolower(possible.tail)))
    stop(paste0(
      deparse(substitute(tail.error.control)),
      " must be one or more of the following:\n \"",
      paste(possible.tail, collapse = "\", \""),
      "\"."
    ))
  
  #Check method
  if (is.null(method) ||
      !is.character(method) ||
      !all(tolower(method) %in% tolower(possible.method)))
    stop(paste0(
      deparse(substitute(method)),
      " must be one or more of the following:\n \"",
      paste(possible.method, collapse = "\", \""),
      "\"."
    ))
  #Check score.type
  if (is.null(score.type) ||
      !is.character(score.type) ||
      !all(tolower(score.type) %in% tolower(possible.score)))
    stop(paste0(
      deparse(substitute(score.type)),
      " must be one or more of the following:\n \"",
      paste(possible.score, collapse = "\", \""),
      "\"."
    ))
  
  #Fix correct values of method
  if (any(tolower(method) == "all")) {
    method = possible.method[possible.method != "all"]
  } else
    method = possible.method[tolower(possible.method) %in% tolower(method)]
  
  #Fix correct values of tail.error.control
  if (any(tolower(tail.error.control) == "both")) {
    tail.error.control = possible.tail[possible.tail != "both"]
  } else
    tail.error.control = possible.tail[tolower(possible.tail) %in% tolower(tail.error.control)]
  
  #Fix correct values of score.type
  if (any(tolower(score.type) == "both")) {
    score.type = possible.score[possible.score != "both"]
  } else
    score.type = possible.score[tolower(possible.score) %in% tolower(score.type)]
  
  # Users may pass in a string for the verbose argument
  if (verbose == TRUE)
    txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # SPLIT If the user passed indices for the split, use them----
  if (!is.null(split)) {
    i1 = split
    # Otherwise make a random split
  } else {
    if (!is.null(seed))
      set.seed(seed)
    i1 = sample(1:n, floor(n * rho))
  }
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n", txt, n1, n2))
    cat(sprintf("%sTraining on first part ...\n", txt))
  }
  
  # Train on first part----
  q.probs = sort(q.probs)
  q.train = c()
  if (tune_quantile_CV)
  {
    CV_tuned = tune_CV_q(
      X = x[i1, ],
      Y = y[i1],
      score.type = score.type,
      train.quant.fun = train.quant.fun,
      predict.quant.fun = predict.quant.fun,
      alpha = alpha,
      coverage_factor = coverage_factor,
      tau_by = tau_by,
      max.n.tau = max.n.tau
    )
    if ("basic" %in% score.type)
      q.probs = CV_tuned$basic
    if ("optimal" %in% score.type)
    {
      taus_lo = CV_tuned[["optimal"]][, 1]
      n.taus = length(taus_lo)
      taus = c(taus_lo, CV_tuned[["optimal"]][, 2])
      q.train = unique(c(q.train, taus)) |>
        sort()
    }
  } else{
    if ("optimal" %in% score.type)
    {
      taus_lo = seq(0, alpha, by = tau_by)
      if (length(taus_lo) > max.n.tau)
        taus_lo = seq(0, alpha, length.out = max.n.tau)
      n.taus = length(taus_lo)
      taus = c(taus_lo, 1 - (alpha - taus_lo))
      q.train = unique(c(q.train, taus)) |>
        sort()
    }
  }
  q.train = unique(c(q.train, q.probs)) |>
    sort()
  
  if ("CQR-m" %in% method)
  {
    q.train = unique(c(q.train, 0.5)) |>
      sort()
  }
  
  #set.seed(123)
  out = train.quant.fun(x[i1, , drop = F], y[i1], q = q.train)
  fit = matrix(predict.quant.fun(out, q = q.probs, x), nrow = n)
  if (verbose) {
    cat(sprintf(
      "%sComputing residuals and quantiles on second part ...\n",
      txt
    ))
  }
  
  # Get score and quantiles, and prediction on second----
  eps = 1e-6
  if ("basic" %in% score.type)
  {
    pred.calib = matrix(predict.quant.fun(out, x[i2, , drop = F], q = q.probs),
                        nrow = n2,
                        ncol = 2)
    pred = matrix(predict.quant.fun(out, q = q.probs, x0), nrow = nx0)
    res = cbind(pred.calib[, 1] - y[i2], y[i2] - pred.calib[, 2])
  }
  
  if ("optimal" %in% score.type)
  {
    pos_col_calib.optim = small.interval(x[i2, , drop = F],
                                         taus = taus,
                                         out = out,
                                         predict.quant.fun = predict.quant.fun)
    pos_col_pred.optim = small.interval(x0,
                                        taus = taus,
                                        out = out,
                                        predict.quant.fun = predict.quant.fun)
    
    pred.calib.optim = pos_col_calib.optim[["Int"]]
    pred.optim = pos_col_pred.optim[["Int"]]
    res.optim = cbind(pred.calib.optim[, 1] - y[i2], y[i2] - pred.calib.optim[, 2])
  }
  
  if ("CQR-m" %in% method) {
    med.calib = predict.quant.fun(out, x[i2, , drop = F], q = 0.5)
    pred.med = predict.quant.fun(out, x0, q = 0.5)
  }
  
  if (any(c("UACQR-s", "UACQR-p") %in% method))
  {
    q.boot = c()
    if ("optimal" %in% score.type) {
      q.boot = taus
      pos_col_calib.optim = pos_col_calib.optim[["pos_col"]]
      pos_col_pred.optim = pos_col_pred.optim[["pos_col"]]
    }
    if ("basic" %in% score.type) {
      q.boot = unique(c(q.boot, q.probs)) |>
        sort()
      pos_col_q.basic = which(q.boot %in% q.probs)
      pos_col_calib = matrix(rep(pos_col_q.basic, each = n2),
                             ncol = 2,
                             byrow = FALSE)
      pos_col_pred = matrix(rep(pos_col_q.basic, each = nx0),
                            ncol = 2,
                            byrow = FALSE)
    }
    #correction of position in case q.probs are not in taus
    if (all(c("basic", "optimal") %in% score.type) &
        any(!q.probs %in% taus))
    {
      pos_col_remove.optim = which(q.boot == q.probs[which(!q.probs %in% taus)])
      for (j in pos_col_remove.optim)
      {
        i_add = which(pos_col_calib.optim >= j)
        pos_col_calib.optim[i_add] = pos_col_calib.optim[i_add] + 1
        
        i_add = which(pos_col_pred.optim >= j)
        pos_col_pred.optim[i_add] = pos_col_pred.optim[i_add] + 1
      }
    }
    gc()
    if (!is.null(predict.quant.bootstrap.fun))
    {
      bootstrap = predict.quant.bootstrap.fun(
        out,
        x.train = x[i1, , drop = F],
        y.train = y[i1],
        q = q.boot,
        list_newx = list(x[i2, , drop = F], x0)
      )
    } else {
      #bootstrap
      B = 100
      bootstrap = B.bootstrap_prediction(
        B = B,
        x.train = x[i1, , drop = F],
        y.train = y[i1],
        q = q.boot,
        list_newx = list(x[i2, , drop =
                             F], x0),
        train.quant.fun = train.quant.fun,
        predict.quant.fun = predict.quant.fun
      )
    }
    
    pred.calib.boot = bootstrap[[1]]
    pred.boot = bootstrap[[2]]
    rm(bootstrap)
    gc()
    
    if ("UACQR-s" %in% method) {
      sd.boot_f = function(hat_boot, pos_col)
      {
        sapply(1:nrow(hat_boot), function(i)
          apply(hat_boot[i, pos_col[i, ], ], MARGIN = 1, sd, na.rm = T)) |>
          t()
      }
      if ("basic" %in% score.type) {
        sd.boot.calib = sd.boot_f(hat_boot = pred.calib.boot, pos_col = pos_col_calib)
        sd.boot.pred = sd.boot_f(hat_boot = pred.boot, pos_col = pos_col_pred)
      }
      if ("optimal" %in% score.type) {
        sd.boot.calib.optim = sd.boot_f(hat_boot = pred.calib.boot, pos_col = pos_col_calib.optim)
        sd.boot.pred.optim = sd.boot_f(hat_boot = pred.boot, pos_col = pos_col_pred.optim)
      }
    }
    if ("UACQR-p" %in% method) {
      p.boot_f = function(hat_boot, pos_col) {
        lapply(1:nrow(hat_boot), function(i) {
          apply(cbind(hat_boot[i, pos_col[i, ], ], -Inf, Inf), MARGIN = 1, function(col.j) {
            col.j[is.na(col.j)] = median(col.j[is.finite(col.j)], na.rm = T)
            sort(col.j)
          })
        })
      }
      if ("basic" %in% score.type) {
        pred.calib.boot.sorted = p.boot_f(hat_boot = pred.calib.boot, pos_col = pos_col_calib)
        pred.boot.sorted = p.boot_f(hat_boot = pred.boot, pos_col = pos_col_pred)
        score.boot_p = sapply(seq_len(n2), function(i) {
          q.b_i = pred.calib.boot.sorted[[i]]
          y_i = y[i2][i]
          c(sum(q.b_i[, 1] > y_i), sum(q.b_i[, 2] < y_i))
        }) |>
          t()
      }
      if ("optimal" %in% score.type) {
        pred.calib.boot.sorted.optim = p.boot_f(hat_boot = pred.calib.boot, pos_col = pos_col_calib.optim)
        pred.boot.sorted.optim = p.boot_f(hat_boot = pred.boot, pos_col = pos_col_pred.optim)
        score.boot_p.optim = sapply(seq_len(n2), function(i) {
          q.b_i = pred.calib.boot.sorted.optim[[i]]
          y_i = y[i2][i]
          c(sum(q.b_i[, 1] > y_i), sum(q.b_i[, 2] < y_i))
        }) |>
          t()
      }
    }
    rm(pred.calib.boot, pred.boot)
    gc()
  }
  #intervals----
  calculate.intervals = function(technique,
                                 tail,
                                 score.variety,
                                 exact_coverage) {
    #type of residues and predictions to be used
    switch(score.variety, "basic" = {
      R = res
      Hat.calib = pred.calib
      Hat.test = pred
      if (technique == "UACQR-s") {
        SD.boot_calib = sd.boot.calib
        SD.boot_test = sd.boot.pred
      } else if (technique == "UACQR-p") {
        P.score = score.boot_p
        P.boot.test_sorted = pred.boot.sorted
      }
    }, "optimal" = {
      R = res.optim
      Hat.calib = pred.calib.optim
      Hat.test = pred.optim
      if (technique == "UACQR-s") {
        SD.boot_calib = sd.boot.calib.optim
        SD.boot_test = sd.boot.pred.optim
      } else if (technique == "UACQR-p") {
        P.score = score.boot_p.optim
        P.boot.test_sorted = pred.boot.sorted.optim
      }
    })
    #scores
    switch(
      technique,
      "CQR" = {
        score = R
      },
      "CQR-m" = {
        score = cbind(R[, 1] / (med.calib - Hat.calib[, 1] + eps),
                      R[, 2] / (Hat.calib[, 2] - med.calib + eps))
      },
      "CQR-r" = {
        diff.q = Hat.calib[, 2] - Hat.calib[, 1]
        score = apply(R, MARGIN = 2, function(j)
          j / (diff.q + eps))
      },
      "UACQR-s" = {
        score = cbind(R[, 1] / (SD.boot_calib[, 1] + eps), R[, 2] / (SD.boot_calib[, 2] + eps))
      },
      "UACQR-p" = {
        score = P.score
      }
    )
    if (technique != "UACQR-p")
      score = rbind(score, c(Inf, Inf))
    #quantile of interest
    switch(tail, "single" = {
      score =  apply(score, MARGIN = 1, max)
      
      q = smoothed_cp(
        score = score,
        significance = 1 - alpha,
        n.test = nrow(Hat.test),
        exact_coverage = exact_coverage
      )
      q = cbind(q, q)
      #q = apply(score, MARGIN = 1, max) |>
      #  quantile(probs = 1 - alpha) |>
      #  rep(2)
    }, "separate" = {
      q = cbind(
        smoothed_cp(
          score = score[, 1],
          significance =  1 - alpha / 2,
          n.test = nrow(Hat.test),
          exact_coverage = exact_coverage
        ),
        smoothed_cp(
          score = score[, 2],
          significance =  1 - alpha / 2,
          n.test = nrow(Hat.test),
          exact_coverage = exact_coverage
        )
      )
      
      #q = c(quantile(score[,1], probs = 1 - q.probs[1]),
      #      quantile(score[,2], probs = q.probs[2]))
    })
    #interval
    switch(
      technique,
      "CQR" = {
        I.lo = Hat.test[, 1] - q[, 1]
        I.up = Hat.test[, 2] + q[, 2]
      },
      "CQR-m" = {
        I.lo = Hat.test[, 1] - q[, 1] * (pred.med - Hat.test[, 1] + eps)
        I.up = Hat.test[, 2] + q[, 2] * (Hat.test[, 2] - pred.med + eps)
      },
      "CQR-r" = {
        I.lo = Hat.test[, 1] - q[, 1] * (Hat.test[, 2] - Hat.test[, 1] + eps)
        I.up = Hat.test[, 2] + q[, 2] * (Hat.test[, 2] - Hat.test[, 1] + eps)
      },
      "UACQR-s" = {
        I.lo = Hat.test[, 1] - q[, 1] * (SD.boot_test[, 1] + eps)
        I.up = Hat.test[, 2] + q[, 2] * (SD.boot_test[, 2] + eps)
      },
      "UACQR-p" = {
        #q = round(q) #PROVA AL POSTO DI round(q), ceiling(q) floor(q)
        B = nrow(P.boot.test_sorted[[1]]) - 2 #since +Inf and -Inf have been added
        Hat_p = sapply(seq_len(nx0), function(i) {
          q.b_i = P.boot.test_sorted[[i]]
          c(q.b_i[B + 2 - q[i, 1], 1], q.b_i[q[i, 2] + 1, 2])
        }) |>
          t()
        I.lo = Hat_p[, 1]
        I.up = Hat_p[, 2]
      }
    )
    
    name.pos = paste0(
      gsub(
        x = technique,
        pattern = "-",
        replacement = "."
      ),
      ifelse(tail == "single", "", "_sep"),
      #paste0("_",tail)),
      ifelse(score.variety == "basic", "", "_optim")
    ) #paste0("_",score.variety)))
    
    Int = list()
    Int[[name.pos]] = list(lo = I.lo, up = I.up)
    return(Int)
  }
  
  #CQR combinations to be calculated---
  tab = expand.grid(
    method,
    tail.error.control,
    score.type,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    as.matrix()
  
  if (nrow(tab) > 1)
    tab = tab[order(tab[, 1], -xtfrm(tab[, 2]), tab[, 3]), ]
  #interval calculation
  Intervals = apply(tab, MARGIN = 1, function(X)
    calculate.intervals(
      technique = X[1],
      tail = X[2],
      score.variety = X[3],
      exact_coverage = exact_coverage
    )) |>
    unlist(recursive = FALSE)
  
  if (length(Intervals) == 1)
  {
    Intervals = unlist(Intervals, recursive = FALSE)
    names(Intervals) = c("lo", "up")
  }
  return(list(
    Intervals = Intervals,
    fit = fit,
    split = i1
  ))
}
