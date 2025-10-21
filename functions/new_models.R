#######NEW MODELS FUNCTIONS###########
#------------
#A) Functions of new regression models that can be used with the conformalInference library
#------------
#library(devtools)
#install_github(repo="ryantibs/conformal", subdir="conformalInference")
#library(conformalInference)
#xgb.funs------
xgb.funs = function(nrounds = 50,
                    eta = 0.3,
                    gamma = 0,
                    lambda_L2 = 1,
                    alpha_L1 = 0,
                    max_depth = 6,
                    min_child_weight = 1,
                    varfrac = 1,
                    obsfrac = 1,
                    nthreads = 1)
  #si possono togliere
{
  if (!require("xgboost", quietly = TRUE)) {
    stop("Package xgboost not installed (required here)!")
  }
  conformalInference:::check.pos.int(nrounds)
  conformalInference:::check.num.01(eta)
  conformalInference:::check.pos.num(gamma)
  conformalInference:::check.pos.num(lambda_L2)
  conformalInference:::check.pos.num(alpha_L1)
  conformalInference:::check.pos.int(max_depth)
  conformalInference:::check.pos.int(min_child_weight)
  conformalInference:::check.pos.int(nthreads)
  conformalInference:::check.num.01(varfrac)
  conformalInference:::check.num.01(obsfrac)
  
  train.fun = function(x, y, out = NULL) {
    dmat = xgb.DMatrix(data = x, label = y)
    return(
      xgboost(
        data = dmat,
        nrounds = nrounds,
        eta = eta,
        gamma = gamma,
        lambda = lambda_L2,
        alpha = alpha_L1,
        max_depth = max_depth,
        min_child_weight = min_child_weight,
        colsample_bytree = varfrac,
        subsample = obsfrac,
        nthread = nthreads
      )
    )
  }
  predict.fun = function(out, newx) {
    return(predict(out, xgb.DMatrix(data = newx)))
  }
  active.fun = function(out) {
    return(list(which(
      out$feature_names %in%
        xgb.importance(model = out)$Feature
    )))
  }
  return(list(
    train.fun = train.fun,
    predict.fun = predict.fun,
    active.fun = active.fun
  ))
}
#knn.funs-----
kknn.funs = function(k = 5,
                     distance = 2,
                     kernel = "epanechnikov",
                     scale = F)
{
  if (!require("kknn", quietly = TRUE)) {
    stop("Package kknn not installed (required here)!")
  }
  poss.kernel = c(
    "rectangular",
    "triangular",
    "epanechnikov",
    "biweight",
    "triweight",
    "cos",
    "inv",
    "gaussian",
    "rank",
    "optimal"
  )
  conformalInference:::check.pos.int(k)
  conformalInference:::check.pos.num(distance)
  conformalInference:::check.bool(scale)
  if (!kernel %in% poss.kernel)
    stop(paste(
      deparse(substitute(kernel)),
      "must be one of:\n",
      paste(poss.kernel, collapse = ", ")
    ))
  
  train.fun = function(x, y, out = NULL) {
    return(out = data.frame(y = y, x))
  }
  predict.fun = function(out, newx) {
    df.new = data.frame(newx)
    if (is.null(colnames(newx)))
      colnames(df.new) = setdiff(colnames(out), "y")
    mod = kknn(
      y ~ .,
      train = out,
      test = df.new,
      k = k,
      distance = distance,
      kernel = kernel,
      scale = scale
    )
    return(predict(mod))
  }
  active.fun = function(out) {
    return(NULL)
  }
  return(list(
    train.fun = train.fun,
    predict.fun = predict.fun,
    active.fun = active.fun
  ))
}
#svm.funs-----
svm.funs = function(cost = 1,
                    epsilon = 0.1,
                    gamma = 1,
                    degree = 3,
                    coef0 = 0,
                    kernel = "linear",
                    scale = F)
{
  param.set_list$SVM$par.set
  if (!require("e1071", quietly = TRUE)) {
    stop("Package e1071 not installed (required here)!")
  }
  poss.kernel = c("radial", "polynomial", "linear", "sigmoid")
  conformalInference:::check.pos.int(degree)
  conformalInference:::check.num(cost)
  conformalInference:::check.num(epsilon)
  conformalInference:::check.num(gamma)
  conformalInference:::check.num(coef0)
  conformalInference:::check.bool(scale)
  if (!kernel %in% poss.kernel)
    stop(paste(
      deparse(substitute(kernel)),
      "must be one of:\n",
      paste(poss.kernel, collapse = ", ")
    ))
  
  train.fun = function(x, y, out = NULL) {
    return(
      svm(
        x = x,
        y = y,
        kernel  = kernel,
        cost    = cost,
        epsilon = epsilon,
        gamma   = gamma,
        degree  = degree,
        coef0   = coef0,
        scale   = scale
      )
    )
  }
  
  predict.fun = function(out, newx) {
    return(predict(out, newdata = newx))
  }
  active.fun = function(out) {
    return(NULL)
  }
  return(list(
    train.fun = train.fun,
    predict.fun = predict.fun,
    active.fun = active.fun
  ))
}
#rf_ranger.funs------
rf_ranger.funs = function (ntree = 500, varfrac = 0.333, replace = TRUE,
                           obsfrac = ifelse(replace, 1, 0.632), nodesize = 5,
                           nthreads = 0, keep.inbag = F, importance = T) 
{
  if (!require("ranger", quietly = TRUE)) {
    stop("Package ranger not installed (required here)!")
  }
  conformalInference:::check.pos.num(nthreads)
  conformalInference:::check.int(nthreads)
  conformalInference:::check.pos.int(ntree)
  conformalInference:::check.num.01(varfrac)
  conformalInference:::check.bool(replace)
  conformalInference:::check.num.01(obsfrac)
  conformalInference:::check.pos.int(nodesize)
  conformalInference:::check.bool(keep.inbag)
  conformalInference:::check.bool(importance)
  train.fun = function(x, y, out = NULL) {
    p = ncol(x)
    return(ranger::ranger(x = x ,y=y, num.trees = ntree, mtry = floor(varfrac*p),
                          min.node.size = nodesize, replace=replace, sample.fraction = obsfrac,
                          num.threads = nthreads, keep.inbag = keep.inbag,
                          importance = ifelse(importance, "impurity", "none") ))
  }
  predict.fun = function(out, newx) {
    var.names = out$forest$independent.variable.names
    newx = as.data.frame(newx)
    colnames(newx) = var.names
    return(predict(out, newx)$predictions)
  }
  
  active.fun = function(out) {
    return(list(which(out$variable.importance>0)))
  }
  return(list(train.fun = train.fun, predict.fun = predict.fun, 
              active.fun = active.fun))
}
#--------------
#B) Functions of quantile regression models that can be used with conformal.quant.split function
#-------------
#quant.rf_ranger.funs------
quant.rf_ranger.funs = function(ntree = 500,
                                varfrac = 0.333,
                                replace = TRUE,
                                obsfrac = ifelse(replace, 1, 0.632),
                                nodesize = 5,
                                nthreads = 0,
                                keep.inbag = T,
                                importance = T)
{
  if (!require("ranger", quietly = TRUE)) {
    stop("Package ranger not installed (required here)!")
  }
  # Check arguments----
  #if (is.null(q) || !length(q)==2 || !is.numeric(q) || any(q<0) || any(q>1))
  #  stop(paste(deparse(substitute(q)), "must be a vector of two numbers between 0 and 1"))
  #so nthreads must be a integer in [0,+Inf)
  conformalInference:::check.pos.num(nthreads)
  conformalInference:::check.int(nthreads)
  conformalInference:::check.pos.int(ntree)
  conformalInference:::check.num.01(varfrac)
  conformalInference:::check.bool(replace)
  conformalInference:::check.num.01(obsfrac)
  conformalInference:::check.pos.int(nodesize)
  conformalInference:::check.bool(keep.inbag)
  conformalInference:::check.bool(importance)
  #train----
  train.fun = function(x, y, q = NULL, out = NULL) {
    p = ncol(x)
    return(
      ranger(
        x = x ,
        y = y,
        quantreg = T,
        num.trees = ntree,
        mtry = floor(varfrac * p),
        min.node.size = nodesize,
        replace = replace,
        sample.fraction = obsfrac,
        num.threads = nthreads,
        keep.inbag = keep.inbag,
        importance = ifelse(importance, "impurity", "none")
      )
    )
  }
  #predict----
  predict.fun = function(out, q, newx) {
    return(
      predict(
        out,
        newx,
        type = "quantiles",
        quantiles = q,
        num.threads = nthreads
      )$predictions
    )
  }
  #bagging shortcut----
  predict.bagging.fun = function(out,
                                 x.train,
                                 y.train,
                                 q,
                                 list_newx = NULL,
                                 newx = NULL) {
    if (is.null(newx)) {
      n.new_dataset = length(list_newx)
      end_tab = sapply(list_newx, nrow) |>
        cumsum()
      start_tab = c(1, head(end_tab, -1) + 1)
      newx = do.call(rbind, list_newx)
    } else{
      list_newx = NULL
      n.new_dataset = 1
    }
    
    n.newx = nrow(newx)
    n.train = nrow(x.train)
    n.tree = out$num.trees
    n.q = length(q)
    
    in.bag = out$inbag.counts
    # For each tree, the leaves belonging to each observation in the training set are identified.
    nodes_train = predict(out, data = x.train, type = "terminalNodes")$predictions
    # For each tree, the leaves belonging to each new observation are identified.
    nodes_new = predict(out, data = newx, type = "terminalNodes")$predictions
    # quantiles estimated for each new observation for each tree boostrap
    if (n.new_dataset == 1) {
      quant_by.tree = array(
        NA_real_,
        dim = c(n.newx, n.q, n.tree),
        dimnames = list(rownames(newx), paste0(q * 100, "%"), NULL)
      )
    } else {
      quant_by.tree = lapply(seq_len(n.new_dataset), function(j)
        array(
          NA_real_,
          dim = c(nrow(list_newx[[j]]), n.q, n.tree),
          dimnames = list(rownames(list_newx[[j]]), paste0(q * 100, "%"), NULL)
        ))
    }
    for (b in seq_len(n.tree))
    {
      # The bootstrap sample corresponding to the tree is obtained
      pos.b = rep(seq_len(n.train), times = in.bag[[b]])
      # The training Ys are grouped by leaf of belonging of the tree
      y_by_leaf = split(y.train[pos.b], nodes_train[pos.b, b])
      # Only the leaves used for predictions are considered
      y_by_leaf = y_by_leaf[as.numeric(names(y_by_leaf)) %in% nodes_new[, b]]
      #For each leaf of interest, the quantiles are calculated.
      tab = vapply(
        y_by_leaf,
        FUN = quantile,
        probs = q,
        FUN.VALUE = numeric(n.q),
        na.rm = TRUE
      ) |>
        t() |>
        rbind("-1" = c(NA, NA))
      leaf_id = as.character(nodes_new[, b])
      # If a new observation ends up in a leaf without observations from the training set, you get NA
      leaf_id[which(!leaf_id %in% rownames(tab))] = "-1"
      if (n.new_dataset == 1L) {
        quant_by.tree[, , b] = tab[leaf_id, ]
      } else {
        tab = tab[leaf_id, ]
        for (j in seq_len(n.new_dataset))
          quant_by.tree[[j]][, , b] = tab[start_tab[j]:end_tab[j], ]
      }
    }
    
    return(quant_by.tree)
  }
  #active----
  active.fun = function(out) {
    return(which(importance(out) > 0))
  }
  return(
    list(
      train.fun = train.fun,
      predict.fun = predict.fun,
      predict.bootstrap.fun = predict.bagging.fun,
      active.fun = active.fun
    )
  )
}
#quant.rf.funs-----
quant.rf.funs = function(ntree = 500,
                         varfrac = 0.333,
                         replace = TRUE,
                         obsfrac = ifelse(replace, 1, 0.632),
                         nodesize = 5,
                         nthreads = 1,
                         keep.inbag = T)
{
  if (!require("quantregForest", quietly = TRUE)) {
    stop("Package quantregForest not installed (required here)!")
  }
  # Check arguments-----
  #if (is.null(q) || !length(q)==2 || !is.numeric(q) || any(q<0) || any(q>1))
  #  stop(paste(deparse(substitute(q)), "must be a vector of two numbers between 0 and 1"))
  conformalInference:::check.pos.int(nthreads)
  conformalInference:::check.pos.int(ntree)
  conformalInference:::check.num.01(varfrac)
  conformalInference:::check.bool(replace)
  conformalInference:::check.num.01(obsfrac)
  conformalInference:::check.pos.int(nodesize)
  conformalInference:::check.bool(keep.inbag)
  #train----
  train.fun = function(x, y, q = NULL, out = NULL) {
    n = nrow(x)
    p = ncol(x)
    return(
      quantregForest(
        x = x,
        y = y,
        nthreads = nthreads,
        ntree = ntree,
        mtry = floor(varfrac * p),
        replace = replace,
        sampsize = floor(obsfrac * n),
        nodesize = nodesize,
        keep.inbag = keep.inbag
      )
    )
  }
  #prediction----
  predict.fun = function(out, q, newx) {
    return(predict(out, newx, what = q))
  }
  #bagging shortcut----
  predict.bagging.fun = function(out, x.train, y.train, q, newx) {
    # Transform the model into a randomForest object
    out.rf = out
    class(out.rf) = "randomForest"
    n.newx = nrow(newx)
    n.train = nrow(x.train)
    n.tree = out.rf$ntree
    n.q = length(q)
    
    in.bag = out$inbag
    # For each tree, the leaves belonging to each observation in the training set are identified.
    nodes_train = attr(predict(out.rf, x.train, nodes = TRUE), "nodes")
    # For each tree, the leaves belonging to each new observation are identified.
    nodes_new = attr(predict(out.rf, newx , nodes = TRUE), "nodes")
    # quantiles estimated for each new observation for each tree boostrap
    quant_by.tree = array(
      NA_real_,
      dim = c(n.newx, n.q, n.tree),
      dimnames = list(rownames(newx), paste0(q * 100, "%"), NULL)
    )
    for (b in seq_len(n.tree))
    {
      # the bootstrap sample corresponding to the tree is obtained
      pos.b = rep(seq_len(n.train), times = in.bag[[b]])
      # The training Ys are grouped by leaf of belonging of the tree
      y_by_leaf = split(y.train[pos.b], nodes_train[pos.b, b])
      # only the leaves used for predictions are considered
      y_by_leaf = y_by_leaf[as.numeric(names(y_by_leaf)) %in% nodes_new[, b]]
      #For each leaf of interest, the quantiles are calculated.
      tab = vapply(
        y_by_leaf,
        FUN = quantile,
        probs = q,
        FUN.VALUE = numeric(n.q),
        na.rm = TRUE
      ) |>
        t() |>
        rbind("-1" = c(NA, NA))
      leaf_id = as.character(nodes_new[, b])
      # if a new observation ends up in a leaf without observations from the training set, you get NA
      leaf_id[which(!leaf_id %in% rownames(tab))] = "-1"
      quant_by.tree[, , b] = tab[leaf_id, ]
    }
    return(quant_by.tree)
  }
  active.fun = function(out) {
    return(list(which(out$importance > 0)))
  }
  return(
    list(
      train.fun = train.fun,
      predict.fun = predict.fun,
      predict.bootstrap.fun = predict.bagging.fun,
      active.fun = active.fun
    )
  )
}
