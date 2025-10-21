#########CQR_helper_functions#######

#Functions used within the “conformal.quant.split” function
#smoothed_cp-----
#if necessary and if “exact_coverage = T” applies smoothed CP to obtain exact
#coverage on the calibration set (score)
smoothed_cp = function(score,
                       significance = 1 - alpha,
                       n.test,
                       exact_coverage = T)
{
  k = ceiling(significance * length(score))
  score_sorted = sort(score)
  q.star = score_sorted[k]
  smaller_score = score_sorted[score_sorted < q.star]
  
  #if: no smoothed CP
  if (!sum(score <= q.star) > k |
      length(smaller_score) == 0 | !exact_coverage)
    return(cbind(rep(q.star, n.test)))
  
  #else: smoothed CP:
  n = length(score)
  n_less = sum(score < q.star)
  n_equal = sum(score == q.star)
  tau = (significance * n  - n_less) / n_equal #(n=n2+1 beacuse was added +Inf to scores)
  q.prev = max(smaller_score)
  
  #tie‑breaking
  u = runif(n.test)
  ifelse(u <= tau, q.star, q.prev)
  
  return(cbind(ifelse(u <= tau, q.star, q.prev)))
}
#f.pos-----
#Given the length of multiple intervals, select the position of the interval with
#the minimum length. If there are several intervals with the minimum length,
#select the most central one.
f.pos = function(xi)
{
  pos.int = which(xi == min(xi))
  #if there are multiple intervals with minimum width, the most “central” one is considered
  if (length(pos.int) > 1)
    pos.int = pos.int[which.min(abs(length(xi) / 2 - pos.int))]
  return(pos.int)
}
#small.interval-----
#For each row of X, calculate the relative “optimal” interval.
small.interval = function(X,
                          taus,
                          out = NULL,
                          predict.quant.fun = NULL,
                          q.taus = NULL)
{
  N = nrow(X)
  #number of taus couple
  n.taus = length(taus) / 2
  #predict
  if (is.null(q.taus))
    q.taus = matrix(predict.quant.fun(out, q = taus, X),
                    nrow = N,
                    ncol = 2 * n.taus)
  t.len = q.taus[, (n.taus + 1):(n.taus * 2)] - q.taus[, 1:n.taus]
  t.pos = apply(t.len, 1, f.pos)
  
  rows = 1:N
  coords_lo = cbind(rows, t.pos)
  coords_hi = cbind(rows, t.pos + n.taus)
  return(list(
    Int = cbind(q.taus[coords_lo], q.taus[coords_hi]),
    pos_col = cbind(t.pos, t.pos + n.taus)
  ))
}
#tune_CV_q-----
#is tuned with 5-CV CQR quantile levels, and/or the distance of the optimal-CQR
#quantile levels.
tune_CV_q = function(X = x[i1, ],
                     Y = y[i1],
                     score.type = score.type,
                     train.quant.fun = train.quant.fun,
                     predict.quant.fun = predict.quant.fun,
                     alpha = 0.1,
                     coverage_factor = 0.95,
                     tau_by = tau_by,
                     max.n.tau = 50)
{
  target_coverage = (1 - alpha) * coverage_factor
  range_less = 0.25
  range_add = 0.1
  vals_by = 0.01
  all_q = c()
  grid_len = seq(
    from = max(target_coverage - range_less, 0.01),
    to = min(target_coverage + range_add, .99),
    by = vals_by
  ) |>
    round(3) |>
    append(1 - alpha) |>
    unique() |>
    sort(decreasing = T)
  if ("basic" %in% score.type)
  {
    grid_q_low = round((1 - grid_len) / 2, 3)
    grid_q_hi = round(1 - grid_q_low, 3)
    grid_q = cbind(grid_q_low, grid_q_hi)
    n.q.couple = nrow(grid_q)
    
    res_q = matrix(0, ncol = 2, nrow = n.q.couple)
    colnames(res_q) = c("coverage", "width")
    
    all_q.basic = as.vector(grid_q) |>
      sort()
    all_q = append(all_q, all_q.basic)
  }
  
  if ("optimal" %in% score.type)
  {
    q_taus = lapply(grid_len, function(len_i) {
      taus_lo = seq(0, 1 - len_i, by = tau_by)
      if (length(taus_lo) > max.n.tau)
        taus_lo = seq(0, 1 - len_i, length.out = max.n.tau)
      taus_hi = taus_lo + len_i
      return(round(cbind(taus_lo, taus_hi), 3))
    })
    res_l = matrix(0, ncol = 2, nrow = length(grid_len))
    colnames(res_l) = c("coverage", "width")
    
    all_q.optimal = unlist(q_taus) |>
      round(3) |>
      unique() |>
      sort()
    
    all_q = append(all_q, all_q.optimal) |>
      unique() |>
      sort()
  }
  
  pos.q.basic = match(all_q.basic, all_q)
  pos.q.optim = lapply(q_taus, function(q_i) {
    match(as.vector(q_i), all_q)
  })
  
  #cv
  K = 5
  folds = fold_assignment(nrow(X), K)
  weights.folds = table(folds) / nrow(X)
  
  for (k.f in seq_len(K))
  {
    id_f = which(folds == k.f)
    mod = train.quant.fun(x = X[-id_f, ], y = Y[-id_f], q = all_q)
    x_test = X[id_f, ]
    y_test = Y[id_f]
    pred_all = predict.quant.fun(mod, q = all_q, x_test)
    w_fold = weights.folds[k.f]
    if ("basic" %in% score.type) {
      pred = pred_all[, pos.q.basic]
      for (i in seq_len(n.q.couple)) {
        q_lo = pred[, i]
        q_hi = pred[, n.q.couple * 2 + 1 - i]
        res_q[i, "width"] = res_q[i, "width"] + mean(q_hi - q_lo) * w_fold
        res_q[i, "coverage"] = res_q[i, "coverage"] + mean(y_test >= q_lo &
                                                             y_test <= q_hi) * w_fold
      }
    }
    if ("optimal" %in% score.type) {
      for (i in seq_len(length(grid_len))) {
        taus.i = as.vector(q_taus[[i]])
        pred = pred_all[, pos.q.optim[[i]]]
        tau_Int = small.interval(X = x_test,
                                 taus = taus.i,
                                 q.taus = pred)$Int
        res_l[i, "width"] = res_l[i, "width"] + mean(tau_Int[, 2] - tau_Int[, 1]) * w_fold
        res_l[i, "coverage"] = res_l[i, "coverage"] + mean(y_test >= tau_Int[, 1] &
                                                             y_test <= tau_Int[, 2]) * w_fold
      }
      
    }
  }
  res = list()
  if ("basic" %in% score.type) {
    pos.best = which(res_q[, "coverage"] >= target_coverage)
    pos.best = pos.best[which.min(res_q[pos.best, "width"])] |>
      tail(1) # if there are multiple intervals with min length
    if (length(pos.best) == 0)
      pos.best = 1
    res[["basic"]] = grid_q[pos.best, ]
    names(res[["basic"]]) = c("q_lo", "q_up")
  }
  if ("optimal" %in% score.type) {
    pos.best = which(res_l[, "coverage"] >= target_coverage)
    pos.best = pos.best[which.min(res_l[pos.best, "width"])] |>
      tail(1) # if there are multiple intervals with min length
    if (length(pos.best) == 0)
      pos.best = 1
    res[["optimal"]] = q_taus[[pos.best]]#if there are multiple intervals with min length
  }
  return(res)
  
}

#B.bootstrap_prediction-----
#bootstrap predictions are obtained to apply UACQR methods. This function will
#be used if no ad hoc function is specified based on the underlying model.
B.bootstrap_prediction = function(B = 100,
                                  x.train = x[i1, , drop = F],
                                  y.train = y[i1],
                                  q = q.boot,
                                  list_newx = list(x[i2, , drop = F], x0),
                                  train.quant.fun = train.quant.fun,
                                  predict.quant.fun = predict.quant.fun)
{
  n.q = length(q)
  n.train = nrow(x.train)
  n.new_dataset = length(list_newx)
  end_tab = sapply(list_newx, nrow) |>
    cumsum()
  start_tab = c(1, head(end_tab, -1) + 1)
  newX = do.call(rbind, list_newx)
  
  bootstrap = lapply(seq_len(n.new_dataset), function(j)
    array(
      NA_real_,
      dim = c(nrow(list_newx[[j]]), n.q, B),
      dimnames = list(rownames(list_newx[[j]]), paste0(q * 100, "%"), NULL)
    ))
  for (b in seq_len(B))
  {
    iboot = sample(seq_len(n.train), replace = T, size = n.train)
    mod_b = train.quant.fun(x.train[iboot, ], y.train[iboot], q = q)
    tab = predict.quant.fun(mod_b, q, newx = newX)
    for (j in seq_len(n.new_dataset))
      bootstrap[[j]][, , b] = tab[start_tab[j]:end_tab[j], ]
  }
  return(bootstrap)
}

