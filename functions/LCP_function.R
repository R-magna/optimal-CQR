#####LCP##########
#Application of Localized Conformal Prediction.

#The functions were obtained by adapting the relative R functions of:
#Guan, Leying. "Localized conformal prediction: A generalized inference framework for conformal prediction."
#Biometrika 110.1 (2023): 33-50.
#https://github.com/LeyingGuan/LCP


# install.packages("devtools")
# devtools::install_github("LeyingGuan/LCP")
library(LCP)

#pairwise_dist-----
pairwise_dist = function(A, B = NULL) {
  A = as.matrix(A)
  if (is.null(B)){
    B = A } else{ 
    B = as.matrix(B) }
  # (||a||^2 + ||b||^2 - 2 a·b)^{1/2}, clamp numerica a >=0
  An = rowSums(A * A)
  Bn = rowSums(B * B)
  D2 = outer(An, Bn, "+") - 2 * tcrossprod(A, B)
  sqrt(pmax(D2, 0))
}

#lcp_lcr (LCP regression score):----------
lcp_lcr = function(x,
                   y,
                   x0,
                   train.fun,
                   predict.fun,
                   alpha = 0.1,
                   rho = 0.2,
                   h = NULL, # se NULL, fa tuning su griglia
                   tune_grid = NULL, # override griglia; altrimenti automatica
                   eps_inf = 0.05, # max % di PI infiniti ammessa nel tuning
                   seed = NULL,
                   i1 = NULL,
                   verbose = FALSE) {
  require(LCP)
  
  set.seed(seed)
  x = as.matrix(x)
  y = as.numeric(y)
  x0 = as.matrix(x0)
  n = nrow(x)
  stopifnot(length(y) == n)
  
  if (verbose == TRUE)
    txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  # split
  if (is.null(i1))
    i1 = sample(seq_len(n), floor(n * rho), replace = F)
  
  Xfit = x[i1, , drop = FALSE]
  yfit = y[i1]
  Xcal = x[-i1, , drop = FALSE]
  ycal = y[-i1]
  n1 = length(i1)
  n2 = n - n1
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n", txt, n1, n2))
    cat(sprintf("%sTraining on first part ...\n", txt))
  }
  # fit model
  out = train.fun(x = Xfit, y = yfit)
  mu_cal = predict.fun(out, newx = Xcal) |>
    as.vector()
  mu_te  = predict.fun(out, newx = x0) |>
    as.vector()
  
  if (verbose) {
    cat(sprintf(
      "%sComputing residuals and quantiles on second part ...\n",
      txt
    ))
  }
  #score
  Vcal = abs(ycal - mu_cal)
  
  # distance matrices (between calib-calib and test-calib)
  D_cal     = pairwise_dist(Xcal, Xcal)     # n_cal x n_cal
  D_testcal = pairwise_dist(x0 , Xcal)     # n_test x n_cal
  
  # ordered for LCP
  ord = order(Vcal)
  H_cal_ord  = D_cal[ord, ord]
  Vcal_ord   = Vcal[ord]
  H_new      = D_testcal[, ord, drop = FALSE]
  H_newT     = t(H_new)
  
  # LCP constructor given h, returns delta(x0)
  build_with_h = function(hval) {
    mod = LCPmodule$new(
      H = H_cal_ord,
      V = Vcal_ord,
      h = hval,
      alpha = alpha,
      type = "distance"
    )
    mod$lower_idx()
    mod$cumsum_unnormalized()
    mod$LCP_construction(Hnew = H_new, HnewT = H_newT)
    mod$band_V  #semi-wide vector (delta) may contain Inf!
  }
  
  chosen_h = h
  deltas = NULL
  #if not given h, tuning on grid
  if (is.null(chosen_h)) {
    base_scale = median(D_cal[D_cal > 0])
    if (is.na(base_scale) || base_scale == 0)
      base_scale = 1
    if (is.null(tune_grid))
      tune_grid = base_scale * c(0.25, 0.5, 1, 2, 4)
    
    best = list(score = Inf,
                h = tune_grid[1],
                deltas = NULL)
    for (hh in tune_grid) {
      dlt = build_with_h(hh)
      inf_rate = mean(is.infinite(dlt))
      #average length of finished PIs (2*delta) penalizes if above threshold infinite "eps_inf"
      avg_len = if (all(is.infinite(dlt)))
        Inf
      else
        mean(2 * dlt[is.finite(dlt)])
      penalized = if (inf_rate <= eps_inf)
        avg_len
      else
        avg_len + 1e6 * (inf_rate - eps_inf)
      if (penalized < best$score)
        best = list(score = penalized,
                    h = hh,
                    deltas = dlt)
    }
    chosen_h = best$h
    deltas   = best$deltas
  } else {
    deltas = build_with_h(chosen_h)
  }
  
  # test intervals: [muhat(x_te) - delta, muhat(x_te) + delta]
  lower = mu_te - deltas
  upper = mu_te + deltas
  
  data.frame(
    center = mu_te,
    delta  = deltas,
    lower  = lower,
    upper  = upper,
    width  = 2 * deltas,
    is_infinite = !is.finite(deltas),
    h = chosen_h,
    alpha = alpha,
    row.names = NULL
  )
}

#LCP (quantile-based score )
lcp_lcqr = function(x,
                    y,
                    x0,
                    alpha = 0.1,
                    rho = 0.2,
                    train.quant.fun,
                    predict.quant.fun,
                    q.probs = c(alpha / 2, 1 - alpha / 2),
                    h = NULL,
                    tune_grid = NULL,
                    eps_inf = 0.05,
                    seed = NULL,
                    i1 = NULL,
                    verbose = FALSE) {
  require(LCP)
  
  set.seed(seed)
  x  = as.matrix(x)
  y  = as.numeric(y)
  x0 = as.matrix(x0)
  n  = nrow(x)
  stopifnot(length(y) == n)
  
  if (verbose == TRUE)
    txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  
  # split
  if (is.null(i1))
    i1 = sample(seq_len(n), floor(n * rho), replace = FALSE)
  
  Xfit = x[i1, , drop = FALSE]
  yfit = y[i1]
  Xcal = x[-i1, , drop = FALSE]
  ycal = y[-i1]
  n1   = length(i1)
  n2   = n - n1
  
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n", txt, n1, n2))
    cat(sprintf("%sTraining quantile model on first part ...\n", txt))
  }
  # target quantiles 
  q.probs = sort(q.probs)
  # fit qr model
  out = train.quant.fun(x = Xfit, y = yfit, q = q.probs)
  # prediction quantiles
  Qcal = predict.quant.fun(out, newx = Xcal, q = q.probs)
  Qte  = predict.quant.fun(out, newx = x0 , q = q.probs)
  
  ql_cal = as.numeric(Qcal[, 1])
  qh_cal = as.numeric(Qcal[, 2])
  ql_te  = as.numeric(Qte [, 1])
  qh_te  = as.numeric(Qte [, 2])
  
  if (verbose) {
    cat(sprintf("%sComputing CQR residual scores on second part ...\n", txt))
  }
  
  # score CQR: V_i = max{ ql(x_i) - y_i, y_i - qh(x_i) }
  Vcal = pmax(ql_cal - ycal, ycal - qh_cal)
  
  # distance matrices (between calib-calib and test-calib)
  D_cal     = pairwise_dist(Xcal, Xcal)   # n_cal x n_cal
  D_testcal = pairwise_dist(x0 , Xcal)    # n_test  x n_cal
  
  # ordered for LCP
  ord        = order(Vcal)
  H_cal_ord  = D_cal[ord, ord]
  Vcal_ord   = Vcal[ord]
  H_new      = D_testcal[, ord, drop = FALSE]
  H_newT     = t(H_new)
  
  # LCP constructor given h, returns delta(x0)
  build_with_h = function(hval) {
    mod = LCPmodule$new(
      H = H_cal_ord,
      V = Vcal_ord,
      h = hval,
      alpha = alpha,
      type = "distance"
    )
    mod$lower_idx()
    mod$cumsum_unnormalized()
    mod$LCP_construction(Hnew = H_new, HnewT = H_newT)
    mod$band_V
  }
  
  # tuning of h if not fixed
  chosen_h = h
  deltas = NULL
  if (is.null(chosen_h)) {
    base_scale = median(D_cal[D_cal > 0])
    if (is.na(base_scale) || base_scale == 0)
      base_scale = 1
    if (is.null(tune_grid))
      tune_grid = base_scale * c(0.25, 0.5, 1, 2, 4)
    
    best = list(score = Inf,
                h = tune_grid[1],
                deltas = NULL)
    for (hh in tune_grid) {
      dlt = build_with_h(hh)
      inf_rate = mean(is.infinite(dlt))
      avg_len  = if (all(is.infinite(dlt)))
        Inf
      else
        mean((qh_te - ql_te)[is.finite(dlt)] + 2 * dlt[is.finite(dlt)])
      penalized = if (inf_rate <= eps_inf)
        avg_len
      else
        avg_len + 1e6 * (inf_rate - eps_inf)
      if (penalized < best$score)
        best = list(score = penalized,
                    h = hh,
                    deltas = dlt)
    }
    chosen_h = best$h
    deltas   = best$deltas
  } else {
    deltas = build_with_h(chosen_h)
  }
  
  # final intervals: [ q_lo(x0) - Δ(x0), q_hi(x0) + Δ(x0) ]
  lower = ql_te - deltas
  upper = qh_te + deltas
  
  data.frame(
    ql = ql_te,
    qh = qh_te,
    delta = deltas,
    lower = lower,
    upper = upper,
    width = upper - lower,
    center = (ql_te + qh_te) / 2,
    is_infinite = !is.finite(deltas),
    h = chosen_h,
    alpha = alpha,
    row.names = NULL
  )
}