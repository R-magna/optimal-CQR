#b.opt.dcp-----
#For each row of Q.yx, calculate the optimal bhat in [0, alpha] for Optimal DCP
b.opt.dcp = function(Q.yx, taus, alpha)
{
  b.grid = taus[taus <= alpha]
  u.grid = b.grid + 1 - alpha
  
  # Interpolazione su (b.grid + 1 - alpha) come in dcp.opt (functions_final.R) da https://github.com/kwuthrich/Replication_DCP
  Q.yx.u = sapply(seq_len(nrow(Q.yx)), function(t) {
    approx(x = taus, y = Q.yx[t, ], xout = u.grid, rule = 2)$y
  })
  if (length(b.grid) == 1) {
    Q.yx.u = matrix(Q.yx.u, ncol = 1)
  } else {
    Q.yx.u = t(Q.yx.u)
  }
  
  t.len = Q.yx.u - Q.yx[, seq_along(b.grid), drop = FALSE]
  t.pos = apply(t.len, 1, f.pos)
  return(b.grid[t.pos])
}





########DCP########
#Function to apply Split Distributional Conformal Prediction (DCP) and Optimal DCP
#(Chernozhukov, Wuthrich, Zhu, 2021) with code shared at the link https://github.com/kwuthrich/Replication_DC

#possible DCP variations:
# - tail.error.control = ("single", "separate", "both")
# - score.type = ("basic", "optimal", "both")
# - exact_coverage = (T, F) (if is T, and if it is necessary smoothed CP is applied)

#DCP.split-----
DCP.split = function(x,
                     y,
                     x0,
                     train.quant.fun,
                     predict.quant.fun,
                     alpha = 0.1,
                     rho = 0.5,
                     tail.error.control = "single",
                     score.type = "both",
                     tau_by = 0.005,
                     max.n.tau = 50,
                     taus = NULL,
                     exact_coverage = T,
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
  
  possible.tail = c("single", "separate", "both")
  possible.score = c("basic", "optimal", "both")
  
  # Check input arguments ----
  # Check rho
  if (is.null(rho) ||
      length(rho) != 1 || !is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("rho must be a number in between 0 and 1")
  
  # Check alpha
  if (is.null(alpha) ||
      length(alpha) != 1 || !is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be a number in between 0 and 1")
  
  # Build or check taus grid
  if (is.null(taus)) {
    if (is.null(tau_by) ||
        length(tau_by) != 1 || !is.numeric(tau_by) || tau_by <= 0 || tau_by >= 1)
      stop("tau_by must be a number in between 0 and 1")
    taus = round(seq(0, 1, by = tau_by), 6)
  } else {
    if (!is.numeric(taus) || any(taus <= 0) || any(taus >= 1))
      stop("taus must be a numeric vector with values strictly between 0 and 1")
    taus = sort(unique(taus))
  }
  
  # Check tail.error.control
  if (is.null(tail.error.control) ||
      !is.character(tail.error.control) ||
      !all(tolower(tail.error.control) %in% tolower(possible.tail)))
    stop(paste0(
      deparse(substitute(tail.error.control)),
      " must be one or more of the following:\n \"",
      paste(possible.tail, collapse = "\", \""),
      "\"."
    ))
  
  # Check score.type
  if (is.null(score.type) ||
      !is.character(score.type) ||
      !all(tolower(score.type) %in% tolower(possible.score)))
    stop(paste0(
      deparse(substitute(score.type)),
      " must be one or more of the following:\n \"",
      paste(possible.score, collapse = "\", \""),
      "\"."
    ))
  
  # Fix correct values of tail.error.control
  if (any(tolower(tail.error.control) == "both")) {
    tail.error.control = possible.tail[possible.tail != "both"]
  } else
    tail.error.control = possible.tail[tolower(possible.tail) %in% tolower(tail.error.control)]
  
  # Fix correct values of score.type
  if (any(tolower(score.type) == "both")) {
    score.type = possible.score[possible.score != "both"]
  } else
    score.type = possible.score[tolower(possible.score) %in% tolower(score.type)]
  
  if ("optimal" %in% score.type && !any(taus <= alpha))
    stop("For optimal DCP, taus must contain at least one value <= alpha")
  
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
  q.probs = c(alpha / 2, 1 - alpha / 2)
  q.train = sort(unique(c(taus, q.probs)))
  
  out = train.quant.fun(x[i1, , drop = F], y[i1], q = q.train)
  fit = matrix(predict.quant.fun(out, q = q.probs, x), nrow = n)
  
  if (verbose) {
    cat(sprintf(
      "%sComputing ranks, conformity scores and quantiles on second part ...\n",
      txt
    ))
  }
  
  # Get predicted quantile grid, rearrangement, and ranks on calibration & test----
  tQ.calib = matrix(
    predict.quant.fun(out, x[i2, , drop = F], q = taus),
    nrow = n2,
    ncol = length(taus)
  )
  # Rearrangement (monotone quantile curves)
  Q.calib = t(apply(tQ.calib, 1, FUN = sort))
  # Empirical CDF ranks \hat{U}_t = \hat{F}(Y_t, X_t) on calibration set
  u.hat = rowMeans(Q.calib <= y[i2])
  
  tQ.test = matrix(
    predict.quant.fun(out, x0, q = taus),
    nrow = nx0,
    ncol = length(taus)
  )
  Q.test = t(apply(tQ.test, 1, FUN = sort))
  
  if ("basic" %in% score.type)
  {
    # Two-column rank deviations around 0.5: max is |u.hat - 0.5|
    res = cbind(0.5 - u.hat, u.hat - 0.5)
    center.test = rep(0.5, nx0)
  }
  
  if ("optimal" %in% score.type)
  {
    # Estimate optimal shift \hat{b}(X, \alpha) in [0, \alpha]
    bhat.calib = b.opt.dcp(Q.yx = Q.calib, taus = taus, alpha = alpha)
    bhat.test  = b.opt.dcp(Q.yx = Q.test,  taus = taus, alpha = alpha)
    
    center.calib.optim = bhat.calib + (1 - alpha) / 2
    center.test.optim  = bhat.test  + (1 - alpha) / 2
    
    # Two-column rank deviations around bhat + (1 - alpha)/2: max is |u.hat - bhat - (1 - alpha)/2|
    res.optim = cbind(center.calib.optim - u.hat, u.hat - center.calib.optim)
  }
  
  #intervals----
  calculate.intervals = function(tail,
                                 score.variety,
                                 exact_coverage) {
    # Type of rank residuals and centers to be used
    switch(score.variety, "basic" = {
      score = res
      center = center.test
    }, "optimal" = {
      score = res.optim
      center = center.test.optim
    })
    
    # Append Inf so that k = ceiling((1 - alpha) * (n2 + 1)) in smoothed_cp
    score = rbind(score, c(Inf, Inf))
    
    # Quantile of interest
    switch(tail, "single" = {
      score = apply(score, MARGIN = 1, max)
      
      q = smoothed_cp(
        score = score,
        significance = 1 - alpha,
        n.test = nx0,
        exact_coverage = exact_coverage
      )
      q = cbind(q, q)
    }, "separate" = {
      q = cbind(
        smoothed_cp(
          score = score[, 1],
          significance = 1 - alpha / 2,
          n.test = nx0,
          exact_coverage = exact_coverage
        ),
        smoothed_cp(
          score = score[, 2],
          significance = 1 - alpha / 2,
          n.test = nx0,
          exact_coverage = exact_coverage
        )
      )
    })
    
    # Interval construction by inverting rank condition on taus grid
    Hat_dcp = sapply(seq_len(nx0), function(i) {
      cond = (center[i] - taus <= q[i, 1]) & (taus - center[i] <= q[i, 2])
      if (!any(cond)) {
        idx_closest = which.min(abs(taus - center[i]))
        c(Q.test[i, idx_closest], Q.test[i, idx_closest])
      } else {
        ci = Q.test[i, cond]
        c(min(ci), max(ci))
      }
    }) |>
      t()
    
    I.lo = Hat_dcp[, 1]
    I.up = Hat_dcp[, 2]
    
    name.pos = paste0(
      "DCP",
      ifelse(tail == "single", "", "_sep"),
      ifelse(score.variety == "basic", "", "_optim")
    )
    
    Int = list()
    Int[[name.pos]] = list(lo = I.lo, up = I.up)
    return(Int)
  }
  
  # DCP combinations to be calculated---
  tab = expand.grid(
    tail.error.control,
    score.type,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) |>
    as.matrix()
  
  if (nrow(tab) > 1)
    tab = tab[order(-xtfrm(tab[, 1]), tab[, 2]), , drop = FALSE]
  
  # Interval calculation
  Intervals = apply(tab, MARGIN = 1, function(X)
    calculate.intervals(
      tail = X[1],
      score.variety = X[2],
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

