############stability_selection_lasso###########
#f_stability.selection_lasso----
f_stability.selection_lasso = function(x.train,
                                       y.train,
                                       B = 10,
                                       lambda.min.ratio = 1e-4,
                                       n.lambda = 200,
                                       vect.lambda = NULL)
{
  require(glmnet)
  x.train = as.matrix(x.train)
  N = nrow(x.train)
  if (!is.factor(y.train)) {
    family = "gaussian"
    f.cal.err = function(x , y, pred)
      #mse
    {
      y.hat = pred[, x]
      return(mean((y.hat - y)^2))
    }
    f.cal.j.pos = function(mod, x)
    {
      j.pos = which(mod$beta[, x] != 0)
      vect.pos = rep(0, mod$dim[1])
      vect.pos[j.pos] = 1
      return(vect.pos)
    }
  }
  if (is.null(vect.lambda)) {
    lambda.max = max(abs(crossprod(
      scale(x.train, center = T, scale = T), y.train
    ))) / N
    lambda.min = lambda.max * lambda.min.ratio
    
    vect.lambda = seq(
      from = log(lambda.min),
      to = log(lambda.max * 1.01),
      length.out = n.lambda
    ) |>
      exp()
  } else {
    n.lambda = length(vect.lambda)
  }
  
  name.col = colnames(x.train)
  res = matrix(0, ncol = length(name.col) + 3, nrow = n.lambda)
  colnames(res) = c("lambda", "err", "num.selected", name.col)
  res = as_tibble(res)
  res[, "lambda"] = sort(vect.lambda, decreasing = T)
  s.vect = paste0("s", 0:(n.lambda - 1))
  
  for (b in 1:B)
  {
    fold.i = fold_assignment(N, K = 2)
    #data split
    X1 = x.train[fold.i == 1, ]
    y1 = y.train[fold.i == 1]
    
    X2 = x.train[fold.i == 2, ]
    y2 = y.train[fold.i == 2]
    
    #model estimation for each lambda
    mod = glmnet(
      X1,
      y1,
      family = family,
      alpha = 1,
      lambda = vect.lambda,
      maxit = 500000
    )
    #model predictions and error
    pred = predict(mod, newx = X2, type = "response")
    res[, "err"] = res[, "err"] + sapply(s.vect, f.cal.err, y = y2, pred = pred)
    
    #saving selected variables
    res[, name.col] = res[, name.col] + t(sapply(s.vect, f.cal.j.pos, mod =
                                                   mod))
    
    #same thing for the second data partition
    mod = glmnet(
      X2,
      y2,
      family = family,
      alpha = 1,
      lambda = vect.lambda,
      maxit = 500000
    )
    pred = predict(mod, newx = X1, type = "response")
    res[, "err"] = res[, "err"] + sapply(s.vect, f.cal.err, y = y1, pred = pred)
    res[, name.col] = res[, name.col] + t(sapply(s.vect, f.cal.j.pos, mod =
                                                   mod))
  }
  
  res$num.selected = rowSums(res[, name.col])
  res[, -1] = res[, -1] / B / 2
  return(res)
}