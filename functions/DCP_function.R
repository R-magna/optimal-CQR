#######DCP#############
#Application of Distributional Conformal Prediction.

#The functions were obtained by adapting the relative R functions of:
#Chernozhukov, V., W¨uthrich, K., & Zhu, Y. (2021). Distributional conformal prediction.
#Proceedings of the National Academy of Sciences, 118 (48), e2107794118.
#https://github.com/kwuthrich/Replication_DCP/tree/main

#dcp.qr-----
dcp.qr = function(Y0,
                  X0,
                  Y1,
                  X1,
                  Y.test,
                  X.test,
                  taus,
                  alpha.sig,
                  score.type = "both")
{
  require(quantreg)
  possible.score = c("basic", "optimal", "both")
  
  if (is.null(score.type) ||
      !is.character(score.type) ||
      !all(tolower(score.type) %in% tolower(possible.score)))
    stop(paste0(
      deparse(substitute(score.type)),
      " must be one or more of the following:\n \"",
      paste(possible.score, collapse = "\", \""),
      "\"."
    ))
  #Fix correct values of score.type
  if (any(tolower(score.type) == "both")) {
    score.type = possible.score[possible.score != "both"]
  } else
    score.type = possible.score[tolower(possible.score) %in% tolower(score.type)]
  
  
  X0 = as.matrix(X0)
  X1 = as.matrix(X1)
  X.test = as.matrix(X.test)
  Y0 = as.numeric(Y0)
  Y1 = as.numeric(Y1)
  
  
  beta.qr = matrix(NA, dim(X0)[2] + 1, length(taus))
  for (t in 1:length(taus)) {
    beta.qr[, t] = rq.fit.br(cbind(1, X0), Y0, tau = taus[t])$coefficients
  }
  Intervals = list()
  if ("basic" %in% score.type)
  {
    tQ.yx   = cbind(1, X1) %*% beta.qr
    Q.yx    = t(apply(tQ.yx, 1, FUN = sort))
    u.hat   = rowMeans((Q.yx <= matrix(Y1, length(Y1), length(taus))))
    cs      = abs(u.hat - 0.5)
    
    tQ.test   = cbind(1, X.test) %*% beta.qr
    Q.test    = t(apply(tQ.test, 1, FUN = sort))
    u.test    = rowMeans((Q.test  <= matrix(Y.test, length(Y.test), length(taus))))
    cs.test   = abs(u.test - 0.5)
    
    k         = ceiling((1 - alpha.sig) * (1 + length(Y1)))
    threshold = sort(cs)[k]
    
    ci.grid   = abs(taus - 0.5)
    #cov.qr    = (cs.test <= threshold)
    
    lb = ub = rep(NA, length(Y.test))
    for (i in 1:length(Y.test)) {
      ci     = Q.test[i, (ci.grid <= threshold)]
      ub[i]  = max(ci)
      lb[i]  = min(ci)
    }
    Intervals[["DCP.qr"]] = list(lo = lb, up = ub)
    #leng.qr = ub-lb
    #leng.qr[which(leng.qr==-Inf)] = NA
    
  }
  if ("optimal" %in% score.type)
  {
    XXX = rbind(X1, X.test)
    YYY = c(Y1, Y.test)
    tQ.yx   = cbind(1, XXX) %*% beta.qr
    Q.yx    = t(apply(tQ.yx, 1, FUN = sort))
    u.hat = rowMeans((Q.yx <= matrix(YYY, length(YYY), dim(beta.qr)[2])))
    
    bhat = rep(NA, length(YYY))
    b.grid = taus[taus <= alpha.sig]
    for (t in 1:length(YYY)) {
      leng = rep(NA, length(b.grid))
      leng.test = rep(NA, length(b.grid))
      for (b in 1:length(b.grid)) {
        Q.yx.u = approx(
          x = taus,
          y = Q.yx[t, ],
          xout = (b.grid[b] + 1 - alpha.sig),
          rule = 2
        )$y
        leng[b] = Q.yx.u - Q.yx[t, b]
      }
      bhat[t] = b.grid[which.min(leng)]
    }
    
    ind.test = (length(Y1) + 1):length(YYY)
    
    cs.opt = abs(u.hat - bhat - (1 - alpha.sig) / 2)
    
    k           = ceiling((1 - alpha.sig) * (1 + length(Y1)))
    threshold   = sort(cs.opt[-ind.test])[k]
    
    cov.opt   = (cs.opt[ind.test] <= threshold)
    
    # leng.opt = NULL
    ub = c()
    lb = c()
    for (t in ind.test) {
      ci.grid = abs(taus - bhat[t] - (1 - alpha.sig) / 2)
      ci = Q.yx[t, (ci.grid <= threshold)]
      ub = append(ub, max(ci))
      lb = append(lb, min(ci))
      # leng.opt = c(leng.opt,ub-lb)
    }
    Intervals[["DCP.qr_optim"]] = list(lo = lb, up = ub)
    #leng.opt[which(leng.opt==-Inf)] = NA
  }
  return(Intervals) #cov.qr=cov.qr,leng.qr=leng.qr
}



#dcp.dr-----
dcp.dr = function(Y0,
                  X0,
                  Y1,
                  X1,
                  Y.test,
                  X.test,
                  ys,
                  taus,
                  alpha.sig)
{
  X0 = as.matrix(X0)
  X1 = as.matrix(X1)
  X.test = as.matrix(X.test)
  Y0 = as.numeric(Y0)
  Y1 = as.numeric(Y1)
  
  beta.dr = matrix(NA, dim(X0)[2] + 1, length(ys))
  for (y in 1:length(ys)) {
    beta.dr[, y] = glm.fit(cbind(1, X0), (Y0 <= ys[y]), family = binomial(link =
                                                                            "logit"))$coefficients
  }
  tF.yx   = plogis(cbind(1, X1) %*% beta.dr)
  F.yx    = t(apply(tF.yx, 1, FUN = sort))
  
  cs = rep(NA, length(Y1))
  for (t in 1:length(Y1)) {
    u.hat = approx(
      x = ys,
      y = F.yx[t, ],
      xout = Y1[t],
      rule = 2
    )$y
    cs[t] = abs(u.hat - 0.5)
  }
  
  tF.test   = plogis(cbind(1, X.test) %*% beta.dr)
  F.test    = t(apply(tF.test, 1, FUN = sort))
  
  cs.test = rep(NA, length(Y.test))
  for (t in 1:length(Y.test)) {
    u.hat = approx(
      x = ys,
      y = F.test[t, ],
      xout = Y.test[t],
      rule = 2
    )$y
    cs.test[t] = abs(u.hat - 0.5)
  }
  
  k           = ceiling((1 - alpha.sig) * (1 + length(Y1)))
  threshold   = sort(cs)[k]
  ci.grid   = abs(taus - 0.5)
  #cov.dr    = (cs.test <= threshold)
  
  lb = ub = rep(NA, length(Y.test))
  for (i in 1:length(Y.test)) {
    ci     = ys[abs(F.test[i, ] - 0.5) <= threshold]
    ub[i]  = max(ci)
    lb[i]  = min(ci)
  }
  #leng.dr = ub-lb
  #leng.dr[which(leng.dr==-Inf)] = NA
  return(list(lo = lb, up = ub))
}
