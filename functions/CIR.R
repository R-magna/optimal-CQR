########CIR########
#Function to apply CIR an CTI.
#The function was written maintaining a structure similar to that of the
#“conformal.pred.split” function (for Slit CP) of "conformalInference" package.
#https://github.com/ryantibs/conformal/tree/master/conformalInference

#possible variations:
# - method = ("CIR", "CTI")
# - exact_coverage = (T, F) (if is T, and if it is necessary smoothed CP is applied)
# - T.intervals = numeric (defines the number of subintervals defined by estimating T.intervals+1 conditional quantiles)

#CIR-----
CIR.split = function(
    x,
    y,
    x0,
    train.quant.fun,
    predict.quant.fun,
    alpha = 0.1,
    rho = 0.5,
    method = c("CIR","CIR_mio","CTI"),
    T.intervals = 100,
    exact_coverage = T,
    split = NULL,
    seed = NULL,
    verbose = FALSE,
    calculate_fit = T)
{
  require(data.table)
  # Set up data----
  x = as.data.table(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = as.data.table(x0, ncol = p)
  nx0 = nrow(x0)
  
  possible.method = c("CIR","CIR_mio","CTI","all") #CHR, CTI, CIR+, CIR++ !!!CIR++ è una mia variazione
  # Check input arguments ----
  #conformalInference:::check.args(x=x,y=y,x0=x0,alpha=alpha,train.fun=train.quant.fun,
  #                                predict.fun=predict.quant.fun)
  # Check rho
  if (is.null(rho) ||
      length(rho) != 1 || !is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("rho must be a number in between 0 and 1")

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

  #Fix correct values of method
  if (any(tolower(method) == "all")) {
    method = possible.method[possible.method != "all"]
  } else
    method = possible.method[tolower(possible.method) %in% tolower(method)]
  
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
  q.probs = seq(0,1, by=1/T.intervals)
  Intervals = list()
  #set.seed(123)
  out = train.quant.fun(x[i1, , drop = F], y[i1], q = q.probs)
  if(calculate_fit)
    fit = matrix(predict.quant.fun(out, q = q.probs, x), nrow = n)
  if (verbose) {
    cat(sprintf(
      "%sComputing residuals and quantiles on second part ...\n",
      txt
    ))
  }
  
  # Get score and quantiles, and prediction on second----
  pred.calib = matrix(predict.quant.fun(out, x[i2, , drop = F], q = q.probs),
                      nrow = n2,
                      ncol = T.intervals+1)
  length_int.calib = pred.calib[, -1] - pred.calib[, -ncol(pred.calib)]
  y_interval.calib =  lapply(seq_len(n2), function(i) {
    which(pred.calib[i,-ncol(pred.calib)] <= y[i2][i] & pred.calib[i,-1] >= y[i2][i])
  })
  
  #prediction intervals
  pred.test = matrix(predict.quant.fun(out, x0, q = q.probs),
                     nrow = nx0,
                     ncol = T.intervals+1)
  length_int.test = pred.test[, -1] - pred.test[, -ncol(pred.test)]
  
  if(any( c("CIR", "CIR_mio") %in% method))
  {
    #ATTENZIONE: !!!uguale al codice CIR!!!
    #in caso di piu intervalli di stessa lunghezza minima si prende il primo
    pos.start_calib = apply(X = length_int.calib, MARGIN = 1, which.min)
    pos_start.test = apply(X = length_int.test, MARGIN = 1, which.min)
  }
  ###CIR_mio
  if("CIR_mio" %in% method)
  {
    score.calib_continuous =  sapply(seq_len(n2), function(i) {
      I.start = pos.start_calib[i]
      I.len = length_int.calib[i,]
      I.y = y_interval.calib[[i]]
      
      k = 1
      C.k_min = I.start
      C.k_max = I.start
      I.last = I.start
      
      if(length(I.y)==0)
      {
        return(T.intervals)
      } else if(I.start %in% I.y)
      {
        return(1)
      }
      while (k < T.intervals)
      {
        if(C.k_min==1)
        {
          k = k + min(I.y) - C.k_max
          C.k_max = min(I.y)
          I.last = min(I.y)
          break
        } else if(C.k_max==T.intervals)
        {
          k = k + max(I.y) - C.k_min
          C.k_min = max(I.y)
          I.last = max(I.y)
          break
        } else
        {
          k = k+1
          if(I.len[C.k_min-1]<=I.len[C.k_max+1])
          {
            C.k_min = C.k_min - 1
            I.last = C.k_min
          } else
          {
            C.k_max = C.k_max + 1
            I.last = C.k_max
          }
        }
        if(any(I.y %in% C.k_min:C.k_max))
        {
          break
        }
      }
      if(I.last == C.k_min)
      {
        fraz.score = (pred.calib[i,C.k_min+1] - y[i2][i]) / (pred.calib[i,C.k_min+1] - pred.calib[i,C.k_min])
      } else if(I.last == C.k_max)
      {
        fraz.score = (y[i2][i] - pred.calib[i,C.k_max]) / (pred.calib[i,C.k_max+1] - pred.calib[i,C.k_max])
      } else
      {
        print("errore: I.last non è un sottointervallo estremo")
      }
      if(fraz.score<0 | fraz.score>1)
      {
        print("errore: fraz.score<0 | fraz.score>1")
      }
      k = k-1 + fraz.score
      return(k)
    })
    score.calib_continuous = c(score.calib_continuous, Inf)
    q_hat_continuous = smoothed_cp(
      score = score.calib_continuous,
      significance = 1 - alpha,
      n.test = nx0,
      exact_coverage = exact_coverage
    )
    Int.predict =  sapply(seq_len(nx0), function(i) {
      I.start = pos_start.test[i]
      I.len = length_int.test[i,]
      k.left = ceiling(q_hat_continuous[i,]) - 1
      fraz.q_hat_continuous = q_hat_continuous[i,] - k.left
      
      C.k_min = I.start
      C.k_max = I.start
      I.last = I.start
      
      while (k.left > 0)
      {
        if(C.k_min==1)
        {
          C.k_max = C.k_max + k.left
          k.left = 0 
          I.last = C.k_max
        } else if(C.k_max==T.intervals)
        {
          C.k_min = C.k_min - k.left
          k.left = 0 
          I.last = C.k_min
        } else
        {
          k.left = k.left - 1
          if(I.len[C.k_min-1]<=I.len[C.k_max+1])
          {
            C.k_min = C.k_min - 1
            I.last = C.k_min
          } else
          {
            C.k_max = C.k_max + 1 
            I.last = C.k_max
          }
        }
      }
      lo = pred.test[i,C.k_min]
      up = pred.test[i,C.k_max+1]
      
      I.last.len = length_int.test[i,I.last]
      if(fraz.q_hat_continuous > 0 & I.last.len>0)
      {
        delta = (1-fraz.q_hat_continuous) * I.last.len
        if (I.last == C.k_min)
        { 
          lo = lo + delta
        } else if (I.last == C.k_max)
        {
          up = up - delta
        }
      }
      return(c(lo,up))
    })|>
      t()
    Intervals[["CIR_mio"]] = list(lo = Int.predict[,1],up = Int.predict[,2])
  }
  
  #ATTENZIONE: !!!uguale al codice CIR!!!
  #1) nel caso in cui la lunghezza dell intervallo di sx è uguale a quella dell intervallo di dx
  #si sceglie l intervalli di sx
  #2) nel caso in cui la vera y è minore del quantile stimato di livello 0, o maggiore del quantile
  #stimato di livello 1, si uniscono tutti gli intervalli, e quindi si ha score = T
  if("CIR" %in% method)
  {
    if("CIR_mio" %in% method)
    {
      score.calib = ceiling(score.calib_continuous)
    } else
    {
      score.calib =  sapply(seq_len(n2), function(i) {
        I.start = pos.start_calib[i]
        I.len = length_int.calib[i,]
        I.y = y_interval.calib[[i]]
        
        k = 1
        C.k_min = I.start
        C.k_max = I.start
        
        if(length(I.y)== 0)
        {
          return(T.intervals)
        }
        else if(I.start %in% I.y)
        {
          return(1)
        }
        while (k < T.intervals)
        {
          if(C.k_min==1)
          {
            k = k + min(I.y)- C.k_max
            break
          } else if(C.k_max==T.intervals)
          {
            k = k + max(I.y) - C.k_min
            break
          } else
          {
            k = k+1
            if(I.len[C.k_min-1]<=I.len[C.k_max+1])
            {
              C.k_min = C.k_min - 1
            } else
            {
              C.k_max = C.k_max + 1 
            }
          }
          if(any(I.y %in% C.k_min:C.k_max))
          {
            break
          }
        }
        return(k)
      })
    }
    score.calib = c(score.calib, Inf)
    q_hat = smoothed_cp(
      score = score.calib,
      significance = 1 - alpha,
      n.test = nx0,
      exact_coverage = exact_coverage
    )
    Int.predict =  sapply(seq_len(nx0), function(i) {
      I.start = pos_start.test[i]
      I.len = length_int.test[i,]
      k.left = q_hat[i,] - 1
      C.k_min = I.start
      C.k_max = I.start
      
      while (k.left > 0)
      {
        if(C.k_min==1)
        {
          C.k_max = C.k_max + k.left
          k.left = 0 
        } else if(C.k_max==T.intervals)
        {
          C.k_min = C.k_min - k.left
          k.left = 0 
        } else
        {
          k.left = k.left - 1
          if(I.len[C.k_min-1]<=I.len[C.k_max+1])
          {
            C.k_min = C.k_min - 1
          } else
          {
            C.k_max = C.k_max + 1 
          }
        }
      }
      return(pred.test[i,c(C.k_min,C.k_max+1)])
    })|>
      t()
    Intervals[["CIR"]] = list(lo = Int.predict[,1],up = Int.predict[,2])
  }
  
  if ("CTI" %in% method)
  {
    score.calib =  sapply(seq_len(n2), function(i) {
      I.len = length_int.calib[i, ]
      I.y = y_interval.calib[[i]]
      return(min(I.len[I.y]))
    })
    score.calib = c(score.calib, Inf)
    
    table(round(score.calib,2))
    q_hat = smoothed_cp(
      score = score.calib,
      significance = 1 - alpha,
      n.test = nx0,
      exact_coverage = exact_coverage
    )
    
    Int.predict =  lapply(seq_len(nx0), function(i) {
      I.len = length_int.test[i, ]
      t = q_hat[i, ]
      I.edges = pred.test[i, ]
      C.I.pos = which(I.len <= t)
      C_x = list()
      
      #BORDERLINE CASE no subinterval with length < t 
      if (length(C.I.pos) == 0)
      {
        return(C_x) #Returns an empty list
      }
      
      C.start = I.edges[C.I.pos[1]]
      C.end = I.edges[tail(C.I.pos, 1) + 1]
      
      pos_disjointed = which(diff(C.I.pos) > 1)
      
      C_l = C.start
      for (j in pos_disjointed)
      {
        C_u = I.edges[C.I.pos[j] + 1]
        C_x = append(C_x, list(c(C_l, C_u)))
        C_l = I.edges[C.I.pos[j + 1]]
      }
      C_u = C.end
      C_x = append(C_x, list(c(C_l, C_u)))
      return(C_x)
    })
    Intervals[["CTI"]] = Int.predict
    
  }
  
  
  result = list(
    Intervals = Intervals,
    split = i1
  )
  if(calculate_fit)
    result[["fit"]] = fit
  
  return(result)
}


