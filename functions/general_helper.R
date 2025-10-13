#########general_helper_functions###########

#fold_assignment------
#Assign each observation to a fold in (1,..,K) in a fair manner
fold_assignment = function(n.obs, K)
{
  ii = rep(NA, n.obs)
  obs_for_fold = floor(n.obs / K)
  for (fold in seq_len(K))
  {
    ii[sample(which(is.na(ii)), replace = F, size = obs_for_fold)] = fold
  }
  #assigning missing observations to folds
  ii[is.na(ii)] = sample(1:K, replace = F, size = sum(is.na(ii)))
  return(ii)
}

#winkler_score------------
#calculates the Winkler score of the prediction intervals
winkler_score = function(lo, up, y.true, alpha)
{
  width = up - lo
  res = (lo - y.true) * (y.true < lo) + (up - y.true) * (y.true > up)
  return(width + 2 / alpha * res)
}
#create_tab_out.cp----
#given the output of a CP procedure, constructs the results table
create_tab_out.cp = function(out.cp,
                             y.true,
                             alpha = 0.1,
                             y_min = NULL,
                             y_max = NULL)
{
  if (!require("dplyr", quietly = TRUE)) {
    stop("Package dplyr not installed (required here)!")
  }
  tab = tibble(LO = as.vector(out.cp$lo),
               UP = as.vector(out.cp$up),
               y = y.true)
  if (!is.null(y_min))
    tab = mutate(tab,
                 LO = ifelse(LO < y_min, y_min, LO),
                 UP = ifelse(UP < y_min, y_min, UP))
  if (!is.null(y_max))
    tab = mutate(tab,
                 LO = ifelse(LO > y_max, y_max, LO),
                 UP = ifelse(UP > y_max, y_max, UP))
  tab = mutate(
    tab,
    coverage = y <= UP & y >= LO,
    is.out.dx = ifelse(y > UP, T, F),
    is.out.sx = ifelse(y < LO, T, F),
    width = abs(UP - LO),
    w.score = width + 2 / alpha * case_when(is.out.dx ~ abs(y - UP), is.out.sx ~ abs(LO - y), T ~ 0)
  ) |>
    select(-y)
  return(tab)
}

