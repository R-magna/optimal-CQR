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
  if(all(c("lo","up") %in% names(out.cp)))
  {
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
  } else #CASE CTI
  {
     tab = sapply(seq_len(length(out.cp)), function(i) {
      Ints = out.cp[[i]]
      y.i = y.true[i]
      width.i = 0
      coverage.i = F
      is.out.dx.i = NA
      is.out.sx.i = NA
     
      for(I in Ints)
      {
        lo.i = I[1]
        up.i = I[2]
        if (!is.null(y_min))
        {
          lo.i = ifelse(lo.i < y_min, y_min, lo.i)
          up.i = ifelse(up.i < y_min, y_min, up.i)
        }
        if (!is.null(y_max))
        {
          lo.i = ifelse(lo.i > y_max, y_max, lo.i)
          up.i = ifelse(up.i > y_max, y_max, up.i)
        }
        width.i = width.i + up.i - lo.i
        if(y.i >= lo.i & y.i <= up.i)
        {
          coverage.i = T
          is.out.dx.i = F
          is.out.sx.i = F
        }
      }
      I.edges = unlist(Ints)
      delta.i = min(abs(I.edges - y.i))
      if(!coverage.i)
      {
        if(y.i > max(I.edges))
        {
          is.out.dx.i = T
          is.out.sx.i = F
        } else if(y.i < min(I.edges))
        {
          is.out.dx.i = F
          is.out.sx.i = T
        }
      }
      
      w.score.i = width.i + 2 /alpha * ifelse(coverage.i, 0, delta.i)
      c(
        LO = NA,
        UP = NA,
        coverage = coverage.i,
        is.out.dx = is.out.dx.i,
        is.out.sx = is.out.sx.i,
        width = width.i,
        w.score = w.score.i
      )
    })|>
       t()|>
       as_tibble()
  }
  
  return(tab)
}

#write_AVG.files-----
#write average result of comparison CQR as a CSV file
write_AVG.files = function(paths = NULL)
{
  require(dplyr)
  require(data.table)
  if (is.null(paths))
    paths = list.dirs("comparison_CQR/", recursive = F)
  paths = sub("/$", "", paths)
  
  for (path in paths)
  {
    type = strsplit(path, "/")[[1]][2]
    path_avg = paste0(path, "/AVG_result/")
    if (!dir.exists(path_avg))
      dir.create(path_avg, recursive = TRUE)
    for (path_dataset in setdiff(list.dirs(path, recursive = F), paste0(path, "/AVG_result")))
    {
      dt.name = strsplit(x = path_dataset, split = "/")[[1]] |>
        tail(1)
      path.file.avg = paste0(path_avg, dt.name, "._.avg", ".csv")
      if (!file.exists(path.file.avg))
      {
        tab_avg = c()
        for (file.name in list.files(path_dataset,
                                     pattern = "\\.csv$",
                                     recursive = F))
        {
          rep = ((strsplit(file.name, "rep_")[[1]][2]) |>
                   strsplit(split = "\\."))[[1]][1]
          dt = fread(paste0(path_dataset, "/", file.name))
          tab_avg = rbind(tab_avg,
                          cbind(
                            rep = rep,
                            select(
                              dt,
                              dataset,
                              type,
                              is.out.dx,
                              is.out.sx,
                              coverage,
                              width,
                              w.score
                            ) |>
                              group_by(dataset, type) |>
                              summarise_all(mean)
                          ))
        }
        fwrite(tab_avg, file = path.file.avg)
      }
    }
    #create finals tab
    path_final = paste0(path_avg, "final_tab/")
    if (!dir.exists(path_final))
      dir.create(path_final, recursive = TRUE)
    
    tab_avg = c()
    for (file.name in list.files(path_avg, pattern = "\\.csv$", recursive = F))
    {
      dt.name = strsplit(x = file.name, split = "\\._.")[[1]][1]
      dt = fread(paste0(path_avg, file.name))
      if (!"dataset" %in% colnames(dt))
        dt =  cbind(dataset = dt.name, dt)
      tab_avg = rbind(tab_avg, dt)
    }
    
    fwrite(tab_avg,
           file = paste0(path_final, type, "._.", "all_rep_result", ".csv"))
  }
}

