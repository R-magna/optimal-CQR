#####IMPORT DATASETS######
#function to import and transform benchmark datasets.
#It is an adaptation in R, with some modifications, of the Python function:
#"GetDataset" of https://github.com/msesia/cqr-comparison.git of the paper
#Sesia, M., & Candès, E. J. (2020). A comparison of some conformal quantile regression methods.
#Stat, 9(1), e261.
#(is also added the Superconductivity Data)

#import_dataset------
import_dataset = function(data.name,
                          base.path = "datasets/",
                          factor_to_dummy = T) {
  require(readr)
  require(dplyr)
  if (factor_to_dummy)
    require(fastDummies)
  
  name.col.target = c()
  name.col.remove = c()
  name.col.factor = c()
  
  if (data.name == "bike")
  {
    data = read_csv(paste0(base.path, "bike.csv"))
    
    name.col.target = "cnt"
    name.col.remove = c("instant", "dteday", "casual", "registered")
    name.col.factor = c("season",
                        "yr",
                        "weathersit",
                        "holiday",
                        "weekday",
                        "workingday")
    
  } else if (data.name == "bio")
  {
    data = read_csv(paste0(base.path, "CASP.csv"))
    name.col.target = "RMSD"
    
  } else if (data.name == "blog")
  {
    data = read_csv(paste0(base.path, "blogData_train.csv"), col_names = FALSE)
    
    n.uniq_j = apply(data, MARGIN = 2, function(x)
      length(unique(x)))
    
    #Remove variables with identical values (zero variability)
    name.col.remove = names(which(n.uniq_j == 1))
    data = select(data, -all_of(name.col.remove))
    name.col.remove = c()
    
    #Remove identical variables (keep only 1 variable per group of identical variables)
    sigs = sapply(data, function(col)
      paste(col, collapse = "\r"))
    
    groups_euqal.var = split(names(sigs), sigs) |>
      unname() |>
      Filter(
        f = function(x)
          length(x) > 1
      )
    
    name.col.remove = sapply(groups_euqal.var, `[`, -1) |>
      unlist()
    data = select(data, -all_of(name.col.remove))
    name.col.remove = c()
    
    #Remove a dummy variable for each type
    dummies_j = list(paste0("X", 263:269), paste0("X", 270:276))
    for (facotor.group in seq_len(length(dummies_j)))
    {
      var.f = dummies_j[[facotor.group]]
      #check
      if (all(rowSums(data[, var.f]) == 1))
      {
        name.col.remove = c(name.col.remove, var.f[1])
      }
    }
    
    name.col.target = "X281"
    
  } else if (data.name == "community")
  {
    data = read_csv(paste0(base.path, "communities.data"), col_names = FALSE)
    lines = readLines(paste0(base.path, "communities.names"))
    
    colnames(data) = lines[grepl(x = lines, pattern = "^@attribute")] |>
      strsplit(" ") |>
      sapply(`[`, 2)
    
    data[data == "?"] = NA
    
    #non-informative variables are removed
    name.col.remove = c("state", "county", "community", "communityname", "fold")
    
    data = select(data, -all_of(name.col.remove)) |>
      mutate(across(everything(), as.numeric))
    name.col.remove = c()
    
    #the only missing value is imputed
    data$OtherPerCap[is.na(data$OtherPerCap)] = mean(data$OtherPerCap, na.rm = T)
    
    n.na_j = apply(data, MARGIN = 2, function(x)
      sum(is.na(x)))
    
    #all variables with missing values are removed
    name.col.remove = colnames(data)[n.na_j > 0]
    name.col.target = "ViolentCrimesPerPop"
    
  } else if (data.name == "concrete")
  {
    data = read_csv(paste0(base.path, "concrete.csv"))
    
    colnames(data) = strsplit(colnames(data), "\\(") |>
      sapply(`[`, 1) |>
      trimws(which = "right") |>
      gsub(pattern = " ", replacement = "_")
    
    name.col.target = colnames(data)[ncol(data)]
  } else if (data.name %in% paste0("facebook_", 1:2))
  {
    n.varaint = strsplit(data.name, split = "_")[[1]][2]
    
    data = read_csv(paste0(base.path, "Features_Variant_", n.varaint, ".csv"),
                    col_names = FALSE)
    
    #Remove variables with identical values (zero variability)
    n.uniq_j = apply(data, MARGIN = 2, function(x)
      length(unique(x)))
    name.col.remove = colnames(data)[n.uniq_j == 1]
    
    #Remove a dummy variable for each type
    dummies_j = list(paste0("X", 40:46), paste0("X", 47:53))
    for (facotor.group in seq_len(length(dummies_j)))
    {
      var.f = dummies_j[[facotor.group]]
      #check
      if (all(rowSums(data[, var.f]) == 1))
      {
        name.col.remove = c(name.col.remove, var.f[1])
      }
    }
    
    name.col.target = colnames(data)[ncol(data)]
    
  } else if (data.name == "homes")
  {
    data = read_csv(paste0(base.path, "kc_house_data.csv"))
    
    name.col.target = "price"
    name.col.remove = c("id", "date")
    
  } else if (data.name %in% paste0("meps_", 19:21))
  {
    data = read_csv(paste0(base.path, data.name , "_reg.csv"))[, -1]
    
    colnames(data) = colnames(data) |>
      gsub(pattern = "=", replacement = "_") |>
      gsub(pattern = "-", replacement = ".")
    
    data = relocate(data, RACE, .after = last_col())
    
    #Remove a dummy variable for each type
    uniq.var = which((
      colnames(data) |>
        gsub(patter = "[0-9_.]", replacement = "") |>
        table()
    ) > 1) |>
      names()
    name.col.remove = c()
    for (var_u in uniq.var)
    {
      sub.dt = select(data, starts_with(var_u))
      #check
      if (all(rowSums(sub.dt) == 1))
      {
        name.col.remove = append(name.col.remove, colnames(sub.dt)[1])
      }
    }
    data = select(data, -all_of(name.col.remove))
    name.col.remove = c()
    #Remove identical variables (keep only 1 variable per group of identical variables)
    sigs = sapply(data, function(col)
      paste(col, collapse = "\r"))
    
    groups_euqal.var = split(names(sigs), sigs) |>
      unname() |>
      Filter(
        f = function(x)
          length(x) > 1
      )
    
    name.col.remove = sapply(groups_euqal.var, `[`, -1) |>
      unlist()
    
    name.col.target = "UTILIZATION_reg"
    
  } else if (data.name == "star")
  {
    data = read_csv(paste0(base.path, "STAR.csv"))[, -1]
    data = na.omit(data)
    
    encoding.var = list(
      gender     = c(female = 0, male = 1),
      ethnicity  = c(
        cauc = 0,
        afam = 1,
        asian = 2,
        hispanic = 3,
        amindian = 4,
        other = 5
      ),
      star       = c(
        regular = 0,
        small = 1,
        `regular+aide` = 2
      ),
      lunch      = c(free = 0, `non-free` = 1),
      school     = c(
        `inner-city` = 0,
        suburban = 1,
        rural = 2,
        urban = 3
      ),
      degree     = c(
        bachelor = 0,
        master = 1,
        specialist = 2,
        `master+` = 3,
        phd = 3
      ),
      ladder     = c(
        level1 = 0,
        level2 = 1,
        level3 = 2,
        apprentice = 3,
        probation = 4,
        pending = 5,
        noladder = 5,
        notladder = 6
      ),
      tethnicity = c(
        cauc = 0,
        afam = 1,
        asian = 2
      )
    )
    for (root.var in names(encoding.var))
    {
      poss.var = c(root.var, paste0(root.var, "k"), paste0(root.var, 1:3))
      name.var = poss.var[poss.var %in% colnames(data)]
      
      mapping = encoding.var[[root.var]]
      data[, name.var] = apply(data[, name.var], 2, function(x)
        mapping[unlist(x)] |>
          as.numeric())
    }
    
    grade.col = sapply(c("read", "math"), function(i)
      c(paste0(i, "k"), paste0(i, 1:3))) |>
      as.vector()
    
    data = mutate(data, grade = rowSums(data[, grade.col]))
    
    name.col.target = "grade"
    name.col.remove = grade.col
  } else if (data.name == "super")
  {
    data = read_csv(paste0(base.path, "super_train.csv"))
    molecules = read_csv(paste0(base.path,"super_unique_m.csv"))
    data = mutate(
      data,
      G_Iron = as.numeric(molecules$Fe > 0),
      G_Cuprate = as.numeric(molecules$Cu > 0 &
                               molecules$O > 0)
    ) |>
      relocate(critical_temp, .after = last_col())
    rm(molecules)

    name.col.target = "critical_temp"
  }
  
  Y = data[[name.col.target]]
  X = select(data, -all_of(c(name.col.remove, name.col.target))) |>
    mutate(across(all_of(name.col.factor), factor)) |>
    mutate(across(where(is.character), factor))
  
  if (factor_to_dummy & any(sapply(X, is.factor)))
    X = dummy_columns(X,
                      remove_first_dummy = T,
                      remove_selected_columns = T)
  return(list(X = X, Y = Y))
}

#fast.check.df-----
fast.check.df = function(df)
{
  n.uniq_j = apply(df, MARGIN = 2, function(x)
    length(unique(x)))
  frac_na_j = colMeans(is.na(df))
  frac_na_i = rowMeans(is.na(df))
  pos.duplicate_j = which(duplicated(as.list(df)) |
                            duplicated(as.list(df), fromLast = TRUE))
  pos.duplicate_i = which(duplicated(df) |
                            duplicated(df, fromLast = TRUE))
  
  if (any(n.uniq_j == 1))
    print(paste(sum(n.uniq_j == 1), "variables with only 1 value"))
  if (any(n.uniq_j < 10))
    print("possible factor variables present")
  if (anyNA(df)) {
    print(paste(sum(is.na(df)), "NAs present"))
    print(paste(max(frac_na_j), "maximum fraction of NA for variables"))
    print(paste(max(frac_na_i), "maximum fraction of NA for observations"))
  }
  if (length(pos.duplicate_j))
  {
    print(paste(length(pos.duplicate_j), "equal variables"))
  }
  if (length(pos.duplicate_i))
  {
    print(paste(length(pos.duplicate_i), "equal obserbation"))
  }
  return(
    list(
      n.uniq_j = n.uniq_j,
      frac_na_j = frac_na_j,
      frac_na_i = frac_na_i,
      pos.duplicate_j = pos.duplicate_j,
      pos.duplicate_i = pos.duplicate_i
    )
  )
}