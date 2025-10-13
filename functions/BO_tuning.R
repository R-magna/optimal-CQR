#########BO_tuning############
#Bayesian Optimization using mlr and mlrMBO libraries
#f_mlrMBO.bayesian.hyperpar.tuning----
f_mlrMBO.bayesian.hyperpar.tuning = function(data_mod,
                                             y,
                                             mod = "regr.ranger",
                                             k.fold.cv = 10,
                                             param.SET,
                                             par.fix = NULL,
                                             seed = 123,
                                             time.budget = 3600)
{
  require(mlr)
  require(mlrMBO)
  
  data_mod = cbind(data_mod, y = y)
  
  #def of task
  if (grepl("classif\\.", mod)) {
    #classification task
    if (is.factor(data_mod[, "y"]))
      data_mod[, "y"] = factor(data_mod[, "y"])
    task = makeClassifTask(data = data_mod, target = "y")
    predict.type = "prob"
    measures = list(acc, auc)
  } else if (grepl("regr\\.", mod)) {
    #regression task
    if (is.numeric(data_mod[, "y"]))
      data_mod[, "y"] = as.numeric(data_mod[, "y"])
    task = makeRegrTask(data = data_mod, target = "y")
    predict.type = "response"
    measures = list(rmse)
  }
  
  #param.SET
  d = length(ParamHelpers::getParamIds(param.SET))
  #def of learner
  if (length(par.fix) > 0) {
    lrn = makeLearner(mod, predict.type = predict.type, par.vals = par.fix)
  } else {
    lrn = makeLearner(mod, predict.type = predict.type)
  }
  #def tune control
  mbo.ctrl = makeMBOControl()
  mbo.ctrl = setMBOControlTermination(mbo.ctrl,
                                      iters = max(d * 10, 30),
                                      time.budget = time.budget)
  mbo.ctrl = setMBOControlInfill(mbo.ctrl, crit = makeMBOInfillCritEI())
  design_n = max(10, 3 * d)
  design = generateDesign(n = design_n, par.set = param.SET)
  ctrl = makeTuneControlMBO(mbo.control = mbo.ctrl, mbo.design = design)
  #def of cv
  rdesc = makeResampleDesc("CV", iters = k.fold.cv)
  set.seed(seed)
  #tuning
  tuned_params = tuneParams(
    learner = lrn,
    task = task,
    control = ctrl,
    par.set = param.SET,
    resampling = rdesc,
    measures = measures,
    show.info = T
  )
  return(tuned_params)
}
