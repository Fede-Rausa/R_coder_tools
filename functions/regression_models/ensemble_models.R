
### default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}


######### boosting --------------

# model_fun should be a function from the script supervised_regressions.R or clustering_models.R


boost_model = function(train, test, yname, 
                       error_fun = error_MAE, model_fun,  
                       nrounds = 20, lr = 0.1, stocastic = F, 
                       tr_prp = 0.5, params=list()){
  
  ## function to perform least squares boosting
  ## model_fun should be a function that takes as input train, test, yname, 
  ## error_fun and params, and returns a list with an attribute predict

  train0 = train
  n = nrow(train)
  train_preds = numeric(n)
  res = train0[,yname] 
  ##attention! Never start by substracting the mean. It won't work
  
  models = list()
  boost_error = numeric(nrounds)
  tau = round(n*tr_prp)
  
  
  for (m in 1:nrounds){
    
    if (stocastic){
      train1 = train0[sample(1:n, tau, replace=F),]
    }else{
      train1 = train0
    }
    
    models[[m]] = model_fun(train=train1, test=train1,
                            yname=yname, error_fun=error_fun,
                            params=params)$predict
    new_preds = models[[m]](train0)
    res = res - lr*new_preds
    train0[,yname] = res
    boost_error[m] = mean(abs(res))
  }
  
  
  pred = function(df){
    preds = numeric(nrow(df))
    for (m in 1:nrounds){
      preds = preds + lr*models[[m]](df)
    }
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  return(list(predict=pred,
              train_error=train_error,
              test_error = test_error,
              train_residuals = train_residuals,
              test_residuals = test_residuals,
              boost_error = boost_error))
  
}


############ bagging ---------------


bagging_model = function(train, test, yname, model_fun, 
                         error_fun= = error_MAE,
                         nrounds=10, replace=T,
                         nobs = NULL, nvars=NULL, 
                         params=list()){
  
  p = ncol(train) - 1
  n = nrow(train)
  if (is.null(nobs)){
    nobs=n
  }
  
  id_y = which(colnames(train)==yname)
  models = list()
  if (!is.null(nvars)){
    nx = min(nvars, p)
  }
  
  for (m in 1:nrounds){
    id_obs = sample(1:n, nobs, replace=replace)
    train0 = train[id_obs,]
    if (is.null(nvars)){
      train1 = train0
    }else{
      inputs = train0[,-id_y]
      id_x = sample(1:p, nx)
      train1 = cbind(train0[,id_x])
      train1[,yname] = train0[,yname]
    }
    
    models[[m]] = model_fun(train=train1, test=train1,
                            yname=yname, error_fun=error_fun,
                            params=params)$predict
  }
  
  
  pred = function(df){
    np = nrow(df)
    mat = matrix(NA, nrow=np, ncol=nrounds)
    
    for (m in 1:nrounds){
      mat[,m] = models[[m]](df)
    }
    preds = rowMeans(mat, na.rm = T)
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  return(list(predict=pred,
              train_error=train_error,
              test_error = test_error,
              train_residuals = train_residuals,
              test_residuals = test_residuals))
  
}




######### bagging + boosting or boosting + bagging --------------

bagboost_model = function(train, test, yname, error_fun = error_MAE, 
                          model_fun, params=list(),
                          first_bagging = T,
                          nrounds_bag=10, replace=T,
                          nobs = NULL, nvars=NULL,
                          nrounds_boost=20, 
                          lr=0.1, stocastic = F, 
                          tr_prp=0.5){
  
  if (first_bagging){
    
    my_model_fun = function(train, test, yname, error_fun, params){
      f = bagging_model(train, test, 
                        yname, model_fun, error_fun,
                        nrounds_bag, replace, 
                        nobs, nvars,params)
      return(f)
    }
    
    boosting = boost_model(train, test, yname, error_fun, my_model_fun,
                           nrounds_boost, lr, stocastic, tr_prp, params)
    
    return(boosting)
    
  }else{
    
    my_model_fun = function(train, test, yname, error_fun, params){
      f = boost_model(train, test, yname, error_fun, model_fun,
                      nrounds_boost, lr, stocastic, tr_prp, params)
      return(f)
    }
    
    bagging = bagging_model(train, test, 
                            yname, my_model_fun, error_fun,
                            nrounds_bag, replace, 
                            nobs, nvars,params)
    return(bagging)
  }
}


########## clustering + supervised --------------------------

ensemble_cluster = function(train, test, yname, error_fun,
                            cluster_fun, super_fun,
                            params=list()){
  
  #cluster_fun is a function that returns an object with the
  #attribute pred_clu, or a function that return directly pred_clu function
  #it should also keep an attribute predict that returns the mean of
  #the y inside the cluster
  
  #super_fun is a classical function for supervised model
  #and it should return an object with the attribute predict
  #or directly the function
  
  #both should receive train, yname and params
  
  
  model_clu = cluster_fun(train=train, test=test, yname=yname, 
                          params=params, error_fun=error_fun)
  pred_clu = model_clu$pred_clu
  
  clusters0 = pred_clu(train)
  c1 = unique(clusters0)
  k = length(c1)
  
  smodels = list()
  
  for (i in 1:k){
    train0 = train[clusters0==c1[i],]
    smodels[[i]] = super_fun(train=train0, test=test, yname=yname, 
                             params=params, error_fun=error_fun)
  }
  
  pred = function(df){
    clusters = pred_clu(df)
    preds = numeric(length(clusters))
    for (i in 1:k){
      filter = (clusters==i)
      if (sum(filter)>0){
        preds[filter] = smodels[[i]]$predict(df[filter,])
      }
    }
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  train_clu = pred_clu(train)
  test_clu = pred_clu(test)
  
  return(list(model_sup=smodels,
              model_clu=model_clu,
              predict=pred,
              pred_clu = pred_clu,
              train_clu = train_clu,
              test_clu = test_clu,
              train_error=train_error,
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
  
}




