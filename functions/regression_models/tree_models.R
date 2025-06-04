

########## rpart CART tree -------------

library(rpart)

model_tree = function(train, test, yname, 
                     error_fun, params=list(
                       method='anova', minbucket=20,
                       cp=0.01, xval=1
                     )){
  #mibucket is the minimum number of observations that should fall in a leaf
  N = nrow(train)
  
  if(is.null(params$method)){
    params$method='anova'
  }
  if(is.null(params$minsplit)){
    params$minsplit=20
  }
  if (is.null(params$minbucket)){
    params$minbucket= ceiling(params$minsplit/3)
  }
  if(is.null(params$cp)){
    params$cp=0.01
  }
  if(is.null(params$xval)){
    params$xval=1
  }
  
  f = formula(paste0(yname, '~.'))
  model = rpart(f, data = train, method=params$method, weights=params$w,
                control = rpart.control(
                  minsplit = params$minsplit,
                  minbucket = params$minbucket,
                  cp = params$cp,
                  xval=params$xval))
  
  pred = function(df){
    preds = predict(model, df)
    return(preds)
  }
  
  
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  
  
  id0 = round(train_preds)
  id1 = unique(id0)
  dict = rank(id1)
  names(dict) = id1
  
  
  pred_clu = function(df){
    indexes = as.character(round(pred(df)))
    clusters = dict[indexes]
    return(clusters)
  }
  
  
  id0[id0<0] = 1

  train_clu = dict[id0]
  test_clu = pred_clu(test)
  
  return(list(model=model, 
              predict=pred,
              pred_clu = pred_clu,
              train_clu = train_clu,
              test_clu = test_clu,
              train_error=train_error, 
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
}


########## rf ------------

library(randomForest)

model_rf = function(train, test, yname, 
                    error_fun, params=list()){
  
  
  
  if (is.null(params$ntree)){
    params$ntree = 100
  }
  if (is.null(params$mtry)){
    params$mtry = floor(sqrt(ncol(train)-1))
  }
  if (is.null(params$importance)){
    params$importance = TRUE
  }
  
  model = randomForest(as.formula(paste0(yname, '~.')), 
                       data=train, 
                       ntree=params$ntree, 
                       mtry=params$mtry, 
                       importance=params$importance,
                       na.action = na.omit, 
                       replace = F)
  
  pred = function(df){
    preds = predict(model, newdata=df)
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  return(list(model=model, 
              predict=pred, 
              train_error=train_error, 
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
}


########### xgb ---------------

library(xgboost)
model_xgb = function(train, test, yname, 
                    error_fun, params=list()){
  
  
  
  if (is.null(params$nrounds)){
    params$nrounds = 100
  }
  if (is.null(params$max_depth)){
    params$max_depth = 6
  }
  if (is.null(params$eta)){
    params$eta = 0.1
  }
  
  model = xgboost(params = params,
                 data = as.matrix(train[, -which(colnames(train) == yname)]),
                 label = train[,yname],
                 nrounds = params$nrounds,
                 early_stopping_rounds = 10,
                 verbose = 0)
  
  
  pred = function(df){
    ddf = as.matrix(df[, -which(colnames(df) == yname)])
    preds = predict(model, ddf)
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  return(list(model=model, 
              predict=pred, 
              train_error=train_error, 
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
}








