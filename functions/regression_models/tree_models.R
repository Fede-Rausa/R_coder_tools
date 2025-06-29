## default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}

######### rpart regression tree ----------------------------


model_rtree = function(train, test, yname, 
                      error_fun=error_MAE, params=list(
                        minbucket=20,
                        cp=0.01, xval=1
                      )){
  
  N = nrow(train)
  library(rpart)
  

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
  if (is.null(params$max_depth)){
    params$max_depth = 30
  }
  
  f = formula(paste0(yname, '~.'))
  model = rpart(f, data = train, method='anova', 
                weights=params$w,
                control = rpart.control(
                  minsplit = params$minsplit,
                  minbucket = params$minbucket,
                  maxdepth = params$max_depth,
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
  id0[id0<0] = 1
  id1 = unique(id0)
  dict = rank(id1)
  names(dict) = id1
  
  
  pred_clu = function(df){
    indexes = as.character(round(pred(df)))
    clusters = dict[indexes]
    return(clusters)
  }
  
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}




########## rf ------------

library(randomForest)
library(ranger)

model_ranger = function(train, test, yname, 
                    error_fun, params=list()){
  
  
  if (is.null(params$ntree)){
    params$ntree = 100
  }
  if (is.null(params$mtry)){
    params$mtry = floor(sqrt(ncol(train)-1))
  }
  if (is.null(params$minbucket)){
    params$minbucket = 1
  }
  if (is.null(params$max.depth)){
    params$max.depth = 0
  }
  if (is.null(params$min.node.size)){
    params$min.node.size = 5
  }
  if (is.null(params$replace)){
    params$replace = T
  }
  if (is.null(params$sample.fraction)){
    params$sample.fraction = ifelse(params$replace, 1, 0.632)
  }
  if (is.null(params$fixsplitvars)){
    params$fixsplitvars = NULL
  }
  
  
  model = ranger(as.formula(paste0(yname, '~.')), 
                 data=train, 
                 num.trees= params$ntree,
                 mtry = params$mtry,
                 min.bucket = params$minbucket,
                 max.depth = params$max.depth,
                 min.node.size = params$min.node.size,
                 sample.fraction = params$sample.fraction,
                 replace = params$replace,
                 always.split.variables = params$fixsplitvars,
                 classification = F)
    
  
  pred = function(df){
    preds = predict(model, df)
    return(preds$predictions)
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}


model_randomForest = function(train, test, yname, 
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}





########## gbm -------------------------

library(gbm)

model_gbm = function(train, test, yname, 
                    error_fun, params=list()){

  if (is.null(params$ntree)){
    params$ntree = 100
  }
  if (is.null(params$max_depth)){
    params$max_depth = 1
  }
  if (is.null(params$lr)){
    params$lr = 0.1
  }
  if (is.null(params$min.node.size)){
    params$min.node.size = 1
  }
  if (is.null(params$sample.fraction)){
    params$sample.fraction = 0.5
  }



  f = as.formula(paste0(yname, '~.'))
  model = gbm(f, train,
              n.trees = params$ntree,
              shrinkage = params$lr,
              interaction.depth = params$max_depth,
              bag.fraction = params$sample.fraction,
              distribution = 'gaussian')

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

  return(list(model=model,
              predict=pred,
              train_preds = train_preds,
              test_preds = test_preds,
              train_error=train_error,
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
}







########### xgb ---------------

library(xgboost)
model_xgb = function(train, test, yname, 
                    error_fun=error_MAE, params=list()){
  
  
  
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}








