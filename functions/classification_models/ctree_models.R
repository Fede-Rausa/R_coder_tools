
### default loss --------

error_misclass = function(preds, y){
  tab = table(preds, y)
  N = length(y)
  acc = sum(diag(tab))/N
  return(1 - acc)
}


####### classif tree ----------

library(rpart)

model_ctree = function(train, test, yname, 
                       error_fun = error_misclass, params=list(
                         method='class', minbucket=20,
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
  
  # if (params$method == 'class'){
  #   train[,yname] = factor(train[,yname], 
  #                           levels=unique(train[,yname]))
  #   test[,yname] = factor(test[,yname], 
  #                         levels=levels(train[,yname]))
  # }
  
  f = formula(paste0(yname, '~.'))
  model = rpart(f, data = train, method='class', 
                weights=params$w,
                control = rpart.control(
                  minsplit = params$minsplit,
                  minbucket = params$minbucket,
                  cp = params$cp,
                  xval=params$xval))
  
  pred = function(df){
    preds = predict(model, df, type='class')
    return(preds)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  id0 = train_preds
  id1 = unique(id0)
  dict = rank(id1)
  names(dict) = id1
  
  
  pred_clu = function(df){
    indexes = as.character(pred(df))
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
              test_error=test_error))
}



