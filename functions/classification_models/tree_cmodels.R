
library(rpart)


####### classif tree ----------

model_ctree = function(train, test, yname, 
                     error_fun, params=list(
                       method='class', minbucket=20,
                       cp=0.01, xval=1
                     )){
  
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
  
  # if (params$method == 'class'){
  #   train[,yname] = factor(train[,yname], 
  #                           levels=unique(train[,yname]))
  #   test[,yname] = factor(test[,yname], 
  #                         levels=levels(train[,yname]))
  # }
  
  f = formula(paste0(yname, '~.'))
  model = rpart(f, data = train, method=params$method, weights=params$w,
                control = rpart.control(
                  minsplit = params$minsplit,
                  minbucket = params$minbucket,
                  cp = params$cp,
                  xval=params$xval))
  
  pred = function(df){
    if (params$method=='anova'){
      preds = predict(model, df)
    }else{
      preds = predict(model, df, type='class')
    }
    return(preds)
  }
  
  
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  if (params$method=='anova'){
    train_residuals = train_preds - train[,yname]
    test_residuals = test_preds - test[,yname]
  }else{
    train_residuals = NULL
    test_residuals = NULL
  }
  

  
  
  if (params$method=='anova'){
    id0 = round(train_preds)
  }else{
    id0 = train_preds
  }
  id1 = unique(id0)
  dict = rank(id1)
  names(dict) = id1
  
  
  pred_clu = function(df){
    if (params$method=='anova'){
      indexes = as.character(round(pred(df)))
    }else{
      indexes = as.character(pred(df))
    }
    clusters = dict[indexes]
    return(clusters)
  }
  
  if (params$method=='anova'){
    id0[id0<0] = 1
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
              test_residuals=test_residuals))
}



