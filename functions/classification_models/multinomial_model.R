
#### multinomial logistic regression ---------------

library(VGAM)

model_multilogit = function(train, test, yname, error_fun,
                            params=list()){
  
  library(VGAM)
  f = paste0(yname, '~.')
  
  model = vglm(f, data=train, family=multinomial)
  
  pred = function(df){
    probsmat = predict(model, df, type='response')
    preds = colnames(probsmat)[apply(probsmat, 1, which.max)]
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
  
  #pred_clu = pred
  #train_clu = train_preds
  #test_clu = test_preds
  
  return(list(model=model, 
              predict=pred,
              pred_clu = pred_clu,
              train_clu = train_clu,
              test_clu = test_clu,
              train_error=train_error, 
              test_error=test_error))
  
}



