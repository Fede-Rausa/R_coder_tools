

GridSearch = function(train, test, yname, error_fun, parList, model_fun){
  
  #example of parList: 
  #parList = list(k=c(1,3,4), j = c('dog', 'cat'), h = seq(1,5,0.5))

  parDF = expand.grid(parList)
  
  N = nrow(parDF)
  test_errors = train_errors = rep(NA, N)
  models = list()
  
  for (i in 1:N){
    params = as.list(parDF[i,])
    newmodel = model_fun(train=train, test=test, 
                         error_fun=error_fun, 
                         yname=yname, params=params)
    test_errors[i] = newmodel$test_error
    train_errors[i] = newmodel$train_error
    models[[i]] = newmodel
  }
  parDF[,'test_error'] = test_errors
  parDF[,'train_error'] = train_errors
  
  best_model = models[[which.min(test_errors)]]
  
  return(list(reportGrid = parDF, 
              best_model = best_model, 
              mod_list = models))
  
}
