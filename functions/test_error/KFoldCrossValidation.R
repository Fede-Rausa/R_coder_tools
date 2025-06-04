

library(caret)

cv_error = function(df, yname, error_fun, model_fun, kf = 3, 
                    seed=42, params=list()){
  
  target = df[,yname]
  
  if (!is.null(params$seed)){
    set.seed(seed)
    folds = createFolds(target, kf)
  }else{
    folds = createFolds(target, kf)
  }
  
  testerror = numeric(kf)
  models = list()
  
  for (i in 1:kf){
    
    idf = folds[[i]]
    trainset = df[-idf,]
    testset = df[idf,]
    models[[i]] = model_fun(train=trainset, test=testset, 
                            yname=yname, error_fun=error_fun,
                            params=params)
    testerror[i] = models[[i]]$test_error
  }
  
  
  
  return(list(cverror=testerror,
              meanerror = mean(testerror),
              models = models))
}


