

library(caret)


source('https://raw.githubusercontent.com/Fede-Rausa/R_coder_tools/refs/heads/main/functions/preprocess/TrainTestSplit.R')


bootstrap_error = function(df, yname, error_fun, model_fun, 
                      niter=3, tr_prp = 0.7, seed=42,
                      params=list()){
  
  
  testerror = numeric(niter)
  models = list()
  
  for (i in 1:niter){
    tts = train_test_split(df=df, tr_prp=tr_prp, yname=yname, seed=i+seed)
    
    models[[i]] = model_fun(tts$train, tts$test, 
                            yname=yname, error_fun=error_fun, 
                            params=params)
    
    testerror[i] = models[[i]]$test_error
  }
  
  
  return(list(bterror=testerror,
              meanerror = mean(testerror),
              models = models))
}




