
###### regression losses --------------

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}

error_MSE = function(preds, y){
  return(mean((preds - y)^2))
}

error_RMSE = function(preds, y){
  return(sqrt(mean((preds - y)^2)))
}

error_R2 = function(preds, y){
  return(1 - sum((preds - y)^2)/sum((y - mean(y))^2))
}

error_MAPE = function(preds, y){
  return(mean(abs((preds - y)/y)))
}


##### classification losses ----------------

error_misclass = function(preds, y){
  tab = table(preds, y)
  N = length(y)
  acc = sum(diag(tab))/N
  return(1 - acc)
}



