
######### ols, wls from lm ---------------

model_wls = function(train, test, yname, 
                     error_fun, params=list()){
  
  N = nrow(train)
  
  if(is.null(params$w)){
    params$w = rep(1, N)
  }
  
  f = formula(paste0(yname, '~.'))
  model = lm(f, data=train, weights=params$w)
  
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
              train_error=train_error, 
              test_error=test_error,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
}



######## glmnet lasso, ridge, elastic net ---------------

library(glmnet)

model_elnet = function(train, test, yname, 
                        error_fun, params=list(
                          alpha=NULL, w = NULL, lambda=0.5,
                          type='lasso', nlambda=1)){
  
  if (is.null(params$alpha)){
    if (params$type=='lasso'){
      params$alpha=1
    }else{
      if (params$type=='ridge'){
        params$alpha=0
      }else{
        if (params$type=='elnet'){
          params$alpha=0.5
        }
      }
    }
  }
  
  if (is.null(params$lambda)){
    params$lambda = 0.01
  }
  
  inputvars = colnames(train)[colnames(train)!=yname]
  input = train[,inputvars]
  output = train[,yname]
  model = glmnet(input, output, weights = params$w,  lambda=params$lambda,
                 alpha=params$alpha, family='gaussian')
  
  pred = function(df){
    input = as.matrix(df[,inputvars])
    preds = predict(model, input)
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



####### wls local regression ----------------------


model_localwls = function(train, test, yname,
                          error_fun, 
                          params=list(
                          coordinates = NULL,
                          rm_coordinates = F,
                          distance='euclidean', 
                          wfun = NULL,
                          p = 2)){
  
  
  if (is.null(params$coordinates)){
    params$coordinates = colnames(train)[which(colnames(train)!=yname)]
  }
  if (is.null(params$rm_coordinates)){
    params$rm_coordinates = F
  }
  if (is.null(params$wfun)){
    wfun = function(d){return(20*(max(d) - d)^10)}
  }else{
    wfun = params$wfun
  }
  
  if (is.null(params$distance)){
    params$distance = 'euclidean'
  }else{
    if(!(params$distance %in% c("euclidean", "maximum", 
                                "manhattan", "canberra", 
                                "binary", "minkowski"))){
    params$distance = 'euclidean'
    }
  }
      
    
  
  if (!params$rm_coordinates){
    varinputs = colnames(train)
  }else{
    varinputs = colnames(train)[colnames(train) %in% params$coordinates]
  }
  
  if (is.null(params$p)){
    params$p = 2
  }
  
  
  coordf = scale(train[,params$coordinates])
  means = attr(coordf,"scaled:center")
  devstd = attr(coordf,"scaled:scale")
  
  traindf = train[,varinputs]
  f = formula(paste0(yname, '~.'))
  n = nrow(train)
  
  X = as.matrix(traindf[,colnames(traindf)!=yname])
  Y = as.matrix(traindf[,yname])
  
  
  # Efficient WLS implementation
  wls_solve <- function(X, y, w) {
    sqrt_w <- sqrt(w)
    
    X_weighted <- X * sqrt_w  # Broadcasting: each row multiplied by its weight
    y_weighted <- y * sqrt_w

    XtX <- crossprod(X_weighted)      # t(X_weighted) %*% X_weighted
    Xty <- crossprod(X_weighted, y_weighted)  # t(X_weighted) %*% y_weighted
    return(solve(XtX, Xty))
  }
  
  
  pred = function(df){
    
    coordf0 = scale(df[,params$coordinates], 
                    center = means, scale = devstd)
    coordf1 = rbind(coordf, coordf0)
    rownames(coordf1) = 1:nrow(coordf1)
    
    dmat = as.matrix(dist(coordf1, method=params$distance, p=params$p))
    dmat0 = dmat[-(1:n),  -((n+1):ncol(dmat))]
    
    predf = df[, varinputs[varinputs!=yname]]
    
    pred_wls = numeric(nrow(df))
    for (i in 1:nrow(df)){
      #D = dmat0[i,]
      #W = wfun(D)
      #B = wls_solve(X, Y, W)
      #row = as.matrix(predf[i,])
      #pred_wls[i] = t(B) %*% t(row)
      pred_wls[i] = t(wls_solve(X, Y, wfun(dmat0[i,]))) %*% t(as.matrix(predf[i,]))
    }

    return(pred_wls)
  }
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]

  return(list(predict=pred, 
             train_error=train_error,
             test_error=test_error,
             train_residuals=train_residuals,
             test_residuals=test_residuals))
}

