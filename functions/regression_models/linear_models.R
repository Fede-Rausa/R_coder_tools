
## default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}


######### ols, wls from lm ---------------

model_wls = function(train, test, yname, 
                     error_fun=error_MAE, params=list()){
  
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}



######## glmnet lasso, ridge, elastic net ---------------

library(glmnet)

model_elnet = function(train, test, yname, 
                        error_fun=error_MAE, params=list(
                          alpha=NULL, w = NULL, lambda=0.5,
                          type='lasso', nlambda=1, standardize=F)){
  
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

  if (is.null(params$standardize)){
    params$standardize = F
  }
  
  if (is.null(params$lambda)){
    params$lambda = 0.01
  }
  
  inputvars = colnames(train)[colnames(train)!=yname]
  input = train[,inputvars]
  output = train[,yname]
  model = glmnet(input, output, weights = params$w,  lambda=params$lambda,
                 alpha=params$alpha, family='gaussian', 
                 standardize=params$standardize)
  
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
              test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
  
}



####### wls local regression ----------------------


model_localwls = function(train, test, yname,
                          error_fun=error_MAE, 
                          params=list(
                          predict_train = FALSE,
                          coordinates = NULL,
                          rm_coordinates = F,
                          distance='euclidean', 
                          wfun = NULL,
                          p = 2,
                          ecount=F)){
  
  
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
  if(is.null(params$ecount)){
    params$ecount = F
  }
  if(is.null(params$predict_train)){
    params$predict_train = F
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
    varinputs = colnames(train)[!(colnames(train) %in% params$coordinates)]
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
  
  X = cbind(rep(1, nrow(traindf)), 
            as.matrix(traindf[,colnames(traindf)!=yname]))
  Y = as.matrix(traindf[,yname])
  
  
  # Efficient WLS implementation

  # Ordinary least squares with Cholesky
  ols_chol <- function(X, y) {
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)
    R <- chol(XtX)
    #backsolve and forwardsolve are the functions
    #that solve uses
    beta_hat <- backsolve(R, forwardsolve(t(R), Xty))
    return(beta_hat)
  }

  ols_QR <- function(X, y) {
    qr_obj <- qr(X)
    qr.coef(qr_obj, y)
  }

  
  wls_solve <- function(X, y, w) {
    sqrt_w <- sqrt(w)
    
    X_weighted <- X * sqrt_w  # Broadcasting: each row multiplied by its weight
    y_weighted <- y * sqrt_w

    #B = ols_chol(X_weighted, y_weighted)
    #return(B)
    
    XtX <- crossprod(X_weighted)      # t(X_weighted) %*% X_weighted
    Xty <- crossprod(X_weighted, y_weighted)  # t(X_weighted) %*% y_weighted
    return(solve(XtX, Xty))
  }

  wls_QR = function(X,y,w){
    sqrt_w <- sqrt(w)
    X_weighted <- X * sqrt_w  # Broadcasting: each row multiplied by its weight
    y_weighted <- y * sqrt_w
    return(ols_QR(X_weighted, y_weighted))
  }

  # basic linear regression for error handling

  olsB = ols_QR(X, Y)
  
  
  pred = function(df, ecount=F){
    
    coordf0 = scale(df[,params$coordinates], 
                    center = means, scale = devstd)
    coordf1 = rbind(coordf, coordf0)
    rownames(coordf1) = 1:nrow(coordf1)
    
    dmat = as.matrix(dist(coordf1, method=params$distance, p=params$p))
    dmat0 = dmat[-(1:n),  -((n+1):ncol(dmat))]
    
    predf = as.matrix(cbind(rep(1, nrow(df)), 
                      df[, varinputs[varinputs!=yname]]))
    err_count = 0
    
    pred_wls = numeric(nrow(df))
    for (i in 1:nrow(df)){
      #D = dmat0[i,]
      #W = wfun(D)
      #B = wls_solve(X, Y, W)
      #row = as.matrix(predf[i,])
      #pred_wls[i] = t(B) %*% t(row)

      B = tryCatch(expr = wls_solve(X, Y, wfun(dmat0[i,])), 
                     error=function(msg){
                     err_count <<- err_count + 1
                     return(tryCatch(
                       expr = wls_QR(X, Y, wfun(dmat0[i,])),
                       error = function(msg){
                         err_count <<- err_count + 1
                         return(olsB)
                       }))})
      
      
      pred_wls[i] = t(B) %*% predf[i,]
    }
    if(params$ecount){
    print(paste0('errors count: ', err_count, '  errors ratio: ', err_count/nrow(df)))  
    }
    return(pred_wls)
  }
  

  if (params$predict_train){
    train_preds = pred(train)
    train_error = error_fun(train_preds, train[,yname])
    train_error = error_fun(train_preds, train[,yname])
    train_residuals = train_preds - train[,yname]
  }else{
    train_preds = NULL
    train_error = NULL
    train_error = NULL
    train_residuals = NULL
  }

  
  test_preds = pred(test)
  
  
  test_error = error_fun(test_preds, test[,yname])
  
  
  test_residuals = test_preds - test[,yname]

  return(list(predict=pred, 
             train_error=train_error,
             test_error=test_error,
             train_residuals=train_residuals,
             test_residuals=test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}




############ group penalized regression -----------------------------


library(grpreg)

model_grpreg = function(train, test, yname, error_fun, 
                        params=list(group = NULL,
                                    penalty='grLasso',
                                    sel = NULL,
                                    lambda=NULL)){
  
  if (is.null(params$group)){
    params$group = 1:(ncol(train)-1)
  }
  
  if (is.null(params$sel)){
    params$sel = colnames(train)[colnames(train)!=yname]
  }
  
  if (length(params$sel) != length(params$group)){
    stop('number of group members and predictors should be equal')
  }
  
  if (is.null(params$lambda)){
    params$lambda = 0.00000000001
  }
  
  
  model = grpreg(train[,params$sel], 
                 train[,yname], 
                 group=params$group, 
                 penalty=params$penalty,
                 lambda=params$lambda)
  
  
  pred = function(df){
    input = as.matrix(df[,params$sel])
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
              train_preds = train_preds,
              test_preds = test_preds,
              train_residuals=train_residuals,
              test_residuals=test_residuals))
  
}






