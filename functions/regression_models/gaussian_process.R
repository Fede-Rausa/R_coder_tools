## default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}


######### gaussian process ---------------

model_gp = function(train, test, yname, 
                     error_fun=error_MAE, params=list()){
  
  N = nrow(train)
  
  if(is.null(params$lambda)){
    params$lambda = 0
  }
  
  if (is.null(params$kernel)){
    rbfk = function(x,x0){
      return(exp(-(sum(abs(x-x0)^2)/2)))
    }
    params$kernel = rbfk 
  }
  
  
  
  
  x = train[,colnames(train)!=yname]
  xnames = colnames(x)
  x = as.matrix(x)
  y = unlist(train[,yname])
  vy = var(y)
  

  
  
  mykernel = params$kernel
  
  
  kmat = matrix(nrow=N, ncol=N)
  for (i in 1:nrow(kmat)){
    for (j in 1:ncol(kmat)){
      kmat[i,j] = mykernel(x[i,], x[j,])
    }
  }
  
  
  A = kmat + diag(N)*params$lambda
  #covmat = solve(A)
  L = chol(A)
  #P = t(solve(L))
  P = t(backsolve(L, diag(nrow(L))))
  covmat = t(P)%*%P
  alpha = covmat %*% y
  
  
  meangp = function(x0){
    m = 0
    for (i in 1:N){
      m=m+ mykernel(x[i,],x0)*alpha[i]
    }
    return(m)
  }
  
  
  vargp = function(x0){
    kvec = numeric(N)
    for (i in 1:N){
      kvec[i] = mykernel(x0,x[i,])
    }
    
    return(mykernel(x0,x0) - t(kvec) %*% covmat %*% kvec)
  }
  
  
  pred = function(df){
    xdf = as.matrix(df[,xnames])
    n = nrow(xdf)
    preds = numeric(n)
    for (i in 1:n){
      preds[i] = meangp(xdf[i,])
    }
    return(preds)
  }
  
  predvar = function(df){
    xdf = as.matrix(df[,xnames])
    n = nrow(xdf)
    vars = numeric(n)
    for (i in 1:n){
      vars[i] = vargp(xdf[i,])
    }
    return(vars)
  }
  
  
  model = list(meangp=meangp, vargp=vargp, 
               A=A, L=L, P=P, covmat=covmat, alpha=alpha, 
               mykernel=mykernel)
  
  
  train_preds = pred(train)
  test_preds = pred(test)
  
  train_error = error_fun(train_preds, train[,yname])
  test_error = error_fun(test_preds, test[,yname])
  
  train_residuals = train_preds - train[,yname]
  test_residuals = test_preds - test[,yname]
  
  return(list(model = model, 
              predict = pred,
              predvar = predvar,
              train_error = train_error, 
              test_error = test_error,
              train_residuals = train_residuals,
              test_residuals = test_residuals,
              train_preds = train_preds,
              test_preds = test_preds))
}

