



##this is the last official versione of the MSR model.
##It can be optimized with SGD methods or BFGS.
##It can handle L1 or L2 penalties, local or global.
##It is always convenient to scale the input variables before training
##a constant is always added

MSR_model = function(train, yname, k=2, lr=0.01, niter=10,
                     momentum=0, decay=0, train_type='BFGS',
                     penalty_L1 = FALSE, penalty_local=FALSE,
                     init_W = NULL, verbose=0, get_errors=FALSE,
                     winit = 0.001){
  
  y = train[,yname]
  dfx = train %>% 
    dplyr::select(-yname) %>% 
    as.matrix() %>%
    scale() 

  mycenters = attr(dfx,"scaled:center")
  myscales = attr(dfx,"scaled:scale")
  
 dfx=dfx %>% 
    as.data.frame() %>%
    mutate(const=1) %>% 
    as.matrix()
  
  xnames = colnames(dfx)
  
  
  N = length(y)
  p = ncol(dfx)
  
  if (is.null(init_W)){
    W = rnorm(p*k) * winit
    W = matrix(W, ncol=k)
  }else{
    if (ncol(init_W)!=k){
      stop(paste0('number of columns of init_W should be k=',k))
    }
    if (nrows(init_W)!=p){
      stop(paste0('number of rows of init_W should be p=',p))
    }
    
    W <- init_W
  }
  
  if (verbose){
    get_errors = TRUE
  }
  
  if (get_errors){
    ERRORS = c()
  }
  
  
  #D = exp(dfx%*%W)  #the D matrix is the focal point of training
  
  
  pred = function(df){ #W is a matrix pxk
    dfx = df %>%
      as.matrix() %>%
      scale(center=mycenters, scale=myscales) %>% 
      as.data.frame() %>%
      mutate(const=1) %>% 
      dplyr::select(xnames) %>% 
      as.matrix()
    
    ##compute pdf matrix
    D = exp(dfx%*%W)
    
    N = nrow(dfx)
    ##compute preds
    yh = numeric(N)
    dens = rowSums(D)
    for (j in 1:k){
      yh = yh + D[,j]*ksmean[j]
    }
    yh = yh/dens
    
    return(yh)
  }
  
  
  
  ##setup gradient function
  gd_gradient = function(W){
    
    #W = matrix(W, ncol=k)
    #k = ncol(W)
    #p = nrow(W)
    #N = nrow(dfx)
    
    ##compute pdf matrix
    D = exp(dfx%*%W)
    
    ##compute means
    ksmean = numeric(k)
    for (j in 1:k){
      d = D[,j]
      ksmean[j] = sum(y*d)/sum(d)
    }
    
    ##compute preds
    yh = numeric(N)
    dens = rowSums(D)
    for (j in 1:k){
      yh = yh + D[,j]*ksmean[j]
    }
    yh = yh/dens
    
    
    if (get_errors){
      er = mean((y-yh)^2)
      ERRORS <<- c(ERRORS, er)
      if (verbose){
        print(paste0('mse: ', er))
        #print(paste0('kms: ', paste(round(ksmean,2), collapse='-')))
      }
    }

    
    
    ##compute gradients
    Wgrad = matrix(0, nrow=nrow(W), ncol=ncol(W))
    res = (y - yh) * (-2)/N
    for (j in 1:k){
      
      nu_j = sum(y*D[,j])
      de_j = sum(D[,j])
      
      de_j_grad = t(dfx)%*%D[,j]
      nu_j_grad = t(dfx)%*%(D[,j]*y)
      
      mj_grad = (nu_j_grad*de_j - de_j_grad*nu_j)/(de_j^2)
      
      
      for (i in 1:N){
        
        nu_i = sum(ksmean * D[i,])
        de_i = sum(D[i,])
        
        dij_grad = dfx[i,]*D[i,j]
        
        yh_i_grad = (de_i*(dij_grad*ksmean[j] + mj_grad*D[i,j]) - nu_i*dij_grad)/(de_i^2)
        
        Wgrad[,j] = Wgrad[,j] + res[i]*yh_i_grad
        
      }
      
    }
    
    ##returns a matrix
    return(Wgrad)
  }
  
  
  if (decay > 0){
    
    if (penalty_L1){
      if (penalty_local){
        penalty_grad = function(W){return(sign(W))}
      }else{
        penalty_grad = function(W){return(sign(W)*sum(abs(W)))}
      }
    }else{
      if (penalty_local){
        penalty_grad = function(W){return(W)}
      }else{
        penalty_grad = function(W){return(sum(W))}
      }
    }
    
    main_gradient = function(W){
      Wgrad = gd_gradient(W)
      Wgrad = Wgrad - decay * penalty_grad(W)
      return(Wgrad)
    }
    
  }else{
    main_gradient = gd_gradient #gets a matrix and returns a matrix
  }
  

  
  
  bfgs_gradient = function(W){ # gets a vector and returns a vector
    W = matrix(W, ncol=k)
    Wgrad = main_gradient(W)
    Wgrad = as.vector(unlist(Wgrad))
    return(Wgrad)
  }
  
  
  
  mseloss = function(W){
    W = matrix(W, ncol=k)
    
    D = exp(dfx%*%W)
    
    km = (t(D)%*%y)/colSums(D)
    
    ypreds = (D%*%km)/rowSums(D)
    
    mse = mean((y-ypreds)^2)
    
    return(mse)
  }
  
  
  ####real training
  
  
  
  if (train_type=='BFGS'){
    
    W = as.vector(unlist(W))
    
    opt = optim(W, fn = mseloss, 
                gr = bfgs_gradient, 
                method='BFGS', 
                control=list(maxit=niter))

    W = opt$par
    W = matrix(W, ncol=k)
    
  }else{
    if (train_type=='SGD'){
      W = matrix(W, ncol=k)
      mom = momentum
      vel = 0
      for (t in 1:niter){
        g = main_gradient(W)
        vel = vel*mom + g*(1-mom)
        W = W - lr*vel
      }
    }else{
      stop('specify train_type as BFGS or SGD')
    }
  }
  
  
  
  #get ksm
  ##compute pdf matrix
  D = exp(dfx%*%W)
  
  ##compute means
  ksmean = numeric(k)
  for (j in 1:k){
    d = D[,j]
    ksmean[j] = sum(y*d)/sum(d)
  }
  
  
  out = list(W = W, predict = pred, ksm = ksmean)
  if (get_errors){
    out[['errors']] <- ERRORS
  }
  
  return(out) 
}







## default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}


model_msr = function(train, test, yname, 
                     error_fun=error_MAE, 
                     params=list()){
  
  N = nrow(train)
  
  if(is.null(params$niter)){
    params$niter = 50
  }
  if(is.null(params$decay)){
    params$decay = 0.001
  }
  if(is.null(params$penalty_L1)){
    params$penalty_L1 = F
  }
  if(is.null(params$penalty_local)){
    params$penalty_local = F
  }
  if(is.null(params$lr)){
    params$lr = 0.01
  }
  if(is.null(params$momentum)){
    params$momentum = 0.3
  }
  if(is.null(params$winit)){
    params$winit = 0.0001
  }  
  if(is.null(params$k)){
    params$k = 2
  }  
  if(is.null(params$verbose)){
    params$verbose = 0
  }  
  if(is.null(params$get_errors)){
    params$get_errors = T
  }  
  if(is.null(params$train_type)){
    params$train_type = 'BFGS'
  }  
  
  
  model = MSR_model(train = train, yname=yname, niter=params$niter,
                 get_errors=params$get_errors, 
                 verbose=params$verbose, decay = params$decay,
                 penalty_L1 = params$penalty_L1, 
                 penalty_local = params$penalty_local, 
                 lr=params$lr, momentum=params$momentum,
                 winit = params$winit, k = params$k, 
                 train_type=params$train_type)
  
  
  
  pred = model$predict
  
  
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






