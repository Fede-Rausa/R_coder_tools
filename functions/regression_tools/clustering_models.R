

###### unscale ------------

unscale = function(scaled_data, scales=NULL, centers=NULL, c0=F){
#scaled_data should be an output of the scale function
#with unscale, setting scales and centers, 
#you can adjust the mean and variance of your distributions

  if (is.null(centers)){
    if (c0){
      p = ncol(scaled_data)
      centers = rep(0, p)
    }else{
      centers = attr(scaled_data, 'scaled:center')
    }
  }
  if (is.null(scale)){
    scales = attr(scaled_data, 'scaled:scale')
  }
  return(scaled_data * scales + centers)
}



predict_cluster = function(df, centers){
  ##df is a dataframe or a matrix
  ##centers is a matrix with the same number of columns
  ##of df
  
  n = nrow(df)
  k = nrow(centers)
  clu = numeric(n)  
  
  for (i in 1:n){
    d = numeric(k)
    for (j in 1:k){
      d[j] = sum((df[i,] - centers[j,])^2)
    }
    clu[i] = which.min(d)
  }
  
  return(clu)
}



####### kmeans -------------------------------------
model_kmeans = function(train, test, yname=NULL,
                        error_fun=error_MAE,
                        params=list(include_y=F, 
                                    scale_w = NULL, 
                                    coordinates=NULL,
                                    median=F,
                                    k = NULL,
                                    seed=NULL)){
  
  
  if (is.null(params$median)){
    params$median = TRUE
  }
  if (is.null(params$coordinates)){
    coordinates = colnames(train)[-which(colnames(train) == yname)]
  }else{
    coordinates = params$coordinates
  }
  if (is.null(params$include_y)){
    params$include_y = F
  }
  
  if (params$include_y){
    coordinates = unique(c(coordinates, yname))
  }
  
  if (is.null(params$scale_w)){
    scale_w = rep(1, length(coordinates))
  }else{
    scale_w = params$scale_w
    if (length(scale_w) != length(coordinates)){
      stop('scale_w must have the same length as coordinates')
    }
  }
  
  
  
  scaled_df = scale(train[, coordinates])
  s_center = attr(scaled_df, "scaled:center")
  s_scale = attr(scaled_df, "scaled:scale")
  names(s_center) = names(s_scale) = names(scale_w) = coordinates

  scaled_df0 = unscale(scaled_df, scales=scale_w, c0=T)
  
  seed = params$seed
  if (is.null(seed)){
    model = kmeans(scaled_df0, centers=params$k)
  }else{
    set.seed(seed)
    model = kmeans(scaled_df0, centers=params$k)
  }
  
  ce = model$centers
  colnames(ce) = coordinates
  coors = coordinates[coordinates!=yname]
  
  
  
  pred_clu = function(df){

    coors0 = colnames(df)[colnames(df) %in% coors] 
    if (length(coors0)==0){
      stop('No coordinates found in the dataframe')
    }
    
    dfmat = as.matrix(unscale(scale(df[, coors0], 
                           center=s_center[coors0], 
                           scale=s_scale[coors0]), 
                           scales=scale_w[coors0], c0=T))
    
    ce0 = ce[,coors0]
    
    clu = predict_cluster(dfmat, ce0)
    
    return(clu)
  }
  

  build_y_dict = function(clusters, yvec, median=T){
    
    k = max(clusters)
    y_dict = numeric(k)
    names(y_dict) = 1:k
    if (is.null(yvec)){
      stop('yvec cannot be NULL')
    }
    
    if (is.numeric(yvec)){
      if (!median){
        for (i in 1:k){
          y_dict[i] = mean(yvec[clusters == i], na.rm=T)
        }
      }else{
        for (i in 1:k){
          y_dict[i] = median(yvec[clusters == i], na.rm=T)
        }
      }
      
    }else{
      for (i in 1:k){
        tab = table(yvec[clusters == i])
        y_dict[i] = names(tab)[which.max(tab)]
      }
    }
    
    return(y_dict)
  }
  
  
  pred_y = function(df, y_dict0=NULL){
    if (is.null(y_dict0)){
      y_dict0 = y_dict
    }
    clu = pred_clu(df)
    preds = y_dict[clu]
    return(preds)
  }
  
  
  train_clu = pred_clu(train)
  test_clu = pred_clu(test)
  
  
  if (!is.null(yname)){
    
    #print('ok')
    y_dict = build_y_dict(model$cluster,  #pred_clu(train),
                          train[, yname], 
                          median=params$median)
    
    train_preds = pred_y(train, y_dict)
    test_preds = pred_y(test, y_dict)
    train_error = error_fun(train_preds, train[, yname])
    test_error = error_fun(test_preds, test[, yname])
    train_residuals = train_preds - train[, yname]
    test_residuals = test_preds - test[, yname]
  }else{
    y_dict = NULL
    train_preds = NULL
    test_preds = NULL
    train_error = NULL
    test_error = NULL
    train_residuals = NULL
    test_residuals = NULL
  }
  
  
  return(list(model=model, 
              predict=pred_y, 
              pred_clu=pred_clu, 
              build_y_dict = build_y_dict,
              centers=model$centers,
              train_clu = train_clu,
              test_clu = test_clu,
              train_residuals=train_residuals,
              test_residuals=test_residuals,
              train_error=train_error,
              test_error=test_error,
              train_preds=train_preds,
              test_preds=test_preds,
              y_dict=y_dict
              ))
}




######## pam ---------------------------------------

library(cluster)

model_pam = function(train, test, yname=NULL,
                        error_fun=error_MAE, 
                        params=list(include_y=F,
                                    scale_w = NULL,
                                    coordinates=NULL,
                                    k=NULL,
                                    median=F,
                                    seed=NULL)){
  
  if (is.null(params$median)){
    params$median = TRUE
  }
  if (is.null(params$coordinates)){
    coordinates = colnames(train)[-which(colnames(train) == yname)]
  }else{
    coordinates = params$coordinates
  }
  if (is.null(params$include_y)){
    params$include_y = F
  }
  
  if (params$include_y){
    coordinates = unique(c(coordinates, yname))
  }
  
  if (!is.null(params$scale_w)){
    scale_w = params$scale_w
  }
  
  if (is.null(params$scale_w)){
    scale_w = rep(1, length(coordinates))
  }else{
    if (length(scale_w) != length(coordinates)){
      stop('scale_w must have the same length as coordinates')
    }
  }

  scaled_df = scale(train[, coordinates])
  s_center = attr(scaled_df, "scaled:center")
  s_scale = attr(scaled_df, "scaled:scale")
  names(s_center) = names(s_scale) = names(scale_w) = coordinates
  coors = coordinates[coordinates!=yname]
  
  
  scaled_df0 = unscale(scaled_df, scales=scale_w, c0=T)
  
  
  seed = params$seed
  if (is.null(seed)){
    model = pam(scaled_df0, k=params$k)
  }else{
    set.seed(seed)
    model = pam(scaled_df0, k=params$k)
  }
  
  
  ce = model$medoids
  colnames(ce) = coordinates
  
  
  pred_clu = function(df){
    
    coors0 = colnames(df)[colnames(df) %in% coors] 
    if (length(coors0)==0){
      stop('No coordinates found in the dataframe')
    }
    
    dfmat = as.matrix(unscale(scale(df[, coors0], 
                                    center=s_center[coors0], 
                                    scale=s_scale[coors0]), 
                              scales=scale_w[coors0], c0=T))
    
    ce0 = ce[,coors0]
    
    clu = predict_cluster(dfmat, ce0)
    
    return(clu)
  }
  

  
  build_y_dict = function(clusters, yvec, median=T){
    
    k = max(clusters)
    y_dict = numeric(k)
    names(y_dict) = 1:k
    if (is.null(yvec)){
      stop('yvec cannot be NULL')
    }
    
    if (is.numeric(yvec)){
      if (!median){
        for (i in 1:k){
          y_dict[i] = mean(yvec[clusters == i], na.rm=T)
        }
      }else{
        for (i in 1:k){
          y_dict[i] = median(yvec[clusters == i], na.rm=T)
        }
      }
      
    }else{
      for (i in 1:k){
        tab = table(yvec[clusters == i])
        y_dict[i] = names(tab)[which.max(tab)]
      }
    }
    
    return(y_dict)
  }
  
  
  pred_y = function(df, y_dict0=NULL){
    if (is.null(y_dict0)){
      y_dict0 = y_dict
    }
    if (yname %in% colnames(df)){
      id_y = which(colnames(df) == yname)
      clu = pred_clu(df[,-id_y])
    }else{
      clu = pred_clu(df)
    }
    
    preds = y_dict0[clu]
    return(preds)
  }
  
  
  train_clu = pred_clu(train)
  test_clu = pred_clu(test)
  
  
  if (!is.null(yname)){
    
    #print('ok')
    y_dict = build_y_dict(model$clustering,  #pred_clu(train),
                          train[, yname], 
                          median=params$median)
  
    
    train_preds = pred_y(train, y_dict)
    test_preds = pred_y(test, y_dict)
    
    
    train_error = error_fun(train_preds, train[, yname])
    test_error = error_fun(test_preds, test[, yname])
    train_residuals = train_preds - train[, yname]
    test_residuals = test_preds - test[, yname]
  }else{
    y_dict = NULL
    train_preds = NULL
    test_preds = NULL
    train_error = NULL
    test_error = NULL
    train_residuals = NULL
    test_residuals = NULL
  }
  
  
  return(list(model=model, 
              predict=pred_y, 
              pred_clu=pred_clu, 
              build_y_dict = build_y_dict,
              centers=model$medoids,
              train_clu = train_clu,
              test_clu = test_clu,
              train_residuals=train_residuals,
              test_residuals=test_residuals,
              train_error=train_error,
              test_error=test_error,
              train_preds=train_preds,
              test_preds=test_preds,
              y_dict=y_dict
  ))
}
