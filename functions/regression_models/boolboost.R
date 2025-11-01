
library(dplyr)


selcol = function(df, cols, rm=F){
  allcols = colnames(df)
  id = allcols %in% cols
  if (rm){
    df0 = df[,!id]
  }else{
    df0 = df[,id]
  }
  return(df0)
}

is.const = function(vec){
  return(all(vec==vec[1]))
}

cutvec = function(vec){
  val = unique(vec)
  vv = val[order(val)]
  cuts = (vv[2:length(vv)] + vv[1:(length(vv)-1)])/2
  return(cuts)
}

and = function(A, B){A*B}

or = function(A,B){A+B-and(A,B)}

not = function(A){1-A}



boolboost = function(df, yname, niter, maxnum=3, noderatio=1, xratio=1, 
                     enforceroot=F, enforcecomb=F, rooteach=0, 
                     enforceand=F, enforceor=F, andeach=0){

##maxnum: number of dummy quantile thresholds for each continuous variable
##noderatio: ratio of nodes randomly extracted
##xratio: ratio of variables randomly extracted
##rooteach: number of iterations after which a root node has to be selected

  
  dfx = df %>% dplyr::select(-yname)
  y = res = df[,yname]
  pp = data_prepro(df, yname, maxnum)
  dummydf = apply_prepro(df, pp$info)
  reports = list()
  
  
  N = nrow(dfx)
  nodedf = data.frame(x0=numeric(N))
  errors = numeric(niter)
  cmd = as.data.frame(matrix(NA, nrow=niter, ncol=10))
  
  rforce = eforce = andforce = orforce = logical(niter) #default always false
  
  
  if (enforceroot){
    rforce[1:niter] = T
  }else{
    
    if (enforcecomb){
      eforce[1:niter] = T
    }
    
    if (andeach>0){
      ids = seq(1, niter, by=andeach)
      eforce[1:niter] = T
      andforce[2] = T
      andforce[ids] = T
      orforce[-ids] = T
    }else{
      if (enforceor){
        eforce[1:niter] = T
        orforce[1:niter] = T
        andforce[1:niter] = F
      }else{
        if (enforceand){
          eforce[1:niter] = T
          andforce[1:niter] = T
          orforce[1:niter] = F
        }
      }
    }

    if (rooteach>0){
      ids = seq(1, niter, by=rooteach)
      rforce[ids] = T
      eforce[-ids] = T
    }

  }
  
  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = niter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  for (t in 1:niter){
    report = explorestep(dummydf, nodedf, res, t, noderatio, xratio, 
                         enforceroot=rforce[t], enforcecomb = eforce[t], 
                         enforceand=andforce[t], enforceor = orforce[t])
    conf = confirmstep(dummydf, nodedf, res, cmd, t, report)
    nodedf = conf$nodedf
    cmd = conf$cmd
    res = conf$res
    errors[t] = mean(res^2)
    reports[[t]] = report
    setTxtProgressBar(pb, t)
  }
  
  close(pb)
  
  predf = function(df){
    dummydf = apply_prepro(df, pp$info)
    nodedf = nodebuilder(dummydf, cmd)$nodedf
    preds = nodepredictor(nodedf, cmd)
    return(preds)
  }
  
  return(list(besterror = errors[t], errors=errors, predict=predf, cmd = cmd, reports=reports))
}


nodebuilder = function(dummydf, cmd){
  nodedf = data.frame(x0=numeric(nrow(dummydf)))
  G = nrow(cmd)
  torm = NA*numeric(G)
  cmd = as.matrix(cmd)
  for (t in 1:G){
    cmdl = cmd[t,]
    torm[t] = cmdl['oldnodename']
    nodedf = nodebuildstep(dummydf, nodedf, cmdl)
  }
  nodedf = nodedf[,-1]
  torm = na.omit(torm)
  tokeep = setdiff(colnames(nodedf), torm)
  return(list(nodedf=nodedf, torm=torm, tokeep=tokeep))
}


nodebuildstep = function(dummydf, nodedf, cmdl){
  cmdl = unlist(cmdl)
  op = cmdl['interaction']
  t = cmdl['currentid']
  xn = cmdl['newvar']
  xv = dummydf[,xn]
  
  if (op=='root'){
    nodedf = cbind(nodedf, xv)
    colnames(nodedf)[ncol(nodedf)] = paste0('V', t)
    
    #nodedf = cbind(nodedf, xv, not(xv))
    #colnames(nodedf)[ncol(nodedf)- c(1,0)] = c(paste0('V', t, 'T'), paste0('V', t, 'F'))
  }else{
    node_n = cmdl['oldnodename']
    node_v = nodedf[,node_n]
    
    if (op=='and'){
      comb = and(node_v, xv)
    }else{  
      if(op=='notand'){
        comb = and(not(node_v), xv)
      }else{
        if(op=='andnot'){
          comb = and(node_v, not(xv))
        }else{
          if(op=='notandnot'){
            comb = and(not(node_v), not(xv))
          }else{
            if (op=='or'){
              comb = or(node_v, xv)
            }else{  
              if(op=='notor'){
                comb = or(not(node_v), xv)
              }else{
                if(op=='ornot'){
                  comb = or(node_v, not(xv))
                }else{
                  if(op=='notornot'){
                    comb = or(not(node_v), not(xv))
                  }
                }
              }
            }
          }
        }
      }
    }
    nodedf = cbind(nodedf, comb)
    colnames(nodedf)[ncol(nodedf)] = paste0('V',t)
    
  }
  
  return(nodedf)
}


nodepredictor = function(nodedf, cmd){
  N = nrow(nodedf)
  G = nrow(cmd)  
  preds = numeric(N)
  cmd = as.matrix(cmd)
  for (t in 1:G){
    cmdl = cmd[t,]
    preds = preds + nodepredstep(nodedf, cmdl)
  }
  
  return(preds)
}


nodepredstep = function(nodedf, cmdl){
  cmdl = unlist(cmdl)
  node_n = cmdl['currentname']
  vb = as.logical(nodedf[,node_n])
  m1 = as.numeric(cmdl['m1'])
  m0 = as.numeric(cmdl['m0'])
  preds = ifelse(vb, m1, m0)
  return(preds)
}


explorestep = function(dummydf, nodedf=NULL, res, t, noderatio=1, xratio=1, 
                       enforceroot=T, enforcecomb=F, enforceand=F, enforceor=F){
  
  #interaction node
  #subsampling old nodes
   if (t > 3){
     g = ncol(nodedf)
     nodedf = nodedf[,sample(1:g, ceiling(g*noderatio))]
   }
  
  #subsampling original variables
  g = ncol(dummydf)
  dummydf = dummydf[,sample(1:g, ceiling(g*xratio))]
  
  
  dnames = colnames(dummydf)
  nodenames = colnames(nodedf)
  h = length(dnames)
  n = length(nodenames)
  
  report = as.data.frame(matrix(NA, ncol= 10, nrow=(h+1)*(n+1)))
  colnames(report) = c(
  'currentid', 'currentname','oldnodeid','oldnodebool',
  'oldnodename', 'interaction','newvar',
  'error', 'm0', 'm1')
  id=0
  
  
  if (t==1 || enforceroot){
    for (i in 1:h){
      
      id = id + 1
      xn = dnames[i]
      vb = as.logical(dummydf[,i])
      m0 = mean(res[!vb])
      m1 = mean(res[vb])
      res0 = res - ifelse(vb, m1, m0)
      mse = mean(res0^2)
      report[id,] = c(t, paste0('V', t), NA, NA,NA, 
                      'root', xn, mse, m0, m1)
    }
    
  }else{
    if (enforcecomb){
      if (enforceand){
        opt = 2:5
      }else{
        if (enforceor){
          opt = 6:9
        }else{
          opt = 2:9
        }
      }
    }else{
        opt = 1:9
    }
    
    for (o in opt){
      if (o==1){ #root
        for (i in 1:h){
          
          id = id + 1
          xn = dnames[i]
          vb = as.logical(dummydf[,i])
          m0 = mean(res[!vb])
          m1 = mean(res[vb])
          res0 = res - ifelse(vb, m1, m0)
          mse = mean(res0^2)
          report[id,] = c(t, paste0('V', t), NA, NA,NA, 
                          'root', xn, mse, m0, m1)
        }
        
      }else{
        if(o==2){  # and 
          
          for (i in 1:h){
            xv = dummydf[,i]
            xn = dnames[i]
            
            for (j in 1:n){
              
              #and operation
              id = id+1
              node_v = nodedf[,j]
              node_n = nodenames[j]
              comb = as.logical(and(node_v, xv))
              m0 = mean(res[!comb])
              m1 = mean(res[comb])
              res0 = res - ifelse(comb, m1, m0)
              mse = mean(res0^2)
              
              report[id,] = c(t, paste0('V', t), j, T,node_n, 
                              'and', xn, mse, m0, m1)
            }
          }        
          
          
        }else{
          if (o == 3){  #notand
            
            for (i in 1:h){
              xv = dummydf[,i]
              xn = dnames[i]
              
              for (j in 1:n){
                
                #and operation
                id = id+1
                node_v = nodedf[,j]
                node_n = nodenames[j]
                comb = as.logical(and(not(node_v), xv))
                m0 = mean(res[!comb])
                m1 = mean(res[comb])
                res0 = res - ifelse(comb, m1, m0)
                mse = mean(res0^2)
                
                report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                'notand', xn, mse, m0, m1)
              }
            }       
            
            
          }else{
            if(o==4){   #andnot
              
              for (i in 1:h){
                xv = dummydf[,i]
                xn = dnames[i]
                
                for (j in 1:n){
                  
                  #and operation
                  id = id+1
                  node_v = nodedf[,j]
                  node_n = nodenames[j]
                  comb = as.logical(and(node_v, not(xv)))
                  m0 = mean(res[!comb])
                  m1 = mean(res[comb])
                  res0 = res - ifelse(comb, m1, m0)
                  mse = mean(res0^2)
                  
                  report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                  'andnot', xn, mse, m0, m1)
                }
              }       
              
            }else{
              if(o==5){ #notandnot
                
                
                for (i in 1:h){
                  xv = dummydf[,i]
                  xn = dnames[i]
                  
                  for (j in 1:n){
                    
                    #and operation
                    id = id+1
                    node_v = nodedf[,j]
                    node_n = nodenames[j]
                    comb = as.logical(and(not(node_v), not(xv)))
                    m0 = mean(res[!comb])
                    m1 = mean(res[comb])
                    res0 = res - ifelse(comb, m1, m0)
                    mse = mean(res0^2)
                    
                    report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                    'notandnot', xn, mse, m0, m1)
                  }
                }
                
                
                
              }else{
                if(o==6){ #or
                  
                  for (i in 1:h){
                    xv = dummydf[,i]
                    xn = dnames[i]
                    
                    for (j in 1:n){
                      
                      #and operation
                      id = id+1
                      node_v = nodedf[,j]
                      node_n = nodenames[j]
                      comb = as.logical(or(node_v, xv))
                      m0 = mean(res[!comb])
                      m1 = mean(res[comb])
                      res0 = res - ifelse(comb, m1, m0)
                      mse = mean(res0^2)
                      
                      report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                      'or', xn, mse, m0, m1)
                      
                    }
                  }
                  
                  
                }else{
                  if(o==7){ #notor
                    
                    for (i in 1:h){
                      xv = dummydf[,i]
                      xn = dnames[i]
                      
                      for (j in 1:n){
                        
                        #and operation
                        id = id+1
                        node_v = nodedf[,j]
                        node_n = nodenames[j]
                        comb = as.logical(or(not(node_v), xv))
                        m0 = mean(res[!comb])
                        m1 = mean(res[comb])
                        res0 = res - ifelse(comb, m1, m0)
                        mse = mean(res0^2)
                        
                        report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                        'notor', xn, mse, m0, m1)
                        
                      }
                    }
                    
                    
                  }else{
                    if(o==8){ #ornot
                      
                      for (i in 1:h){
                        xv = dummydf[,i]
                        xn = dnames[i]
                        
                        for (j in 1:n){
                          
                          #and operation
                          id = id+1
                          node_v = nodedf[,j]
                          node_n = nodenames[j]
                          comb = as.logical(or(node_v, not(xv)))
                          m0 = mean(res[!comb])
                          m1 = mean(res[comb])
                          res0 = res - ifelse(comb, m1, m0)
                          mse = mean(res0^2)
                          
                          report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                          'ornot', xn, mse, m0, m1)
                          
                        }
                      }
                      
                      
                    }else{  #notornot
                      
                      for (i in 1:h){
                        xv = dummydf[,i]
                        xn = dnames[i]
                        
                        for (j in 1:n){
                          
                          #and operation
                          id = id+1
                          node_v = nodedf[,j]
                          node_n = nodenames[j]
                          comb = as.logical(or(not(node_v), not(xv)))
                          m0 = mean(res[!comb])
                          m1 = mean(res[comb])
                          res0 = res - ifelse(comb, m1, m0)
                          mse = mean(res0^2)
                          
                          report[id,] = c(t, paste0('V', t), j, T,node_n, 
                                          'notornot', xn, mse, m0, m1)
                          
                        }
                      }

                      
                      
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  
  return(report)
  
}


confirmstep = function(dummydf, nodedf=NULL, res, cmd, t, report){

  errors = report[,'error']
  
  bestcmd = unlist(report[which.min(errors),])
  cmd[t,] = bestcmd
  
  
  if (t==1){
    xn = bestcmd['newvar']
    nodedf = dummydf[,xn]
    nodedf = data.frame(V1 = nodedf)
    colnames(cmd) = names(bestcmd)
  }else{
    nodedf = nodebuildstep(dummydf, nodedf, bestcmd)
  }

  
  res = res - nodepredstep(nodedf, bestcmd)
  
  return(list(nodedf=nodedf, cmd=cmd, res=res))
}


data_prepro <- function(df, yname, maxnum=10) {
  dfx = df %>% dplyr::select(-yname)
  p = ncol(dfx)
  xnames = colnames(dfx)
  N = nrow(dfx)
  
  isdummy = isnumL = iscat = is1 = rep(FALSE, p)
  #newdf = data.frame(y)
  #colnames(newdf) = yname
  newdf = data.frame(x0=rep(1, N))
  
  preprocess_info = list()
  preprocess_info$variable_info = list()
  
  for (i in 1:p) {
    v = dfx[,i]
    
    name = xnames[i]
    isnum = is.numeric(v)
    values = unique(v)
    L = length(values)
    is2 = L == 2
    is1[i] = L == 1
    
    if (!isnum && !is2) {
      iscat[i] = TRUE
    }
    
    if (isnum && !is2) {
      isnumL[i] = TRUE
    }
    
    if (is2) {
      isdummy[i] = TRUE
    }
    
    var_info = list(name = name, type = NA, extra = NULL)
    
    if (isdummy[i]) {
      ref = values[order(values, decreasing=TRUE)][1]
      v0 = as.numeric(v == ref)
      var_info$type = "dummy"
      var_info$ref = ref
      newdf[[name]] = v0
    } else if (iscat[i]) {
      var_info$type = "categorical"
      var_info$levels = sort(unique(v))
      dummymat = as.data.frame(vec2onehot(v, var_info$levels))
      colnames(dummymat) = paste0(name, '_', colnames(dummymat))
      newdf = cbind(newdf, dummymat)
    } else if (isnumL[i]) {
      var_info$type = "numeric"
      cuts = cutvec(v)
      if (length(cuts) > maxnum) {
        cuts = quantile(cuts, seq(0, 1, length = maxnum + 1))[-c(1, maxnum + 1)]
      }
      var_info$cuts = cuts
      dummymat = sapply(cuts, function(s) as.numeric(v < s))
      colnames(dummymat) = paste0(name, '<', cuts)
      newdf = cbind(newdf, dummymat)
    }
    
    preprocess_info$variable_info[[name]] = var_info
  }
  
  newdf = newdf[,-1]
  return(list(data = newdf, info = preprocess_info))
}


apply_prepro <- function(newdata, preprocess_info) {
  xnames = names(preprocess_info$variable_info)
  N = nrow(newdata)
  newdf = data.frame(x0=numeric(N))
  
  for (name in xnames) {
    v = newdata[[name]]
    var_info = preprocess_info$variable_info[[name]]
    
    if (var_info$type == "dummy") {
      v0 = as.numeric(v == var_info$ref)
      newdf[[name]] = v0
    } else if (var_info$type == "categorical") {
      levels = var_info$levels
      dummymat = as.data.frame(vec2onehot(v, levels))
      colnames(dummymat) = paste0(name, '_', colnames(dummymat))
      newdf = cbind(newdf, dummymat)
    } else if (var_info$type == "numeric") {
      cuts = var_info$cuts
      dummymat = sapply(cuts, function(s) as.numeric(v < s))
      colnames(dummymat) = paste0(name, '<', cuts)
      newdf = cbind(newdf, dummymat)
    }
  }
  newdf = newdf[,-1]
  return(newdf)
}







## default loss ---

error_MAE = function(preds, y){
  return(mean(abs(preds - y)))
}



model_boolboost = function(train, test, yname, 
                     error_fun=error_MAE, params=list()){
  
  N = nrow(train)
  
  if(is.null(params$niter)){
    params$niter= 50
  }
  
  if (is.null(params$maxnum)){
    params$maxnum = 3
  }
  
  
  if (is.null(params$noderatio)){
    params$noderatio = 0.5
  }
  
  
  if (is.null(params$xratio)){
    params$xratio = 0.5
  }
  
  if (is.null(params$enforceand)){
    params$enforceand = T
  }
  
  if (is.null(params$rooteach)){
    params$rooteach = 15
  }  

  
  model = boolboost(train, yname, niter=params$niter, maxnum=params$maxnum, noderatio=params$noderatio, xratio=params$xratio, 
                    enforceand=params$enforceand,   rooteach=params$rooteach)
  
  
  
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









    


