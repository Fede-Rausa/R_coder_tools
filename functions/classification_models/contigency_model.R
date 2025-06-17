model_contigency = function(df, yname, xname){
  x = df[,xname]
  y = df[,yname]
  tab = table(x,y)
  tab
  
  xlev = rownames(tab)
  xpreds = numeric(length(xlev))
  names(xpreds) = xlev
  xprob = matrix(0, nrow=nrow(tab), ncol=ncol(tab))
  rownames(xprob) = rownames(tab)
  colnames(xprob) = colnames(tab)
  for (i in 1:length(xlev)){
    cvet = tab[xlev[i],]
    xprob[i,] = cvet/sum(cvet)
    xpreds[i] = xlev[which.max(cvet)]
  }
  
  pred = function(newdf, prob=F){
    newx = newdf[,xname]
    if (prob){
      preds = xprob[newx,]
    }else{
      preds = xpreds[newx]
    }
    return(preds)
  }
  
  acc = sum(pred(df, prob=F)==y)/nrow(df)
  
  return(list(predict = pred, 
              probmat = xprob,
              xpreds = xpreds,
              acc = acc
  ))
}





model_multicontigency = function(df, yname, xnames, w=NULL){
  
  y = df[,yname]
  ytab = table(y)
  
  if (is.null(w)){
    w = rep(1, length(xnames))
  }
  names(w) = xnames
  
  mlist = list()
  for (xn in xnames){
    mlist[[xn]] = model_contigency(df, yname, xn)
  }
  
  pred = function(newdf){
    ppred = matrix(0, nrow=nrow(newdf), ncol=length(ytab))
    colnames(ppred) = names(ytab)
    
    for (xn in xnames){
      newpred = as.matrix(mlist[[xn]]$predict(newdf, prob=T))
      newpred = newpred * w[xn]
      ppred = ppred + newpred
    }
    colnames(ppred) = names(ytab)
    
    ids = unlist(apply(ppred, 1, which.max))
    preds = colnames(ppred)[ids]
    
    return(preds)
  }
  
  return(list(predict = pred, 
              mlist = mlist,
              w=w))
  
}






