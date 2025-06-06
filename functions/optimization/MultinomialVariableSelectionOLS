


MVS_OLS = function(train_df, test_df, yname, 
                  nsim=100, lr=0.01, lambda=NULL, 
                  probplot=F, rankplot=F, initp = 0.5){
    #MVS: multinomial variable selection
  
  p = ncol(train_df) - 1
  
  if (is.null(lambda)){
    lambda = round(sqrt(p/2))
  }
  
  target = train_df[,yname]
  train_df = train_df %>% dplyr::select(-all_of(yname))
  train_df[,yname] = target
  
  
  pvec = rep(initp, p)
  error = ns =  numeric(nsim)
  matprob = matsel = deltas = matrix(0, ncol=p, nrow=nsim)
  
  sigmoid = function(x){exp(x)/(1+exp(x))}
  logit = function(x){log(x/(1-x))}
  
  for (t in 1:nsim){
    
    #sampling ids
    n = max(rpois(1, lambda), 1)
    
    sel = c(rmultinom(1, n, pvec))
    sel = (sel>0)
    if (!any(sel)){
      sel = (c(rmultinom(1, n, pvec)) > 0)
      if (!any(sel)){
        sel = (c(rmultinom(1, n, pvec)) > 0)
      }
    }
    
    f = formula(paste0(yname, '~.'))
    model = lm(f, data=train_df[,c(sel, TRUE)])
    
    # save results
    matsel[t,] = sel
    error[t] = mean(abs(predict(model, test_df) - test_df[,yname]))
    matprob[t,] = pvec
    ns[t] = n
    
    if (t > 1){
      #update
      
      gain = (error[t] < error[t-1])*1
      loss = (error[t] > error[t-1])*1
      add = (matsel[t,] & !matsel[t-1,])*1
      remove = (!matsel[t,] & matsel[t-1,])*1
      
      delta = abs((error[t-1]-error[t])/error[t-1])
      delta = (delta/n) * (gain*add + loss*remove - gain*remove - loss*add)
      #pvec[sel] = sigmoid(logit(pvec[sel]) + lr*delta)
      pvec= sigmoid(logit(pvec) + lr*delta)
      
      deltas[t,] = delta
    }
  }
  
  if (probplot){
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    
    
    colnames(matprob) = colnames(train)[1:p]
    df_prob = as.data.frame(matprob) %>% 
      mutate(iteration = 1:nsim) %>%
      pivot_longer(cols = -iteration, 
                   names_to = 'variable', 
                   values_to = 'multiprob')
    
    
    prplot = ggplot(df_prob, aes(x = iteration, y = multiprob, color = variable, type=variable)) +
      geom_line(linewidth = 1) +
      labs(title = "Multinomial probabilities for OLS variable selection over time",
           x = "Time",
           y = "Value",
           color = "Variable") + # Set the legend title
      theme_minimal() # Use a minimal theme for cleaner look (optional)
  }
  
  
  if (rankplot){
    
    ranks = rank(pvec)
    errn = numeric(p)
    for (i in 1:p){
      sel = (ranks > (p-i))
      model = lm(y~., data=train_df[,c(sel, TRUE)])
      errn[i] = mean(abs(predict(model, test_df) - test_df[,yname]))
    }
  
    
    
    errn_df = data.frame(MAE=errn, rank=1:p)
    rkeplot = elbowK(errn, ylab='MAE', maxK = p,
                     main = 'Elbow for optimal number of variables')
    
    rkpplot = elbowK(pvec[order(pvec, decreasing=T)], 
                     ylab='multinomial probs', maxK = p,
                     main = 'Elbow for optimal number of variables')
  }else{
    rkeplot = NULL
    rkpplot = NULL
  }
  
  
  
  names(pvec) = colnames(matprob) = colnames(deltas) = colnames(train_df)[1:p]
  return(list(p = pvec[order(pvec, decreasing=T)], 
              ranks = rank(pvec), lr=lr,
              ranksMean = rank(colMeans(matprob)) , 
              MAE = error, ns=ns, probplot=prplot,
              rankplot_error = rkeplot,
              rankplot_prob = rkpplot,
              deltas=deltas, matprob=matprob))
}



elbowK <- function(vet,
                   ylab = 'metric',
                   maxK = NULL,
                   main = 'Elbow method to find optimal k') {
  
  if (is.null(maxK)) {
    maxK <- length(vet)
  }
  Kopt <- 1:maxK
  data <- data.frame(K = Kopt, Metric = vet[Kopt])
  
  a <- data[1, ]
  b <- data[nrow(data), ]
  
  m <- -(abs(b$Metric - a$Metric)) / abs((b$K - a$K))
  q <- a$Metric - m * a$K
  mp <- -1 / m
  
  distances <- numeric(nrow(data))
  projection_points <- data.frame(K_proj = numeric(nrow(data)), Metric_proj = numeric(nrow(data)))
  
  for (i in 1:nrow(data)) {
    c_point <- data[i, ]
    qp <- c_point$Metric - mp * c_point$K
    K_proj <- (qp - q) / (m - mp)
    Metric_proj <- m * K_proj + q
    projection_points[i, ] <- c(K_proj, Metric_proj)
    distances[i] <- sqrt((c_point$K - K_proj)^2 + (c_point$Metric - Metric_proj)^2)
  }
  
  optimal_point <- data[which.max(distances), ]
  optimal_projection <- projection_points[which.max(distances), ]
  bestK <- optimal_point$K
  
  p <- ggplot(data, aes(x = K, y = Metric)) +
    geom_line(color = 'blue') +
    geom_point(color = 'blue') +
    annotate("segment", x = a$K, y = a$Metric, xend = b$K, yend = b$Metric,
             linetype = 'dashed', color = 'red') +
    annotate("segment", x = optimal_point$K, y = optimal_point$Metric,
             xend = optimal_projection$K_proj, yend = optimal_projection$Metric_proj,
             linetype = 'dashed', color = 'red') +
    geom_point(data = optimal_point, aes(x = K, y = Metric),
               shape = 8, color = 'red', size = 3) +
    ggtitle(main) + ylab(ylab)+
    theme_bw()
  
  if (!is.null(maxK) && maxK > 0) {
    p <- p + xlim(0, maxK)
  }
  
  text_x <- quantile(data$K, 0.75)
  text_y <- max(data$Metric)
  p <- p + annotate("text", x = text_x, y = text_y,
                    label = paste('good k =', as.character(bestK)),
                    hjust = 0, vjust = 1)
  
  
  return(list(k=bestK, p=p, d=distances))
}
