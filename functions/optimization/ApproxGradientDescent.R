
######## function for approximated gradient descent --------------



agd = function(init, fun, niter=20, lr=0.1, h=1e-6, 
               init2 = NULL, gnoise=TRUE, eps=NULL, 
               hconstant = T, lrspeed=0, momentum=0,
               restart = F, lrcooling=T,
               type='spsa'){
  
  ##global variables
  p = length(init)
  xmat = gmat = vmat = matrix(0, ncol=p, nrow=niter)
  y = numeric(niter)
  f = fun
  lr0 = lr
  
  if (is.null(eps)){
    eps=h
  }
  h = rep(h,p)
  
  ##second input
  xmat[1,] = init
  if (is.null(init2)){
    if (gnoise){
      xmat[2,] = xmat[1,] + rnorm(p)
    }else{
      xmat[2,] = xmat[1,] + ifelse(rbinom(p, 1, 0.5), -1, 1)
    }
  }else{
    xmat[2,] = init2
  }
  
  y[1] = f(xmat[1,])
  y[2] = f(xmat[2,])
  
  
  #################gradient functions
  
  ##good gradient functions
  if (type=='coordinatewise'){
    grad <- function(t) {
      x = xmat[t,]
      g <- numeric(p)

      for (i in 1:p) {
        x_forward <- x
        x_forward[i] <- x_forward[i] + h[i]
        g[i] <- (f(x_forward) - y[t]) / h[i]
      }
      
      return(g)
    }
  }
  
  
  if (type == 'random') {
    grad = function(t) {
      u = rnorm(p) 
      u = u / sqrt(sum(u^2))
      g = (f(xmat[t,] + h*u) - y[t]) / h * u
      return(g)
    }
  }
  
  if (type == 'fdsa') {
    grad = function(t) {
      delta = ifelse(rbinom(p, 1, 0.5), 1, -1)
      f_plus  = f(xmat[t,] + h*delta)
      f_minus = f(xmat[t,] - h*delta)
      g = (f_plus - f_minus)/(2*h) * delta
      return(g)
    }
  }
  
  ## bad gradient functions
  if (type=='global'){
    grad = function(t){
      g = (f(xmat[t,]+h) - y[t])/h
      return(g)
    }
  }

  if (type=='classic'){
    grad = function(t){
      g = (f(xmat[t,]+h) - f(xmat[t,]-h))/(2*h)
      return(g)
    }
  }

  if (type=='mnemonic'){
    grad = function(t){
      dx = (xmat[t,] - xmat[t-1,])
      g = (yv[t] - yv[t-1]) * dx / sum(dx^2)
      return(g)
    }
  }
  
  

  
  ###############################loop
  for (t in 2:(niter-1)) {
    
    if (!hconstant){
      h = eps * pmax(abs(xmat[t,]), 1)
    }
    
    #lr cooling
    if (lrcooling){
      lr = lr0/(1+lrspeed*t)
    }
    
    #lr restart
    if (restart){
      if (y[t] > y[t-1]) {
        lr <- lr * 0.5
      }
    }

    gmat[t,] = grad(t)
    vmat[t,] = momentum * vmat[t-1,] + (1 - momentum) * gmat[t,]
    xmat[t+1,] <- xmat[t,] - lr * vmat[t,]
    y[t+1] <- f(xmat[t+1,])
  }
  
  return(list(besty=y[which.min(y)], 
              bestx=xmat[which.min(y),],
              y=y, x=xmat, g=gmat))
}
  
  
####example usage ----

# fun <- function(x) { x[1]^2 + x[2]^2}
  
  
# opt = agd(rnorm(2)*10000, niter=30, 
#          type='fdsa', lr=0.2,  gnoise=T,
#          eps=0.1, momentum=0, 
#          hconstant = T,
#          fun=fun, h=0.1)

