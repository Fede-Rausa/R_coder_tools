



PSO = function(init_vec, fun_loss, init_bounds=c(-100,100),
               np=5, lr=0.1, momentum = 0.5, gw = 0.5, niter=100){
#given a function fun_loss that receives as argument just one vector of the length of init_vec
#of real numbers, it searches for the minimum

#niter is the number of iterations
#np is the number of particles
#lr is the learning rate
#momentum is the ratio of the past velocity against the new update
#gw is the ratio of the importance of the global optimum against the particle specific optimum

  
  gopt = init_vec
  gloss = fun_loss(gopt)
  
  best_losses = numeric(niter-1)
  
  p_list = list()
  p_opt = list()
  l_list = list()
  v_list = list()
  for (j in 1:np){
    p_list[[j]] = p_opt[[j]] = init_vec + 
      runif(length(init_vec), 
             min=init_bounds[1], max=init_bounds[2])
    l_list[[j]] = c(gloss, numeric(niter-1))
    v_list[[j]] = p_list[[j]] - gopt
  }
  
  for (t in 1:(niter-1)){
    
    for (j in 1:np){
      
      p0 = p_list[[j]]
      lopt = p_opt[[j]]
      v0 = v_list[[j]]
      w_opt = gopt*gw + lopt*(1-gw)
      v1 = momentum*v0 + (1-momentum)*(p0 - w_opt)
      p1 = p0 - lr * v1
      p_list[[j]] = p1
      
      l1 = fun_loss(p1)
      
      if (l1 < l_list[[j]][t]){
        p_opt[[j]] = p1
        
        if (l1 < gloss){
          gopt = p1
          gloss = l1
        }
      }
      
      l_list[[j]][t+1] = l1
      
    }
    best_losses[t] = gloss
  }
  
  return(list(gopt=gopt, 
              gloss = gloss,
              p_list=p_list, 
              p_opt=p_opt, 
              l_list=l_list, 
              v_list=v_list, 
              best_losses= best_losses))
}





