Simulated_Annealing = function(state, fun_energy, fun_sample, T0=100,
                      alpha=0.99, tmax=10,  max_iter=100, belisle_cooling = T){
  
  ##alpha: when 0 < alpha < 1, alpha is the ratio
  ##at which the temperature decreases over time
  ##alpha is used only if belisle_cooling is FALSE
  ##if belisle is TRUE, then instead of alpha we use tmax
  ##that express the iteration at which the temperature falls
  
  t = T0
  state_list = list()
  energy_progress = numeric(max_iter)
  t_progress = numeric(max_iter)
  best_state <- state
  best_energy <- fun_energy(state)
  
  if (belisle_cooling) {
    update_t = function(t){t / log(((i-1) %/% tmax)*tmax + exp(1))}
  }else{
    update_t = function(t){return(t*alpha)}
  }
  
  
  for (i in 1:max_iter) {
    new_state <- fun_sample(state)
    new_energy <- fun_energy(new_state)
    
    energy_progress[i] <- new_energy
    t_progress[i] = t
    
    if (new_energy < best_energy){
      state <- best_state <- new_state
      best_energy <- new_energy
    }else{
      if (runif(1) < exp((best_energy - new_energy) / t)){
        state <- new_state
      }
    }
    
    state_list[[i]] <- state
    t <- update_t(t)
  }

  return(list(best_state=best_state, 
              best_energy = best_energy, 
              energy_progress = energy_progress,
              temperature = t_progress,
              state_list = state_list))
}

