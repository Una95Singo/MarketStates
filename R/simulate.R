# simulate contains functions that simulate states and processes on states
# una singo, SNGUNA003
# 9 January 2020


#' simulate step_states produces a basic step function. 
#' @param Time length of series 
#' @param K state pattern
#' @param distr distribution of states
#' @return a series of states represented by postive integers >=1
simulate.step_states = function(Time, K, distr){
  states = numeric(Time)
  #states = states + 1
  lower_bound = 1 
  indices = cumsum(distr) * Time
  for (i in 1:length(indices)){
    upper_bound = indices[i]
    states[lower_bound:upper_bound] = states[lower_bound:upper_bound] + K[i]
    lower_bound = upper_bound + 1
  }
  return(states)
}


#' simulate markov model given a transition probability matrix
#' @param Transition transition probabolity matrix
#' @param Time length of 
#' 
simulate.markov_states = function(Transition_Matrix, Time = 50, startingState = 1 ){
  Num_States = nrow(Transition_Matrix) # number of possible states
  states = numeric(Time) 

  states[1] = startingState 
  for(t in 2:Time) {
    # probability vector to simulate next state X_{t+1}
    Conditional_Prob = Transition_Matrix[states[t-1], ]
    ## draw from multinomial and determine state
    states[t] =  which(rmultinom(1, 1, Conditional_Prob) == 1)
  }
  return(states)
}


#' simualte price paths given a state matrix
#' returns a states if no type
#' @param states
#' @param dim dimentionality of prices (number of stocks)
#' @param type  1 - gaussina noise, 2 - geometric brownian motion
#' @param mu a list of vectors
#' @param sigma a list of matricies

simulate.prices = function(states, dim = 1, type = 1, mu = 0.1, sigma = 0.01){
  #states = c(1,1,1,2,2)
  #type = 1
  #mu = 0.1
  #sigma = 0.01
  #mu =  list()
  #mu[[1]] = c(0.1)
  #mu[[2]] = c(-0.1)
  #sigma = list()
  #sigma[[1]] = 0.02
  #sigma[[2]] = 0.01
  
  #mu =  list()
  #mu[[1]] = c(0.1, 0.02)
  #mu[[2]] = c(-0.1, -0.02)
  #sigma = list()
  #sigma[[1]] =  matrix(c(0.01, 0, 0, 0.02), nrow =2)
  #sigma[[2]] =  matrix(c(0.01, 0, 0, 0.02), nrow = 2)
  
  K = length(unique(states))
  Time = length(states)
  returns = matrix(NA, nrow=dim, ncol =Time) # will need to change to matrix if dim>1
  
  if (type==1){
    
      if (K>1){
        for ( i in 1:Time){
          
          returns[,i] = mvrnorm(1, mu = mu[[states[i]]], Sigma = sigma[[states[i]]])
        }
        return(returns)
      }
    
      else{
        
        returns[,i] = mvrnorm(Time, mu = mu, Sigma = sigma)
        return(returns)
      }
  }
  if(type==2){
    return(NA)
  }
  return(states)
}



