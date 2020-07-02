# Viterbi algorithm
# una singo, SNGUNA003
# 9 January 2020
# viterbi function finds the optimal path given a distance matrix of K states over T time points

viterbi <- function(D,      # matrix of distances for each ith observation (T x K)
                    K,      # Number of states to optimise
                    gamma,  # State switching penalty
                    debug = F # activate tracing statements for debugging purposes
                    ){
  
  
  Time = nrow(D)
  if(debug){
    trace(print)
  }
  
  # initialize
  PrevCost = rep(0,K)
  CurrentCost = rep(0, K)
  
  FinalMinVal = NA
  FinalPath = NA
  
  PrevPath = list ()
  CurrentPath = list()
  
  for (i in 1:K){
    PrevPath[[i]] = NA
    CurrentPath[[i]] = NA
  }

  # Viterbi loop
  for ( t in 1:Time){
    for (k in 1:K){
      MinVal = which(PrevCost == min(PrevCost))[1]
      if((PrevCost[MinVal]+gamma) > PrevCost[k]){
        CurrentCost[k] = PrevCost[k] + D[t,k]
        CurrentCost = unlist(CurrentCost)
        CurrentPath[[k]] = c(PrevPath[[k]], k) 
      }
      else{
        CurrentCost[k] = PrevCost[MinVal] + gamma  + D[t,k]
        CurrentCost = unlist(CurrentCost)
        CurrentPath[[k]] = c(PrevPath[[MinVal]], k) 
      }
      #PrevCost = CurrentCost
      #PrevPath = CurrentPath
      
      if(debug){
        print(PrevCost)
        print(PrevPath)
      }
    }
    #update after iterating throuh the states. the paper updates before the state ends.
    PrevCost = CurrentCost
    PrevPath = CurrentPath
  }
  
  
  CurrentCost
  FinalMinVal = which(CurrentCost == min(CurrentCost))[1]
  FinalPath = CurrentPath[[FinalMinVal]]
  FinalPath = na.omit(FinalPath)
  
  # return the final path and cost pair of the viterbi algorithm
  return(list('Final_Path' = FinalPath, 'Final_Cost'= min(CurrentCost)))
}
