# Viterbi algorithm
# una singo, SNGUNA003
# 9 January 2020
# viterbi function finds the optimal path given a matrix of K states over T time points


#' @param D matrix of distances
#' @param K number of states
#' @param gamma temporal consistency parameter
viterbi <- function(D, K, gamma) {
  Time = nrow(D)
  
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



# --- function test
#data = read_excel("data/Test/viterbiDataset.xlsx")
#D = data[, 2:5]
#K = 4 
#gamma = 400 



#vit = viterbi(D=data[, 2:5], K=4, gamma=1)
#vit

