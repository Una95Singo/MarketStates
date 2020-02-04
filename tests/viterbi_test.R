# Test script for the viterbi algorithm
# una singo, SNGUNA003
# 9 January 2020

# import R functions ------------
source("R/viterbi.R")

# libraries ---------------------
library('tidyverse')
library(readxl)


# Test Case 1 -------------------
data = read_excel("data/Test/viterbiDataset.xlsx")

vit = viterbi(D=data[, 2:5], K=4, gamma=400)
vit

# Test Case 2 -------------------



# Viterbi Algorithm using true estimates -----------------------------
data = read_excel("~/Desktop/viterbiDataset.xlsx",sheet = 'Sheet2')

state1.distance = data$`True Distances`[1:10]
state2.distance = data$`True Distances`[11:20]

# Input
D = cbind(state1.distance, state2.distance)
gamma = 1
Time = nrow(D)
K = 2

# initialize
PrevCost = rep(0,K)
CurrentCost = rep(0, K)

FinalMinVal = NA
FinalPath = NA

PrevPath = list ()
PrevPath[[1]] = NA
PrevPath[[2]] = NA


CurrentPath = list()
CurrentPath[[1]] = NA
CurrentPath[[2]] = NA


#PrevCost = D[1,]
#PrevPath[[1]] = 1
#PrevPath[[2]] = 2

# looper
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

FinalPath



# Viterbi Algorithm using true estimates -----------------------------
data = read_excel("~/Desktop/viterbiDataset.xlsx",sheet = 'Sheet2')


K = 2   # number of states          
Time = length(data$Observation)
States = sample(c(1,2), size =Time, replace = T)
ExpectationMaximisation = 1
State.history = matrix(nrow = ExpectationMaximisation, ncol = Time)
Cost.history = c()

for (EM in 1:ExpectationMaximisation) {
  States = sample(c(1,2), size =Time, replace = T)
  # state index
  state1.index = which(States==1)
  state2.index = which(States==2)
  
  # state returns
  state1.return = data$Observation[state1.index]
  state2.return = data$Observation[state2.index]
  
  # State sample statistics
  state1.mu = mean(state1.return)
  state2.mu = mean(state2.return)
  state1.sd = sd(state1.return)
  state2.sd = sd(state2.return)
  
  
  # Mahalanobis distance  ---------------------------------
  state1.distance =  ((data$Observation - state1.mu)^2/state1.sd^2)
  state2.distance =  ((data$Observation - state2.mu)^2/state2.sd^2)
  
  
  # Input
  D = cbind(state1.distance, state2.distance)
  gamma = 1
  Time = nrow(D)
  K = 2
  
  
  # initialize
  PrevCost = rep(0,K)
  CurrentCost = rep(0, K)
  
  FinalMinVal = NA
  FinalPath = NA
  
  PrevPath = list ()
  PrevPath[[1]] = NA
  PrevPath[[2]] = NA
  
  
  CurrentPath = list()
  CurrentPath[[1]] = NA
  CurrentPath[[2]] = NA
  
  
  #PrevCost = D[1,]
  #PrevPath[[1]] = 1
  #PrevPath[[2]] = 2
  
  # looper
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
  
  Cost.history = c(Cost.history, min(CurrentCost))
  State.history[EM,] = unlist(FinalPath)
}


(state1.mu - state2.mu)^2
FinalPath
state1.mu
state2.mu
# get mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


hist(Cost.history)

Final  = which(Cost.history ==getmode(Cost.history))
State.history[Final,]


apply(State.history,2, getmode)

State.history[1:10,]





# Test Case 3 -------------------

