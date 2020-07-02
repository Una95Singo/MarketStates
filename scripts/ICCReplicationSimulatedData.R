# ICC replication on simulated datasets
# Una Singo, SNGUNA003
# 5 June 2020


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/stylisedFacts.R")
source("R/NohAnzatsStateSimulation.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(markovchain)
library(igraph)

set.seed(1)



#----------------------------------------------------------------------------------------
#--------------- 2 state test------------------------------------------------------------
#----------------------------------------------------------------------------------------
# 2 state simulation
marketStates = c("1", "2")
byRow = TRUE
tMatrix = matrix(data = c(0.95, 0.05,
                           0.05, 0.95), byrow = byRow, nrow = 2, dimnames = list(marketStates, marketStates))
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = tMatrix, name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)

# Simulate dateset

twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
#plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)
twoStateData = NohAnzats.simulate( stateSequence = twoStateSequence, g = c(0.1, 0.1), clusters = 2, gaussian = T, stocks = 100)
plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, S0 =1, 'Stock Returns')

a = rowMeans(twoStateData$stockMatrix[, twoStateData$stateSeq == 1 ])
b = rowMeans(twoStateData$stockMatrix[, twoStateData$stateSeq == 2 ])

plot( x = a, y = b)
cor(a,b)

# fit ICC for iteration 1
# find states
par(mfrow = c(1,1))
seg.results = segmentation.procedure(returns =twoStateData$stockMatrix, Time = 2000, gamma=0.0015, iters =25, K = 2, sparse = 2 )
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'distance')

minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

par(mfrow = c(1,2))
plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, title = 'Ground truth market states', S0 = 100)
plotter.series(colMeans(twoStateData$stockMatrix), seg.results$history[[paste("iter:", minIndex[1],"states")]]
               , title = 'Optimal market states', S0 = 100)
ARI = adj.rand.index(twoStateData$stateSeq, seg.results$history[[paste("iter:", minIndex[1],"states")]])
ARI = round(ARI,2)
text(x = 200,y = 100*max(cumsum(colMeans(twoStateData$stockMatrix))), labels=paste('rand index:', ARI ), cex = 1)

# fit markov chain to given state sequence
stateFittedMLE = markovchainFit(data =  seg.results$history[[paste("iter:", minIndex[1],"states")]], method = "mle", name = "Two state estimate")
stateFittedMLE$estimate
ARI.history = c(ARI)
transitionEstimate.history = as(stateFittedMLE$estimate, 'matrix')

for (iter in 2:30){
  # Simulate dateset
  twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
  #plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)
  twoStateData = NohAnzats.simulate( stateSequence = twoStateSequence, g = c(0.1, 0.1), clusters = 2, gaussian = T, stocks = 100)
  #plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, S0 =1, 'Stock Returns')
  
  # find states
  par(mfrow = c(1,1))
  seg.results = segmentation.procedure(returns =twoStateData$stockMatrix, Time = 2000, gamma=0.0015, iters =15, K = 2 )
  plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'distance')
  minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))
  
  par(mfrow = c(1,2))
  plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, title = 'Ground truth market states', S0 = 100)
  plotter.series(colMeans(twoStateData$stockMatrix), seg.results$history[[paste("iter:", minIndex[1],"states")]]
                 , title = 'Optimal market states', S0 = 100)
  
  ARI = adj.rand.index(twoStateData$stateSeq, seg.results$history[[paste("iter:", minIndex[1],"states")]])
  ARI = round(ARI,2)
  text(x = 200,y = 100*max(cumsum(colMeans(twoStateData$stockMatrix))), labels=paste('rand index:', ARI ), cex = 1)
  # fit markov chain to given state sequence
  stateFittedMLE = markovchainFit(data =  seg.results$history[[paste("iter:", minIndex[1],"states")]], method = "mle", name = "Two state estimate")
  stateFittedMLE$estimate
  
  ARI.history = c(ARI.history, ARI)
  transitionEstimate.history = as(stateFittedMLE$estimate, 'matrix') + transitionEstimate.history
}
# estiamte states
transitionEstimate.history * 1/iter
mean(ARI.history)
sd(ARI.history)




#----------------------------------------------------------------------------------------
#--------------- 3 state test------------------------------------------------------------
#----------------------------------------------------------------------------------------
par(mfrow = c(1,1))
# 3 state simulation
marketStates = c("1", "2", '3')
byRow = TRUE
tMatrix = matrix(data = c(0.95, 0.025, 0.025,
                          0.025, 0.95, 0.025,
                          0.025, 0.025, 0.95), byrow = byRow, nrow = length(marketStates), dimnames = list(marketStates, marketStates))
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = tMatrix, name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)

# Simulate dateset
twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
#plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)
twoStateData = NohAnzats.simulate( stateSequence = twoStateSequence, g = c(0.1, 0.1, 0.1), clusters = 3, gaussian = T, stocks = 100)
plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, S0 =1, 'Stock Returns')

a = rowMeans(twoStateData$stockMatrix[, twoStateData$stateSeq == 1 ])
b = rowMeans(twoStateData$stockMatrix[, twoStateData$stateSeq == 2 ])

plot( x = a, y = b)
cor(a,b)

# fit ICC for iteration 1
# find states
par(mfrow = c(1,1))
seg.results = segmentation.procedure(returns =twoStateData$stockMatrix, Time = 2000, gamma=45, iters =20, K = length(marketStates) )
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'distance')

minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

par(mfrow = c(1,2))
plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, title = 'Ground truth market states', S0 = 100)
plotter.series(colMeans(twoStateData$stockMatrix), seg.results$history[[paste("iter:", minIndex[1],"states")]]
               , title = 'Optimal market states', S0 = 100)
ARI = adj.rand.index(twoStateData$stateSeq, seg.results$history[[paste("iter:", minIndex[1],"states")]])
ARI = round(ARI,2)
text(x = 200,y = 100*max(cumsum(colMeans(twoStateData$stockMatrix))), labels=paste('rand index:', ARI ), cex = 1)

# fit markov chain to given state sequence
stateFittedMLE = markovchainFit(data =  seg.results$history[[paste("iter:", minIndex[1],"states")]], method = "mle", name = "Two state estimate")
stateFittedMLE$estimate
ARI.history = c(ARI)
transitionEstimate.history = as(stateFittedMLE$estimate, 'matrix')

for (iter in 2:30){
  # Simulate dateset
  twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
  #plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)
  twoStateData = NohAnzats.simulate( stateSequence = twoStateSequence, g = c(0.1, 0.1), clusters = 2, gaussian = T, stocks = 100)
  #plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, S0 =1, 'Stock Returns')
  
  # find states
  par(mfrow = c(1,1))
  seg.results = segmentation.procedure(returns =twoStateData$stockMatrix, Time = 2000, gamma=0.0015, iters =15, K = 2 )
  plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'distance')
  minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))
  
  par(mfrow = c(1,2))
  plotter.series(colMeans(twoStateData$stockMatrix), twoStateData$stateSeq, title = 'Ground truth market states', S0 = 100)
  plotter.series(colMeans(twoStateData$stockMatrix), seg.results$history[[paste("iter:", minIndex[1],"states")]]
                 , title = 'Optimal market states', S0 = 100)
  
  ARI = adj.rand.index(twoStateData$stateSeq, seg.results$history[[paste("iter:", minIndex[1],"states")]])
  ARI = round(ARI,2)
  text(x = 200,y = 100*max(cumsum(colMeans(twoStateData$stockMatrix))), labels=paste('rand index:', ARI ), cex = 1)
  # fit markov chain to given state sequence
  stateFittedMLE = markovchainFit(data =  seg.results$history[[paste("iter:", minIndex[1],"states")]], method = "mle", name = "Two state estimate")
  stateFittedMLE$estimate
  
  ARI.history = c(ARI.history, ARI)
  transitionEstimate.history = as(stateFittedMLE$estimate, 'matrix') + transitionEstimate.history
}
# estiamte states
transitionEstimate.history * 1/iter
mean(ARI.history)
sd(ARI.history)

