# segmentation test
# una singo, SNGUNA003
# 29 January 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/ICC_EM.R")


# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)

set.seed(1)
# ----------------------------------------------------------
# test 1 
# two state two stock market with zero noise (step function)
# ----------------------------------------------------------
par(mfrow = c(1,1))
simulated = simulate.step_states(130, c(1,2,1), c(0.3,0.5,0.2))
plotter.states(simulated, title = 'simulated states', thickness = 4)
returns =  rbind(simulated*0.2, simulated*0.3)
plotter.series(colMeans(returns),simulated,title = 'two stock - two market state')
# find states



# ICC setup
gamma = 0
sparseMethod = 1
distanceFunction = 1 
K =2

ICC.Output = ICC.cluster(returnsMatrix = returns, sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 100)
seg.results = ICC.Output$OptimalPath
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')
# find the minimum
minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))
minIndex

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$history$`iter: 2 states`, title = 'Optimal market states')

ARI = adj.rand.index(simulated,seg.results$history[[paste("iter:", minIndex[1],"states")]])
ARI = round(ARI,2)
text(x=60,y=0.45, labels=paste('rand index:', ARI ), cex = 1)



# ----------------------------------------------------------
# test 2 markov switching states - two state markov model
# ----------------------------------------------------------

par(mfrow = c(1,1))
# setup transiton matrix for markov states.
Transitions = matrix(c(0.70,0.30,0.10,0.90), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 10, 'probability transition matrix')
simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  130, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)

par(mfrow = c(1,1))
returns =  rbind(simulated*0.2, simulated*0.3)
plotter.series(colMeans(returns),simulated,title = 'two stock - two markov market state')

# find states
par(mfrow = c(1,1))
seg.results = segmentation.procedure(returns = returns, Time = 130, gamma=0.0001, iters =50, K = 2 )
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')

minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$history[[paste("iter:", minIndex[1],"states")]]
, title = 'Optimal market states')


ARI = adj.rand.index(simulated,seg.results$history[[paste("iter:", minIndex[1],"states")]])
ARI = round(ARI,2)

text(x=60,y=0.45, labels=paste('rand index:', ARI ), cex = 1)




# ----------------------------------------------------------
# test 3 markov switching states - three state markov model
# ----------------------------------------------------------
Transitions = matrix(c(0.90,0.05,0.05, 0.05,0.85,0.1, 0.2,0.1,0.7), nrow = 3, ncol = 3 ,byrow = T )
Transitions
plotter.transitions(Transitions, 9, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  130, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)

par(mfrow= c(1,1))
returns =  rbind(simulated*0.2, simulated*0.3)
plotter.series(colMeans(returns),simulated,title = 'two stock - three state market')

# find states
par(mfrow = c(1,1))
seg.results = segmentation.procedure(returns = returns, Time = 130, gamma=0.001, iters =50, K = 3 )
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')

minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$history[[paste("iter:", minIndex[1],"states")]]
               , title = 'Optimal market states')


ARI = adj.rand.index(simulated,seg.results$history[[paste("iter:", minIndex[1],"states")]])
ARI = round(ARI,2)
ARI
text(x=70,y=0.6, labels=paste('ARI:', ARI ), cex = 0.75)




# ----------------------------------------------------------
# test 4 gaussian hidden markov model - two state markov model.
# ----------------------------------------------------------


Transitions = matrix(c(0.70,0.30,0.10,0.90), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 10, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  100, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)

mu =  list()
mu[[1]] = c(0.1, 0.02)
mu[[2]] = c(-0.2, -0.02)
sigma = list()
sigma[[1]] =  matrix(c(0.01, 0, 0, 0.02), nrow =2)
sigma[[2]] =  matrix(c(0.01, 0, 0, 0.02), nrow = 2)
returns = simulate.prices(simulated, mu = mu, sigma = sigma, dim =2 )  

plotter.series(colMeans(returns),simulated, S0 =100, 'Stock Returns')
plotter.series(colMeans(returns),simulated, S0=NA, 'Stock Returns')

# find states

source('R/NohAnzatsStateSimulation.R')


x = NohAnzats.simulate(stateSequence = simulated)


par(mfrow = c(1,1))


# ICC setup
gamma = 0
sparseMethod = 1
distanceFunction = 1 
K =2

ICC.Output = ICC.cluster(returnsMatrix = returns, sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 100)
seg.results = ICC.Output$OptimalPath
seg.results = segmentation.procedure(returns =returns, Time = 100, gamma=0.001, iters = 100, K = 2 , sparse = 1)
plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')

minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states', S0 = 100)
plotter.series(colMeans(returns),seg.results$history[[paste("iter:", minIndex[1],"states")]]
               , title = 'Optimal market states', S0 = 100)

ARI = adj.rand.index(simulated,seg.results$history[[paste("iter:", minIndex[1],"states")]])

ARI = round(ARI,2)
text(x=60,y=0, labels=paste('rand index:', ARI ), cex = 1)




