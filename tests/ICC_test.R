# segmentation test
# una singo, SNGUNA003
# 29 January 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/ICC.R")
source("R/simulate.R")

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
par(mfrow = c(2,2))
seg.results = ICC.cluster(returnsMatrix = returns, gamma = 0, sparseMethod = 1, distanceFunction = 1, K =2)

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$OptimalPath, title = 'Optimal market states')

ARI = adj.rand.index(simulated,seg.results$OptimalPath)
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
seg.results = ICC.cluster(returnsMatrix = returns, gamma = 0, sparseMethod = 1, distanceFunction = 1, K =2)

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$OptimalPath, title = 'Optimal market states')

ARI = adj.rand.index(simulated,seg.results$OptimalPath)
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
seg.results = ICC.cluster(returnsMatrix = returns, gamma = 0, sparseMethod = 1, distanceFunction = 1, K =3)
par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states')
plotter.series(colMeans(returns),seg.results$OptimalPath
               , title = 'Optimal market states')


ARI = adj.rand.index(simulated,seg.results$OptimalPath)
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
par(mfrow = c(1,1))
seg.results = ICC.cluster(returnsMatrix = returns, gamma = 0, sparseMethod = 1, distanceFunction = 1, K =2)

par(mfrow = c(1,2))
plotter.series(colMeans(returns),simulated,title = 'Ground truth market states', S0 = 100)
plotter.series(colMeans(returns),seg.results$OptimalPath
               , title = 'Optimal market states', S0 = 100)

ARI = adj.rand.index(simulated,seg.results$OptimalPath)

ARI = round(ARI,2)
text(x=60,y=0, labels=paste('rand index:', ARI ), cex = 1)




