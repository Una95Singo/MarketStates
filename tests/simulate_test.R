# test script for simulate 
# una singo, SNGUNA003
# 21 January 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)

# test basic step state simulations -----------------------

# test 1 - 2 states 50/50
simulated = simulate.step_states(100, c(1,2), c(0.5,0.5))
plotter.states(simulated, title = 'simulated states', thickness = 4)

# test 2 - 2 states 30/70
simulated = simulate.step_states(130, c(1,2), c(0.3,0.7))
plotter.states(simulated, title = 'simulated states', thickness = 4)

# test 3 - 2 states 30/50/20
simulated = simulate.step_states(130, c(1,2,1), c(0.3,0.5,0.2))
plotter.states(simulated, title = 'simulated states', thickness = 4)

# test 4 - 3 states 30/40/30
simulated = simulate.step_states(100, c(1,2,3), c(0.3,0.4,0.3))
plotter.states(simulated, title = 'simulated states', thickness = 4)

# test 5 - 4 states 20/40/10/30
simulated = simulate.step_states(100, c(1,2,3,4), c(0.2,0.4,0.1,0.3))
plotter.states(simulated, title = 'simulated states', thickness = 4)


# test markov step state simulations ------------------

# test 1 - 2 state markov model
Transitions = matrix(c(0.70,0.30,0.05,0.95), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 10, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  100, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)


# test 2 - 2 state markov model
Transitions = matrix(c(0.99,0.01,0.01,0.99), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 2, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  100, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)


# test 3 - 3 state markov model
Transitions = matrix(c(0.90,0.05,0.05, 0.05,0.85,0.1, 0.2,0.1,0.7), nrow = 3, ncol = 3 ,byrow = T )
Transitions
plotter.transitions(Transitions, 9, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  100, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)

# test noisy return simulations -----------------------

# test 1 single state 1 dimension
simulate.prices(c(1,1,1,1,1), mu = 0.1, sigma = 0.01)  


# test 2 two state model with one dimension
mu =  list()
mu[[1]] = c(0.9)
mu[[2]] = c(-0.1)
sigma = list()
sigma[[1]] = 0.02
sigma[[2]] = 0.01

True_States = c(1,1,1,2,2)
Sample_Returns = simulate.prices(True_States, mu = mu, sigma = sigma)  

plotter.series(colMeans(Sample_Returns),True_States, S0 =100, 'Stock Returns') 
plotter.series(colMeans(Sample_Returns),True_States, S0=NA,'Stock Returns')


# test 3 two state model with two dimensions
mu =  list()
mu[[1]] = c(0.1, 0.02)
mu[[2]] = c(-0.1, -0.02)
sigma = list()
sigma[[1]] =  matrix(c(0.01, 0, 0, 0.02), nrow =2)
sigma[[2]] =  matrix(c(0.01, 0, 0, 0.02), nrow = 2)

True_States = c(1,1,1,2,2)
Sample_Returns = simulate.prices(True_States, mu = mu, sigma = sigma, dim =2 )  

plotter.series(colMeans(Sample_Returns),True_States, S0 =100, 'Stock Returns')
plotter.series(colMeans(Sample_Returns),True_States, S0=NA, 'Stock Returns')

# test 4, random markov states plus two dimensional series
Transitions = matrix(c(0.70,0.30,0.05,0.95), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 10, 'probability transition matrix')

simulated = simulate.markov_states(Transition_Matrix = Transitions,Time =  100, startingState = 1)
plotter.states(simulated, title = 'simulated states', thickness = 4)

mu =  list()
mu[[1]] = c(0.1, 0.02)
mu[[2]] = c(-0.1, -0.02)
sigma = list()
sigma[[1]] =  matrix(c(0.01, 0, 0, 0.02), nrow =2)
sigma[[2]] =  matrix(c(0.01, 0, 0, 0.02), nrow = 2)
Returns = simulate.prices(simulated, mu = mu, sigma = sigma, dim =2 )  

plotter.series(colMeans(Returns),simulated, S0 =100, 'Stock Returns')
plotter.series(colMeans(Returns),simulated, 'Stock Returns')



