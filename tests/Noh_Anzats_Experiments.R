# Noh anzats experimentß test
# una singo, SNGUNA003
# 26 April 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/stylisedFacts.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(timeSeries)
library(ape)
library(Rcpp)
library(umap)
library(emstreeR)
library(randomcoloR)

set.seed(1)


# ----------------------------------------------------------
# Noh anzats horizontal (stock based) clustering simulation
# ----------------------------------------------------------


# step 1 Cluster numbers --------------------------------------------------------
C = 4  # three clusters of stocks
s = 20 # clusters size is 50 each
N = C*s # number of stocks
D = 365  # days
# step 2 Spin labels - generate spin labels using custom similate R script --------------------------------------------------------

spinLabels= simulate.step_states(N, c(1,2,3,4), c(0.25,0.25,0.25,0.25))
par(mfrow= c(1,1))
plotter.states(spinLabels, title = 'Spin Labels viz.', thickness = 4)

# step 3 Random effects --------------------------------------------------------

eta = matrix( rnorm (C*D, 0,1), C, D) # daily effects
eps = matrix( rnorm (C*D, 0,1), N, D) # specific random effects

# step 4 fix Intra-cluster binding strength --------------------------------------------------------

gs = replicate(C,0.2) # same cluster binding strength

# step 5 compute daily returns --------------------------------------------------------

xi = matrix(NA, N, D) # matrix to store daily returns


# double loop to compute daily returns

# double for loop for simualation
for (n in 1:N){
  for (day in 1:D){
    cluster = spinLabels[n]
    xi[n, day] = (sqrt(gs[cluster]) * eta[cluster, day] + eps[n, day])/(sqrt(1+gs[cluster]))
  }
}

# step 6 visualisations ---------------------------------------------------------------

# timeseries plots. 
xi.ts = ts(xi, start =1, end = 365)
ts.plot(xi.ts, plot.type ='single', col = spinLabels)

xi.cummulative = ts(colCumsums(xi.ts))
ts.plot(xi.cummulative, plot.type ='single', col = spinLabels, main ='Cummulative Sum')


plot(cumsum(colMeans(xi)), type='l', main='market return')


# ----------------------------------------------------------
# Noh anzats vertical (day based) clustering simulation 
# ----------------------------------------------------------



# step 1 Cluster numbers --------------------------------------------------------
C = 2  # two clusters of stocks
s = 1000 # clusters size is 50 each
N = C*s # number of days
D = 400  # stocks

# step 2 Spin labels - generate spin labels (daily state assignment) --------------------------------------------------------

# deterministic class labels
spinLabels= simulate.step_states(N, c(1,2), c(0.3,0.7))
par(mfrow= c(1,1))
plotter.states(spinLabels, title = 'Spin Labels viz.', thickness = 4)

# stochastic class labels
# test 1 - 2 state markov model
par(mfrow =  c(1,1))
Transitions = matrix(c(0.99,0.01,0.01,0.99), nrow = 2, ncol = 2 ,byrow = T )
plotter.transitions(Transitions, 10, 'probability transition matrix')
par(mfrow =  c(1,1))
spinLabels = simulate.markov_states(Transition_Matrix = Transitions,Time =  N, startingState = 1)
plotter.states(spinLabels, title = 'simulated states', thickness = 4)


# step 3 Random effects --------------------------------------------------------

etaOne = rnorm(D, 0, 1)
etaTwo = rnorm(D, 0, 1)

#etaOne = rt(D, 3)
#etaTwo = rt(D, 2)
par(mfrow =  c(1,2))
hist(etaOne, main ='State 1 - t-dist 3 df')
hist(etaTwo, main ='State 2 - t-dist 2 df')

eta = matrix( c(etaOne, etaTwo), C, D, byrow = T) # daily effects
#eta = matrix( c(rt(D,3), rt(D,3)) , C, D, byrow = T) # daily effects

eps = matrix(NA, N, D) # specific random effects
# expsilon
for (n in 1:N) {
  eps[n, ] = rnorm(D, 0,1)
}

# test market return
par(mfrow =  c(1,1))
plot(cumsum(colMeans(t(eps))), type='l', main ='Cummulative plot of epsilon', ylab= 'price')


# step 4 fix Intra-cluster binding strength --------------------------------------------------------
gs = c(0.1,0.1) # same cluster binding strength

# step 5 compute daily returns --------------------------------------------------------
xi = matrix(NA, N, D) # matrix to store daily returns

# double for loop for simualation
for (n in 1:N){
  for (d in 1:D){
    cluster = spinLabels[n]
    xi[n, d] = (sqrt(gs[cluster]) * eta[cluster, d] + eps[n, d]) / (sqrt(1+gs[cluster]))
  }
}

# step 6 visualisations ---------------------------------------------------------------

simMarketStocks = t(xi) 
simDailyAvg = colMeans(simMarketStocks)
Time = 1:length(simDailyAvg)

plotter.series(simDailyAvg, spinLabels, S0 =100, 'Stock Returns')
#plotter.series(c(1,simDailyAvg), spinLabels,'Stock Returns')

stylisedFacts.plot(simDailyAvg)

# visualise a sample of 100 stocks
par(mfrow = c(1,1))
x = 1:dim(simMarketStocks)[2]
plot(y = cumsum(simMarketStocks[sample(400, size = 1), ]), x = x, type = 'l', main ='sample simulated stock', ylab = 'cummulative return', ylim = c(-400,400))
for (i in 1:100){
  lines(y = cumsum(simMarketStocks[sample(400, size = 1), ]), x = x, type = 'l', main ='sample simulated stock', ylab = 'cummulative return', col= distinctColorPalette(k = i, altCol = FALSE, runTsne = FALSE))
}

# visualise two states
par(mfrow = c(1,2))
plot(y = rowMeans(simMarketStocks[,spinLabels==1]), x = 1:400, type = 'h')
plot(y = rowMeans(simMarketStocks[,spinLabels==2]), x = 1:400, type = 'h')


mean(simMarketStocks[,spinLabels==1])
sd(simMarketStocks[,spinLabels==1])

mean(simMarketStocks[,spinLabels==2])
sd(simMarketStocks[,spinLabels==2])


write.csv(simMarketStocks, "data/Simulated/NohSim_400by200_gs0.1_etaGaussian.csv")

