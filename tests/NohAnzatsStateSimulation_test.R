# Noh Anzats State Simulation
# Unarine Singo
# 8 June 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/stylisedFacts.R")
source("R/05_ICC.R")
#source("R/NohAnzatsStateSimulation.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(markovchain)
library(depmixS4)
library(mclust)
library(xts)
library(PerformanceAnalytics)

NohAnzats.simulate = function (stateSequence ='', g = c(0.1, 0.1), clusters = 2, stocks = 100, gaussian = T, param){
  # step 1 Cluster numbers --------------------------------------------------------
  #stateSequence = twoStateSequence
  spinLabels = stateSequence
  C = clusters  # two clusters of stocks
  s = max(table(stateSequence)) # clusters size is 50 each
  N = length(spinLabels) # number of days
  D = stocks  # number of stocks
  
  #N = stocks  # number of days
  #D = length(spinLabels)  # number of stocks
  # step 2 simulate cluster effects --------------------------------------------------------
 
  eta = matrix( NA, C, D, byrow = T) 
  if (gaussian){
    for (i in 1:C){
      eta[i, ] = rnorm (D, param[[i]]$mu, param[[i]]$sigma)
    }  
  } else{
    # cluster effects use t-distribution with 3 degrees of freedom
    for (i in 1:C){
      eta[i, ] = rt (D, df = param[[i]]$df)
    }
  }
  
  # step 3 simulate random effects ---------------------------------------------------------
  eps = matrix(NA, N, D) # specific random effects
  # epsilon
  for (n in 1:N) {
    eps[n, ] = rnorm(D, 0, 0.01)
  }
  
  # test market return
  #par(mfrow =  c(1,1))
  #plot(cumsum(colMeans(t(eps))), type='l', main ='Cummulative plot of epsilon', ylab= 'price')
  # step 4 fix Intra-cluster binding strength --------------------------------------------------------
  gs = g # same cluster binding strength
  
  # step 5 compute daily returns --------------------------------------------------------
  xi = matrix(NA, N, D) # matrix to store daily returns

  # double for loop for simualation
  for (n in 1:N){
    for (d in 1:D){
      cluster = spinLabels[n]
      #xi[n, d] = (sqrt(gs[cluster]) * eta[cluster, d] + eps[n, d]) / (sqrt(1+gs[cluster]))
      xi[n, d] = (gs[cluster] * eta[cluster, d]) + (eps[n, d] * sqrt(1+gs[cluster]^2))
    }
  }
  
  return(list( 'stockMatrix' = t(xi), 'stateSeq' = spinLabels, 'eta' = eta, 'eps' = eps, 'intraBinding' =gs))
}



#------------TEST SIMULATION 1
#----------------------------
marketStates = c("1", "2")
byRow = TRUE
tMatrix = matrix(data = c(0.918, 0.082,
                          0.059, 0.941), byrow = byRow, nrow = 2, dimnames = list(marketStates, marketStates))
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = tMatrix, name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 1000, object = twoState, t0 = '1'), 'numeric')
plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)

param = list()
#param[[1]] = list('mu' = 0.001, 'sigma' = 1)
#param[[2]] = list('mu' =  0, 'sigma' = 1)

param[[1]] = list('mu' = 0.00025, 'sigma' = 0.0106)
param[[2]] = list('mu' =  -0.00091, 'sigma' = 0.0402)

g = c(0.001, 0.001)
clusters = 2
stocks = 100
gaussian = T
stateSequence = twoStateSequence


Noh = NohAnzats.simulate(stateSequence = twoStateSequence, g = g, clusters = 2, stocks = 100, gaussian = T, param)
#plot(rowMeans(t(Noh$stockMatrix)), type ='l')
#plot(cumsum(rowMeans(t(Noh$stockMatrix))), type ='l')

par(mfrow =c(1,1))
# visualise identified path
plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
               , title = paste('Simlated market returns'))
# visualise identified path
plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
               , title = paste('Simulated market cummulative returns'), S0 = 100)

#------------------------------------ TEST END
#------------------------------------


#------------TEST SIMULATION 2 with ICC

#------------------------------------ Identify parameters start


# helper functions -----

naSums = function(x){sum(is.na(x))}
#set.seed(1)

# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
GRet
GRet = diff(as.matrix(log(GRet)))
colnames(GRet)

# find estimates
returns = rowMeans(GRet)

NumberOfStates = c(2,3,4,5)
hmmTMatrix = list()
hmmStatePaths = list()

for (i in NumberOfStates){
  # Create and fit the Hidden Markov Model
  hmm = depmix(returns ~ 1, family = gaussian(), nstates = i, data=data.frame(returns=returns))
  hmmfit = fit(hmm, verbose = FALSE)
  post_probs = posterior(hmmfit)
  
  # Output both the true regimes and the 
  # posterior probabilities of the regimes
  layout(1:3)
  plot(100*cumsum(returns)+100, type ='l')
  plot(returns, type ='l')
  matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')
  legend(x='topright', c('Bull','Bear'), fill=1:2, bty='n')
  
  # identified transition states
  summary(hmmfit)
  
  hmmTMatrix[[i]] = matrix(getpars(hmmfit)[(nstates(hmmfit)+1):(nstates(hmmfit)^2+nstates(hmmfit))],
                          byrow=TRUE,nrow=nstates(hmmfit),ncol=nstates(hmmfit))
  
  hmmStatePaths[[i]] = posterior(hmmfit)$state
}



#State sample statistics
State1Returns = returns[hmmStatePaths[[2]]==1]
State2Returns = returns[hmmStatePaths[[2]]==2]

mean(State1Returns)
sd(State1Returns)

mean(State2Returns)
sd(State2Returns)


mean(returns)
sd(returns)



#------------------------------------ Simulate datasets

marketStates = c("1", "2")
byRow = TRUE
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix[[2]], name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 1825, object = twoState, t0 = '1'), 'numeric')

layout(1:2)
plotter.states(hmmStatePaths[[2]], title = 'SnP identified states', thickness = 4)
plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)


param = list()
#param[[1]] = list('mu' = 0.001, 'sigma' = 1)
#param[[2]] = list('mu' =  0, 'sigma' = 1)

param[[1]] = list('mu' = mean(State1Returns),
                  'sigma' = sd(State1Returns))
param[[2]] = list('mu' = mean(State2Returns),
                   'sigma' = sd(State2Returns))

g = c(0.3, 0.3)
clusters = 2
stocks = 100
gaussian = T
stateSequence = twoStateSequence


Noh = NohAnzats.simulate(stateSequence = twoStateSequence, g = g, clusters = 2, stocks = 100, gaussian = T, param)
#plot(rowMeans(t(Noh$stockMatrix)), type ='l')
#plot(cumsum(rowMeans(t(Noh$stockMatrix))), type ='l')

par(mfrow =c(2,1))
# visualise identified path
plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
               , title = paste('Simlated market returns'))
# visualise identified path
plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
               , title = paste('Simulated market cummulative returns'), S0 = 100)



#State sample statistics
SimState1Returns = (Noh$stockMatrix)[twoStateSequence==1]
SimState2Returns = (Noh$stockMatrix)[twoStateSequence==2]

mean(SimState1Returns)
sd(SimState1Returns)

mean(SimState2Returns)
sd(SimState2Returns)

mean((Noh$stockMatrix))
sd((Noh$stockMatrix))

#------------------------------------ Identify market states w simulated datasets

par(mfrow=c(1,1))
# ICC setup
gamma = 16
sparseMethod = 2
distanceFunction = 1 
K =2


ICC.Output = ICC.cluster(returnsMatrix = Noh$stockMatrix, sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 50)

FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches

par(mfrow =c(2,1))
# visualise identified path
plotter.series(colMeans(Noh$stockMatrix),FinalPath
               , title = paste('Optimal market states, gamma:', gamma))

# visualise identified path
plotter.series(colMeans(Noh$stockMatrix),FinalPath
               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)

adjustedRandIndex(twoStateSequence, FinalPath)

#Visualise states
State1 = which(FinalPath == 1)
par(mfrow=c(2,2))
# state 1
means1 = colMeans(t(Noh$stockMatrix)[State1,])
plot(means1, type ='h', ylab='Mu', col ='red', lwd=2, main='mean stock return', ylim = c(min(means1)*1.2, max(means1)*1.2))
sd1 = apply(t(Noh$stockMatrix)[State1,], MARGIN = 2, FUN = sd)
plot(sd1, type ='h', ylab='Sigma', col ='red', lwd=2, main='standard dev', ylim = c(min(sd1), max(sd1)*1.1))
#state 2
means2 = colMeans(t(Noh$stockMatrix)[-State1,])
sd2 = apply(t(Noh$stockMatrix)[-State1,], MARGIN = 2, FUN = sd)
plot(means2, type ='h', ylab='Mu', col ='blue', lwd=2, ylim = c(min(means2)*1.2, max(means2)*1.2))
plot(sd2, type ='h', ylab='Sigma', col ='blue', lwd=2,  ylim = c(min(sd2), max(sd2)*1.1) )


# Sharpe ratios
par(mfrow=c(1,1))

data <- t(Noh$stockMatrix)
dates <-seq(as.Date("2013-02-09"), length = dim(data)[1], by = "days")
Close <- xts(x = data, order.by = dates)


SR1 = SharpeRatio.annualized(Close[State1,])
SR2 = SharpeRatio.annualized(Close[-State1,])
min_y = min(c(SR1, SR2))*1.1
max_y = max(c(SR1, SR2))*1.1          
plot(y=SR1, x=1:100, type ='h', ylab='Mu', col ='red', lwd=4, main='Sharpe ratios', ylim = c(min_y, max_y))
lines(y= SR2, x=1:100, type ='h', ylab='Mu', col ='blue', lwd=2)

par(mfrow = c(1,1))
hist(SR2, col=rgb(0,0,1,0.5), breaks = 10)
hist(SR1, col=rgb(1,0,0,0.5), breaks = 10, add=T)

