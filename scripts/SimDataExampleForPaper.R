# Market Data Simulation tool
# Unarine Singo
# 30 January 2021


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
#source("R/stylisedFacts.R")
source("R/05_ICC.R")
source("R/NohAnzatsStateSimulation.R")

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
library(hash)

# setup
# data storage
states = 2
dimensions = 100
intraClusterStength = seq( from = 0.1, to = 1,length.out = 10)


# identify estimate parameters from real data using a hidden markov model
#helper functions -----
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
GRet = diff(as.matrix(log(GRet)))

# find estimates
returns = rowMeans(GRet)


hmm = depmix(returns ~ 1, family = gaussian(), nstates = states, data=data.frame(returns=returns))
hmmfit = fit(hmm, verbose = FALSE)
post_probs = posterior(hmmfit)

# Output both the true regimes and the 
# posterior probabilities of the regimes
layout(1:3)
plot(100*cumsum(returns)+100, type ='l', ylab = "Cummulative returns")
plot(returns, type ='l', ylab ="Log returns")
matplot(post_probs[,-1], type='l', main='Regime Posterior Probabilities', ylab='Probability')
legend(x='topright', c('State 1','State 2'), fill=1:2, bty='n')

# identified transition states
summary(hmmfit)

hmmTMatrix = matrix(getpars(hmmfit)[(nstates(hmmfit)+1):(nstates(hmmfit)^2+nstates(hmmfit))],
                    byrow=TRUE,nrow=nstates(hmmfit),ncol=nstates(hmmfit))
hmmStatePaths = posterior(hmmfit)$state

#State sample statistics
State1Returns = returns[hmmStatePaths==1]
State2Returns = returns[hmmStatePaths==2]

mean(State1Returns)
mean(State2Returns)
mean(returns)

sd(State2Returns)
sd(State1Returns)
sd(returns)
marketStates = c("1", "2")
byRow = TRUE
#simulate states
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix, name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 1825, object = twoState, t0 = '1'), 'numeric')

param = list()
param[[1]] = list('mu' = mean(State1Returns),
                  'sigma' = sd(State1Returns))
param[[2]] = list('mu' = mean(State2Returns),
                  'sigma' = sd(State2Returns))

param = list()
param[[1]] = list('mu' = 0,
                  'sigma' = 1)
param[[2]] = list('mu' = 0,
                  'sigma' = 1)

for(intra in intraClusterStength){
    g = c(intra, intra)

    stateSequence = twoStateSequence
    Noh = NohAnzats.simulate(stateSequence = twoStateSequence, g = g, clusters = 2, stocks = 100, gaussian = T, param)
    
    par(mfrow =c(2,1))
    # visualise identified path
    plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
                   , title = paste('Simlated market returns', "gs=",intra))
    # visualise identified path
    plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
                   , title = paste('Simulated market cummulative returns', "gs=",intra), S0 = 100)
    
    dev.copy(png,paste('images/NohAnzats examples/Standardexampleplot','_gs',intra,'_','.png',sep =""))
    dev.off()
}
