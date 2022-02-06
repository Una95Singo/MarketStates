# Market Data Simulation
# Unarine Singo
# 3 January 2021


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
StateSets = list('S1'=NA,'S2'=NA, 'S3'=NA, 'S4'=NA, 'S5'=NA)
dimensions = c(10, 25, 50, 75, 100, 200, 300)
No.iters = 1:10
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
#returns = colMeans(Noh$stockMatrix)


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
mean(State2Returns)
mean(returns)


sd(State2Returns)
sd(State1Returns)
sd(returns)



#------------------------------------ 2 State Market Simulation------------------------------------
#--------------------------------------------------------------------------------------------------

marketStates = c("1", "2")
byRow = TRUE
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix[[2]], name = "twoStateMarket")
twoStateSequence = as(rmarkovchain(n = 1825, object = twoState, t0 = '1'), 'numeric')


Datasets = list()
DimHash = hash()
for(dim in dimensions){
  .set(DimHash, keys = dim, values = hash())
  for(intra in intraClusterStength){
    for (it in No.iters){
      #simulate states
      twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix[[2]], name = "twoStateMarket")
      twoStateSequence = as(rmarkovchain(n = 1825, object = twoState, t0 = '1'), 'numeric')
    
      
      param = list()
      param[[1]] = list('mu' = mean(State1Returns),
                        'sigma' = sd(State1Returns))
      param[[2]] = list('mu' = mean(State2Returns),
                        'sigma' = sd(State2Returns))
      
      g = c(intra, intra)
      stocks = Dim
      gaussian = T
      stateSequence = twoStateSequence
      Noh = NohAnzats.simulate(stateSequence = twoStateSequence, g = g, clusters = 2, stocks = dim, gaussian = T, param)
      
      
      par(mfrow =c(2,1))
      # visualise identified path
      plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
                     , title = paste('Simlated market returns'))
      # visualise identified path
      plotter.series(colMeans((Noh$stockMatrix)),twoStateSequence
                     , title = paste('Simulated market cummulative returns'), S0 = 100)
      
      Datasets[[it]] = Noh
    }
    .set(DimHash[[paste(dim)]], keys = intra, values = Datasets)
  }
}

# store dataset into a RDatabase
saveRDS(DimHash, file = "data/Simulated/simulated2StateNohData.rds")


#------------------------------------ Identify market states w simulated datasets

par(mfrow=c(1,1))
# ICC setup
gamma = 0
sparseMethod = 2
distanceFunction = 1 
K =2


ICC.Output = ICC.cluster(returnsMatrix = DimHash[['130']][['0.1']][[1]]$stockMatrix, sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 25)

FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches

par(mfrow =c(2,1))
# visualise identified path
plotter.series(colMeans(DimHash[['130']][['0.1']][[1]]$stockMatrix),FinalPath
               , title = paste('Optimal market states, gamma:', gamma))

# visualise identified path
plotter.series(colMeans(DimHash[['130']][['0.1']][[1]]$stockMatrix),FinalPath
               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)

adjustedRandIndex(DimHash[['130']][['0.325']][[1]]$stateSeq, FinalPath)



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

