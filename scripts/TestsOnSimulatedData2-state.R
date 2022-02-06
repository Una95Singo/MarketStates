# Tests on simulated datasets
# Unarine Singo
# 6 January 2021


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
#source("R/stylisedFacts.R")
source("R/05_ICC.R")
source("R/NohAnzatsStateSimulation.R")
source('R/06_ASCP.R')

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

library(reticulate)



# ---------------------------------------------
# Setup ---------------------------------------
# ---------------------------------------------

# identify estimate parameters from real data using a hidden markov model
#helper functions
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





# ---------------------------------------------
# Simulate 2 state data -----------------------
# ---------------------------------------------
dimensions = 100
states = 2
No.iters = 1:5
intraClusterStength = seq( from = 0.1, to = 1,length.out = 10)


# setup --------------
hmm = depmix(returns ~ 1, family = gaussian(), nstates = states, data=data.frame(returns=returns))
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


# simulate a fixed data set to identify a resonable gamma



DimData = list()

result = matrix(-99, nrow = length(No.iters), ncol = 10)
marketStates = c("1", "2")
byRow = TRUE
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix, name = "twoStateMarket")

param = list()
param[[1]] = list('mu' = mean(State1Returns),
                  'sigma' = sd(State1Returns))
param[[2]] = list('mu' = mean(State2Returns),
                  'sigma' = sd(State2Returns))

twoStateSequence = as(rmarkovchain(n = 600, object = twoState, t0 = '1'), 'numeric') 
stocks = c(100)


# Gamma search - testing
#-------------------------------
# sim data
# gamma search, iterate 100 times and store
it = 100
gamma = c(0, 1, 10, 50, 100)
intra = 1

RandIndex1 = matrix(NA, nrow =length(gamma) , ncol = it )
row.names(RandIndex1) = gamma
for (i in 1:length(gamma)){
  for (j in 1:it){
  print(paste("Gamma iteration:", gamma[i], "of", j ))
  dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra, intra), clusters = 2, stocks = 100, gaussian = T, param)
  
  plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                 , title = paste('Simlated market returns'))
  plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                 , title = paste('Simlated market returns'), S0 = 100)
  
  ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma =gamma[i], K = 2, max.iters = 15)
  
  RandIndex1[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
  }
}

save(RandIndex1, file ="data/GammaSearchScene1.RData")



intra = 0.75

RandIndex2 = matrix(NA, nrow =length(gamma) , ncol = it )
row.names(RandIndex2) = gamma
for (i in 1:length(gamma)){
  for (j in 1:it){
    print(paste("Gamma iteration:", gamma[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra, intra), clusters = 2, stocks = 100, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma =gamma[i], K = 2, max.iters = 15)
    
    RandIndex2[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
  }
}


save(RandIndex2, file ="data/GammaSearchScene2.RData")


intra = 0.5

RandIndex3 = matrix(NA, nrow =length(gamma) , ncol = it )
row.names(RandIndex2) = gamma
for (i in 1:length(gamma)){
  for (j in 1:it){
    print(paste("Gamma iteration:", gamma[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra, intra), clusters = 2, stocks = 100, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma =gamma[i], K = 2, max.iters = 15)
    
    RandIndex3[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
  }
}

save(RandIndex3, file ="data/GammaSearchScene3.RData")

intra = 0.25

RandIndex4 = matrix(NA, nrow =length(gamma) , ncol = it )
row.names(RandIndex2) = gamma
for (i in 1:length(gamma)){
  for (j in 1:it){
    print(paste("Gamma iteration:", gamma[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra, intra), clusters = 2, stocks = 100, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma =gamma[i], K = 2, max.iters = 15)
    
    RandIndex4[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
  }
}
save(RandIndex4, file ="data/GammaSearchScene4.RData")


# plot the gamma search results
layout(1,1)
plot(y=rowMeans(RandIndex1), x=gamma, type='l', main ='Gamma search', xlab = "Gamma", ylab ="Adj. Rand Index", ylim=c(0,1.6), log='x')
points(y=rowMeans(RandIndex1), x=gamma, pch = 15)
abline(v = gamma, col = 'red', lty =3, lwd = 0.5)

lines(y=rowMeans(RandIndex2), x=gamma, col ='blue')
points(y=rowMeans(RandIndex2), x=gamma, pch = 16, col='blue')

lines(y=rowMeans(RandIndex3), x=gamma, col ='green')
points(y=rowMeans(RandIndex3), x=gamma, pch = 17, col='green')

lines(y=rowMeans(RandIndex4), x=gamma, col ='pink')
points(y=rowMeans(RandIndex4), x=gamma, pch = 18, col='pink')

# Add a legend
legend(y = 1.6, x = 1, legend=c("Intra = 1","Intra = 0.75","Intra = 0.5","Intra = 0.25"),
       col=c("black","blue","green", "pink"), lty=1, cex=0.8)

#----------------------------------------

# ICC results - gamma 0, no penalisation 2 state market
#----------------------------------------

# using SnP based hmm

# 100 stocks
it = 100

stockSize = 100
intra = c(0.25, 0.5, 0.75, 1)

ICC.RandIndex.Scene1 = matrix(NA, nrow =length(intra) , ncol = it)
ICC.Mean.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ICC.Std.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ICC.RandIndex.Scene1) = intra
row.names(ICC.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ICC.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

system.time({
for (i in 1:length(intra)){
  for (j in 1:it){
    print(paste("Intracluster iteration:", intra[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma = 0, K = 2, max.iters = 15)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==2]
    
    ICC.RandIndex.Scene1[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
    ICC.Mean.Scene1[i, j] = round(mean(State1Returns),6)
    ICC.Mean.Scene1[i+4, j] = round(mean(State2Returns),6)
    ICC.Std.Scene1[i, j] = round(sd(State1Returns),6)
    ICC.Std.Scene1[i+4, j] =  round(sd(State2Returns), 6)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
    
    Noh.Mean.Scene1[i, j] = round(mean(State1Returns),6)
    Noh.Mean.Scene1[i+4, j] = round(mean(State2Returns),6)
    Noh.Std.Scene1[i, j] = round(sd(State1Returns),6)
    Noh.Std.Scene1[i+4, j] =  round(sd(State2Returns), 6)
  }
}})

save(ICC.RandIndex.Scene1, ICC.Mean.Scene1 , ICC.Std.Scene1, file ="data/ICCScene1.RData")
load(file ="data/ICCScene1.RData")

rowMeans(ICC.RandIndex.Scene1)


# 250 stocks
stockSize = 250
intra = c(0.25, 0.5, 0.75, 1)

ICC.RandIndex.Scene2 = matrix(NA, nrow =length(intra) , ncol = it)
ICC.Mean.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ICC.Std.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ICC.RandIndex.Scene2) = intra
row.names(ICC.Mean.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ICC.Std.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")


for (i in 1:length(intra)){
  for (j in 1:it){
    print(paste("Intracluster iteration:", intra[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma = 0, K = 2, max.iters = 15)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==2]
    
    ICC.RandIndex.Scene2[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
    ICC.Mean.Scene2[i, j] = round(mean(State1Returns),6)
    ICC.Mean.Scene2[i+4, j] = round(mean(State2Returns),6)
    ICC.Std.Scene2[i, j] = round(sd(State1Returns),6)
    ICC.Std.Scene2[i+4, j] =  round(sd(State2Returns), 6)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
    
    Noh.Mean.Scene2[i, j] = round(mean(State1Returns),6)
    Noh.Mean.Scene2[i+4, j] = round(mean(State2Returns),6)
    Noh.Std.Scene2[i, j] = round(sd(State1Returns),6)
    Noh.Std.Scene2[i+4, j] =  round(sd(State2Returns), 6)
  }
}


save(ICC.RandIndex.Scene2, ICC.Mean.Scene2 , ICC.Std.Scene2, file ="data/ICCScene2.RData")
#load(file ="data/ICCScene1.RData")
#save(ICC.RandIndex.Scene2 , file ="data/GammaSearch.RData")
#save(ICC.Mean.Scene2 , file ="data/GammaSearch.RData")
#save(ICC.Std.Scene2 , file ="data/GammaSearch.RData")

# 500 stocks
stockSize = 500
intra = c(0.25, 0.5, 0.75, 1)

ICC.RandIndex.Scene3 = matrix(NA, nrow =length(intra) , ncol = it)
ICC.Mean.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ICC.Std.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ICC.RandIndex.Scene3) = intra
row.names(ICC.Mean.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ICC.Std.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")


for (i in 1:length(intra)){
  for (j in 1:it){
    print(paste("Intracluster iteration:", intra[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma = 0, K = 2, max.iters = 15)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==2]
    
    ICC.RandIndex.Scene3[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
    ICC.Mean.Scene3[i, j] = round(mean(State1Returns),6)
    ICC.Mean.Scene3[i+4, j] = round(mean(State2Returns),6)
    ICC.Std.Scene3[i, j] = round(sd(State1Returns),6)
    ICC.Std.Scene3[i+4, j] =  round(sd(State2Returns), 6)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
    
    Noh.Mean.Scene3[i, j] = round(mean(State1Returns),6)
    Noh.Mean.Scene3[i+4, j] = round(mean(State2Returns),6)
    Noh.Std.Scene3[i, j] = round(sd(State1Returns),6)
    Noh.Std.Scene3[i+4, j] =  round(sd(State2Returns), 6)
  }
}

save(ICC.RandIndex.Scene3, ICC.Mean.Scene3 , ICC.Std.Scene3, file ="data/ICCScene3.RData")
#save(ICC.RandIndex.Scene3 , file ="data/GammaSearch.RData")
#save(ICC.Mean.Scene3 , file ="data/GammaSearch.RData")
#save(ICC.Std.Scene3 , file ="data/GammaSearch.RData")

# 750 stocks

stockSize = 750
intra = c(0.25, 0.5, 0.75, 1)

ICC.RandIndex.Scene4 = matrix(NA, nrow =length(intra) , ncol = it)
ICC.Mean.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ICC.Std.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ICC.RandIndex.Scene4) = intra
row.names(ICC.Mean.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ICC.Std.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")


for (i in 1:length(intra)){
  for (j in 1:it){
    print(paste("Intracluster iteration:", intra[i], "of", j ))
    dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
    
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'))
    plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                   , title = paste('Simlated market returns'), S0 = 100)
    
    ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma = 0, K = 2, max.iters = 15)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==2]
    
    ICC.RandIndex.Scene4[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)
    ICC.Mean.Scene4[i, j] = round(mean(State1Returns),6)
    ICC.Mean.Scene4[i+4, j] = round(mean(State2Returns),6)
    ICC.Std.Scene4[i, j] = round(sd(State1Returns),6)
    ICC.Std.Scene4[i+4, j] =  round(sd(State2Returns), 6)
    
    State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
    State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
    
    Noh.Mean.Scene4[i, j] = round(mean(State1Returns),6)
    Noh.Mean.Scene4[i+4, j] = round(mean(State2Returns),6)
    Noh.Std.Scene4[i, j] = round(sd(State1Returns),6)
    Noh.Std.Scene4[i+4, j] =  round(sd(State2Returns), 6)
  }
}


save(ICC.RandIndex.Scene4, ICC.Mean.Scene4 , ICC.Std.Scene4, file ="data/ICCScene4.RData")
#save(ICC.RandIndex.Scene4 , file ="data/GammaSearch.RData")
#save(ICC.Mean.Scene4 , file ="data/GammaSearch.RData")
#save(ICC.Std.Scene4 , file ="data/GammaSearch.RData")

load(file ="data/ICCScene1.RData")
load(file ="data/ICCScene2.RData")
load(file ="data/ICCScene3.RData")
load(file ="data/ICCScene4.RData")

# plot the gamma search results
layout(1,1)
plot(y=rowMeans(ICC.RandIndex.Scene1), x=intra, type='l', main ='ICC performance results', xlab = "Intra-Cluster Strength", ylab ="Adj. Rand Index", ylim=c(0,1.6))
points(y=rowMeans(ICC.RandIndex.Scene1), x=intra, pch = 15)
abline(v = intra, col = 'red', lty =3, lwd = 0.5)

lines(y=rowMeans(ICC.RandIndex.Scene2), x=intra, col ='blue')
points(y=rowMeans(ICC.RandIndex.Scene2), x=intra, pch = 16, col='blue')

lines(y=rowMeans(ICC.RandIndex.Scene3), x=intra, col ='green')
points(y=rowMeans(ICC.RandIndex.Scene3), x=intra, pch = 17, col='green')

lines(y=rowMeans(ICC.RandIndex.Scene4), x=intra, col ='pink')
points(y=rowMeans(ICC.RandIndex.Scene4), x=intra, pch = 18, col='pink')

# Add a legend
legend(y = 1.6, x = 0.25, legend=c("Stocks = 100","Stocks = 250","Stocks = 500","Stocks = 750"),
       col=c("black","blue","green", "pink"), lty=1, cex=0.8)


round(rowMeans(ICC.RandIndex.Scene1),2)
round(apply(X =ICC.RandIndex.Scene1,MARGIN = 1,FUN = sd),3)
round(rowMeans(ICC.RandIndex.Scene2),2)
round(apply(X =ICC.RandIndex.Scene2,MARGIN = 1,FUN = sd),3)
round(rowMeans(ICC.RandIndex.Scene3),2)
round(apply(X =ICC.RandIndex.Scene3,MARGIN = 1,FUN = sd),2)
round(rowMeans(ICC.RandIndex.Scene4),2)
round(apply(X =ICC.RandIndex.Scene4,MARGIN = 1,FUN = sd),2)




#----------------
# ASPC results - gamma 0, no penalisation 2 state market# 100 stocks
#----------
py_run_file("R/Python/ASPCLionel.py")
it = 1

stockSize = 100
intra = c(0.25, 0.5, 0.75, 1)

ASPC.RandIndex.Scene1 = matrix(NA, nrow =length(intra) , ncol = it)
ASPC.Mean.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ASPC.Std.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene1 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ASPC.RandIndex.Scene1) = intra
row.names(ASPC.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ASPC.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

system.time({
  for (i in 1:length(intra)){
    for (j in 1:it){
      print(paste("Intracluster iteration:", intra[i], "of", j ))
      dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
      
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'))
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'), S0 = 100)
      
      test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix))/ apply(dataNoh$stockMatrix, 1, sd)
    
      #write.csv(x = test, file = paste("R/Python/SimData/",intra[i], "_", j ,".csv",sep = ""))

      #G_ = cor((dataNoh$stockMatrix))
      G_ = cor(test)
      #G_
      
      #py_run_file("R/Python/ASPCLionel.py")
      ASPC.Output =  agglo_spc(G = G_, cn = 2)
 
      adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      State1Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==2]
      
      ASPC.RandIndex.Scene1[i,j] = adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      ASPC.Mean.Scene1[i, j] = round(mean(State1Returns),6)
      ASPC.Mean.Scene1[i+4, j] = round(mean(State2Returns),6)
      ASPC.Std.Scene1[i, j] = round(sd(State1Returns),6)
      ASPC.Std.Scene1[i+4, j] =  round(sd(State2Returns), 6)
      
      State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
      
      Noh.Mean.Scene1[i, j] = round(mean(State1Returns),6)
      Noh.Mean.Scene1[i+4, j] = round(mean(State2Returns),6)
      Noh.Std.Scene1[i, j] = round(sd(State1Returns),6)
      Noh.Std.Scene1[i+4, j] =  round(sd(State2Returns), 6)
    }
  }})

save(ASPC.RandIndex.Scene1, ASPC.Mean.Scene1 , ASPC.Std.Scene1, file ="data/ASPCScene1_2S.RData")


stockSize = 250
intra = c(0.25, 0.5, 0.75, 1)

ASPC.RandIndex.Scene2 = matrix(NA, nrow =length(intra) , ncol = it)
ASPC.Mean.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ASPC.Std.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene2 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ASPC.RandIndex.Scene2) = intra
row.names(ASPC.Mean.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ASPC.Std.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene2) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

system.time({
  for (i in 1:length(intra)){
    for (j in 1:it){
      print(paste("Intracluster iteration:", intra[i], "of", j ))
      dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
      
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'))
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'), S0 = 100)
      
      test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix))/ apply(dataNoh$stockMatrix, 1, sd)
      #test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix)) / apply(dataNoh$stockMatrix, 2, sd)
      #write.csv(x = test, file = paste("R/Python/SimData/",intra[i], "_", j ,".csv",sep = ""))
      
      #G_ = cor((dataNoh$stockMatrix))
      G_ = cor(test)
      #G_
      
      #py_run_file("R/Python/ASPCLionel.py")
      ASPC.Output =  agglo_spc(G = (G_), cn = 2)
      
      adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      State1Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==2]
      
      ASPC.RandIndex.Scene2[i,j] = adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      ASPC.Mean.Scene2[i, j] = round(mean(State1Returns),6)
      ASPC.Mean.Scene2[i+4, j] = round(mean(State2Returns),6)
      ASPC.Std.Scene2[i, j] = round(sd(State1Returns),6)
      ASPC.Std.Scene2[i+4, j] =  round(sd(State2Returns), 6)
      
      State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
      
      Noh.Mean.Scene2[i, j] = round(mean(State1Returns),6)
      Noh.Mean.Scene2[i+4, j] = round(mean(State2Returns),6)
      Noh.Std.Scene2[i, j] = round(sd(State1Returns),6)
      Noh.Std.Scene2[i+4, j] =  round(sd(State2Returns), 6)
    }
  }})

save(ASPC.RandIndex.Scene2, ASPC.Mean.Scene2 , ASPC.Std.Scene2, file ="data/ASPCScene2_2S.RData")

stockSize = 500
intra = c(0.25, 0.5, 0.75, 1)

ASPC.RandIndex.Scene3 = matrix(NA, nrow =length(intra) , ncol = it)
ASPC.Mean.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ASPC.Std.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene3 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ASPC.RandIndex.Scene3) = intra
row.names(ASPC.Mean.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ASPC.Std.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene3) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

system.time({
  for (i in 1:length(intra)){
    for (j in 1:it){
      print(paste("Intracluster iteration:", intra[i], "of", j ))
      dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
      
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'))
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'), S0 = 100)
      
      test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix))/ apply(dataNoh$stockMatrix, 1, sd)
      
      #test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix)) / apply(dataNoh$stockMatrix, 2, sd)
      #write.csv(x = test, file = paste("R/Python/SimData/",intra[i], "_", j ,".csv",sep = ""))
      
      #G_ = cor((dataNoh$stockMatrix))
      G_ = cor(test)
      #G_
      
      #py_run_file("R/Python/ASPCLionel.py")
      ASPC.Output =  agglo_spc(G = (G_), cn = 2)
      
      adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      State1Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==2]
      
      ASPC.RandIndex.Scene3[i,j] = adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      ASPC.Mean.Scene3[i, j] = round(mean(State1Returns),6)
      ASPC.Mean.Scene3[i+4, j] = round(mean(State2Returns),6)
      ASPC.Std.Scene3[i, j] = round(sd(State1Returns),6)
      ASPC.Std.Scene3[i+4, j] =  round(sd(State2Returns), 6)
      
      State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
      
      Noh.Mean.Scene3[i, j] = round(mean(State1Returns),6)
      Noh.Mean.Scene3[i+4, j] = round(mean(State2Returns),6)
      Noh.Std.Scene3[i, j] = round(sd(State1Returns),6)
      Noh.Std.Scene3[i+4, j] =  round(sd(State2Returns), 6)
    }
  }})

save(ASPC.RandIndex.Scene3, ASPC.Mean.Scene3 , ASPC.Std.Scene3, file ="data/ASPCScene3_2S.RData")

stockSize = 750
intra = c(0.25, 0.5, 0.75, 1)

ASPC.RandIndex.Scene4 = matrix(NA, nrow =length(intra) , ncol = it)
ASPC.Mean.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
ASPC.Std.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates
Noh.Mean.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) # store mean estimates
Noh.Std.Scene4 = matrix(NA, nrow =length(intra)*2 , ncol = it) #store sd estimates


row.names(ASPC.RandIndex.Scene4) = intra
row.names(ASPC.Mean.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(ASPC.Std.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

row.names(Noh.Mean.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")
row.names(Noh.Std.Scene4) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2")

system.time({
  for (i in 1:length(intra)){
    for (j in 1:it){
      print(paste("Intracluster iteration:", intra[i], "of", j ))
      dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i]), clusters = 2, stocks = stockSize, gaussian = T, param)
      
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'))
      plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
                     , title = paste('Simlated market returns'), S0 = 100)
      
      test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix))/ apply(dataNoh$stockMatrix, 1, sd)
      #test = (dataNoh$stockMatrix - rowMeans(dataNoh$stockMatrix)) / apply(dataNoh$stockMatrix, 2, sd)
      write.csv(x = test, file = paste("R/Python/SimData/",intra[i], "_", j ,".csv",sep = ""))
      
      #G_ = cor((dataNoh$stockMatrix))
      G_ = cor(test)
      #G_
      
      #py_run_file("R/Python/ASPCLionel.py")
      ASPC.Output =  agglo_spc(G = (G_), cn = 2)
      
      adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      State1Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[ASPC.Output==2]
      
      ASPC.RandIndex.Scene4[i,j] = adj.rand.index(dataNoh$stateSeq, ASPC.Output)
      ASPC.Mean.Scene4[i, j] = round(mean(State1Returns),6)
      ASPC.Mean.Scene4[i+4, j] = round(mean(State2Returns),6)
      ASPC.Std.Scene4[i, j] = round(sd(State1Returns),6)
      ASPC.Std.Scene4[i+4, j] =  round(sd(State2Returns), 6)
      
      State1Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==1]
      State2Returns = colMeans((dataNoh$stockMatrix))[dataNoh$stateSeq==2]
      
      Noh.Mean.Scene4[i, j] = round(mean(State1Returns),6)
      Noh.Mean.Scene4[i+4, j] = round(mean(State2Returns),6)
      Noh.Std.Scene4[i, j] = round(sd(State1Returns),6)
      Noh.Std.Scene4[i+4, j] =  round(sd(State2Returns), 6)
    }
  }})

save(ASPC.RandIndex.Scene4, ASPC.Mean.Scene4 , ASPC.Std.Scene4, file ="data/ASPCScene4_2S.RData")


load(file ="data/ASPCScene1_2S.RData")
load(file ="data/ASPCScene2_2S.RData")
load(file ="data/ASPCScene3_2S.RData")
load(file ="data/ASPCScene4_2S.RData")



# 1. Open jpeg file
#jpeg("images/Findings Chapter/ASPCPerformance3State.jpeg", width = 350, height = "350")
# plot the gamma search results
layout(1,1)
plot(y=rowMeans(ASPC.RandIndex.Scene1), x=intra, type='l', main ='ASPC performance results', xlab = "Intra-Cluster Strength", ylab ="Adj. Rand Index", ylim=c(0,1.6))
points(y=rowMeans(ASPC.RandIndex.Scene1), x=intra, pch = 15)
abline(v = intra, col = 'red', lty =3, lwd = 0.5)

lines(y=rowMeans(ASPC.RandIndex.Scene2), x=intra, col ='blue')
points(y=rowMeans(ASPC.RandIndex.Scene2), x=intra, pch = 16, col='blue')

lines(y=rowMeans(ASPC.RandIndex.Scene3), x=intra, col ='green')
points(y=rowMeans(ASPC.RandIndex.Scene3), x=intra, pch = 17, col='green')

lines(y=rowMeans(ASPC.RandIndex.Scene4), x=intra, col ='pink')
points(y=rowMeans(ASPC.RandIndex.Scene4), x=intra, pch = 18, col='pink')


round(rowMeans(ASPC.RandIndex.Scene1),2)
round(apply(X =ASPC.RandIndex.Scene1,MARGIN = 1,FUN = sd),2)
round(rowMeans(ASPC.RandIndex.Scene2),2)
round(apply(X =ASPC.RandIndex.Scene2,MARGIN = 1,FUN = sd),2)
round(rowMeans(ASPC.RandIndex.Scene3),2)
round(apply(X =ASPC.RandIndex.Scene3,MARGIN = 1,FUN = sd),2)
round(rowMeans(ASPC.RandIndex.Scene4),2)
round(apply(X =ASPC.RandIndex.Scene4,MARGIN = 1,FUN = sd),2)

# Add a legend
legend(y = 1.6, x = 0.25, legend=c("Stocks = 100","Stocks = 250","Stocks = 500","Stocks = 750"),
       col=c("black","blue","green", "pink"), lty=1, cex=0.8)
