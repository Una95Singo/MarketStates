# Error testing
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
library(parallel)
library(foreach)
library(doParallel)



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
states = 4
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


# identified transition states
summary(hmmfit)

hmmTMatrix = matrix(getpars(hmmfit)[(nstates(hmmfit)+1):(nstates(hmmfit)^2+nstates(hmmfit))],
                    byrow=TRUE,nrow=nstates(hmmfit),ncol=nstates(hmmfit))
hmmStatePaths = posterior(hmmfit)$state

#State sample statistics
State1Returns = returns[hmmStatePaths==1]
State2Returns = returns[hmmStatePaths==2]
State3Returns = returns[hmmStatePaths==3]
State4Returns = returns[hmmStatePaths==4]

mean(State1Returns)
mean(State2Returns)
mean(State3Returns)
mean(State4Returns)
mean(returns)


sd(State2Returns)
sd(State1Returns)
sd(State3Returns)
sd(State4Returns)
sd(returns)
layout(1:1)

# simulate a fixed data set to identify a resonable gamma



DimData = list()

result = matrix(-99, nrow = length(No.iters), ncol = 10)
marketStates = c("1", "2", "3", "4")
byRow = TRUE
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = hmmTMatrix, name = "twoStateMarket")

param = list()
param[[1]] = list('mu' = mean(State1Returns),
                  'sigma' = sd(State1Returns))
param[[2]] = list('mu' = mean(State2Returns),
                  'sigma' = sd(State2Returns))

param[[3]] = list('mu' = mean(State3Returns),
                  'sigma' = sd(State3Returns))

param[[4]] = list('mu' = mean(State4Returns),
                  'sigma' = sd(State4Returns))

twoStateSequence = as(rmarkovchain(n = 600, object = twoState, t0 = '1'), 'numeric') 
stocks = c(100)




#----------------------------------------

# ICC results - gamma 0, no penalisation 2 state market
#----------------------------------------

# using SnP based hmm

# 100 stocks
it = 1

stockSize = 100
intra = c(0.25, 0.5, 0.75, 1)

ICC.RandIndex.Scene1 = matrix(NA, nrow =length(intra) , ncol = it)
ICC.Mean.Scene1 = matrix(NA, nrow =length(intra)*4 , ncol = it) # store mean estimates
ICC.Std.Scene1 = matrix(NA, nrow =length(intra)*4 , ncol = it) #store sd estimates
Noh.Mean.Scene1 = matrix(NA, nrow =length(intra)*4 , ncol = it) # store mean estimates
Noh.Std.Scene1 = matrix(NA, nrow =length(intra)*4 , ncol = it) #store sd estimates


row.names(ICC.RandIndex.Scene1) = intra
row.names(ICC.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2", "0.25_3", "0.5_3", "0.75_3", "1_3", "0.25_4", "0.5_4", "0.75_4", "1_4")
row.names(ICC.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2", "0.25_3", "0.5_3", "0.75_3", "1_3", "0.25_4", "0.5_4", "0.75_4", "1_4")

row.names(Noh.Mean.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2", "0.25_3", "0.5_3", "0.75_3", "1_3", "0.25_4", "0.5_4", "0.75_4", "1_4")
row.names(Noh.Std.Scene1) = c("0.25_1", "0.5_1", "0.75_1", "1_1", "0.25_2", "0.5_2", "0.75_2", "1_2", "0.25_3", "0.5_3", "0.75_3", "1_3", "0.25_4", "0.5_4", "0.75_4", "1_4")


dataNoh = NohAnzats.simulate(stateSequence = twoStateSequence, g = c(intra[i], intra[i], intra[i], intra[i]), clusters = 4, stocks = stockSize, gaussian = T, param)

plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
               , title = paste('Simlated market returns'))
plotter.series(colMeans((dataNoh$stockMatrix)),dataNoh$stateSeq
               , title = paste('Simlated market returns'), S0 = 100)


ICC.Output = ICC.cluster(returnsMatrix = dataNoh$stockMatrix, sparseMethod = 2, gamma = 0, K = 4, max.iters = 15)

State1Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==1]
State2Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==2]
State3Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==3]
State4Returns = colMeans((dataNoh$stockMatrix))[ICC.Output$OptimalPath==4]

ICC.RandIndex.Scene1[i,j] = adj.rand.index(dataNoh$stateSeq, ICC.Output$OptimalPath)

sparseMethod = 2
gamma = 0 
K = 4
max.iters = 15
distanceFunction = 1

returnsMatrix = dataNoh$stockMatrix
trace(print)
Time = ncol(returnsMatrix)
optimalViterbi = NA
optimalStateNum = 0
iters = 0
ErrorRate = 0
max.it = max.iters
stop.crit = 0.00001
ds = NA
d.hist= c()
#initialise cluster parameters
states.est = ICC.shuffle(K = K, Time = Time)
theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
#print( det(theta.est$precision[[1]]))
#print( det(theta.est$precision[[2]]))
#det(LoGo.solve(abs(cov(t(returnsMatrix[, which(states.est == 1)])))))
distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
#plot(y = distance.est[1,], x = 1:Time, type='l', col = 'red', main = 'initialisation')
#lines(y = distance.est[2,], x = 1:Time, col = 'blue')


optimalViterbi = viterbi(D = t(distance.est), K=K, gamma=gamma)
states.est = optimalViterbi$Final_Path
states.est

# check if assingment is valid. if not start with random assingment
if (length(unique(states.est)) != K || sum(diff(states.est)!=0)<3 ) {
  print(paste('Produced single state path at iteration', it))
  states.est = ICC.shuffle(K = K, Time = Time)
}

theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = 2)


distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
d.hist =c(d.hist, -1*ICC.TotalDistance(distance.est, gamma = gamma, stateSeq = states.est))

plot(y = distance.est[1,], x = 1:Time, type='l', col = 'red', main = paste('iteration: ', it))
lines(y = distance.est[2,], x = 1:Time, col = 'blue')
if (nrow(distance.est) == 3 ){ lines(y = distance.est[3,], x = 1:Time, col = 'black')}
if (nrow(distance.est) == 4 ){ lines(y = distance.est[3,], x = 1:Time, col = 'yellow')}
if (nrow(distance.est) == 5 ){ lines(y = distance.est[3,], x = 1:Time, col = 'pink')}

