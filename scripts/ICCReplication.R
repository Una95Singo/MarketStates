# ICC replication
# una singo, SNGUNA003
# 9 June 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/ICC.R")
source("R/simulate.R")

# libraries ---------------------
library(readxl)
library(MASS)
library(fossil)
library(tidyverse)
library(lubridate)
library(NetworkToolbox)
library(huge)
library(xts)
library(corrplot)
library(PerformanceAnalytics)
library(timeSeries)
library(xts)


# helper functions -----

naSums = function(x){sum(is.na(x))}

set.seed(1)

# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(1,smaller)]
Ret.dates = GRet[,1]
GRet = GRet[,-1]
GRet = diff(as.matrix(log(GRet)))
colnames(GRet)

# test variables

#hist(abs(diag((GRet - colMeans(GRet)) %*% (LoGo(GRet)) %*% t(GRet- colMeans(GRet)))), xlab = 'distance', main ='Distribution of Mahalanobis distance')

dim(GRet)
par(mfrow = c(1,1))
plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')

GAMMA = 0.005
# Replication test 1 - identify two market states (Bull & Bear) --------------------
ICC.results = ICC.cluster(returnsMatrix = t(GRet), gamma = GAMMA, sparseMethod = 2, distanceFunction = 1, K =2, max.iters = 3)
FinalPath = ICC.results$OptimalPath
# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', GAMMA), S0 = 100)

#Visualise states
State1 = which(FinalPath == 1)
par(mfrow=c(2,2))
# state 1
plot(colMeans(GRet[State1,]), type ='h', ylab='Mu', col ='red', lwd=2, main='mean stock return', ylim = c(-0.0014, 0.0017))
plot(apply(GRet[State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Sigma', col ='red', lwd=2, main='standard dev')
#state 2
plot(colMeans(GRet[-State1,]), type ='h', ylab='Mu', col ='blue', lwd=2, ylim= c(-0.0014, 0.0017))
plot(apply(GRet[-State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Sigma', col ='blue', lwd=2)

# Sharpe ratios
par(mfrow=c(1,1))

data <- GRet
dates <-seq(as.Date("2013-02-09"), length = 1258, by = "days")
Close <- xts(x = data, order.by = dates)


plot(c(SharpeRatio.annualized(Close[State1,])), type ='h', ylab='Mu', col ='red', lwd=2, main='Sharpe ratios')
lines(c(SharpeRatio.annualized(Close[-State1,])), type ='h', ylab='Mu', col ='blue', lwd=2)


plot(x= colMeans(GRet[State1,]), y=apply(GRet[State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch =19, col='red' )
points(x= colMeans(GRet[-State1,]), y=apply(GRet[-State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch = 22, col ='blue')



Gamma = seq(0, 0.01, length=50)
# ----- Model 1 Sparse Precision matrix and temporal consistency -----
par(mfrow = c(1,1))
# Grid search first
viterbiSets = list()
performance.history = matrix(NA, 100, 100)
for (i in 1:length(Gamma)){
  print(paste("gamma", Gamma[i]))
  for (j in 1:100){
    ICC.results = ICC.cluster(returnsMatrix = t(GRet), gamma = Gamma[i], sparseMethod = 2, distanceFunction = 1, K =2, max.iters = 3)
    #viterbiSets[[i]] = ICC.results$OptimalViterbi
    performance.history[i,j] = ICC.results$OptimalViterbi$Final_Cost
  }
}

plot(x = Gamma, y = performance.history, type ='line')
seg.results = segmentation.procedure(returns =t(GRet), Time = 1258, gamma=0.0001, iters = 5, K = 2, sparse = 2 )
minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))
plotter.series(colMeans(t(GRet)),seg.results$history[[paste("iter:", minIndex[1],"states")]]
               , title = 'Optimal market states', S0 = 100)



# -----Test 2, test consistency of state identification -----

# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period

State1SR = list()
State2SR = list()
GAMMA = 0.005
iterations = 100

for (iter in 1:iterations){
  smaller = sample(2:400, size = 100)
  GRet = survivorStocks[, c(1,smaller)]
  GRet = GRet[,-1]
  GRet = diff(as.matrix(log(GRet)))
  par(mfrow = c(1,1))
  #plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')
  # Replication test 1 - identify two market states (Bull & Bear) --------------------
  ICC.results = ICC.cluster(returnsMatrix = t(GRet), gamma = GAMMA, sparseMethod = 2, distanceFunction = 1, K =2, max.iters = 3)
  FinalPath = ICC.results$OptimalPath
  State1 = which(FinalPath == 1)
  
  # Sharpe ratios
  par(mfrow=c(1,1))
  plot(colMeans(GRet[State1,])/apply(GRet[State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Mu', col ='red', lwd=2, main='Sharpe ratios')
  lines(colMeans(GRet[-State1,])/apply(GRet[-State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Mu', col ='blue', lwd=2)
  # store SharpeRatios
  State1SR = colMeans(GRet[State1,])/apply(GRet[State1,], MARGIN = 2, FUN = sd)
  State2SR = colMeans(GRet[-State1,])/apply(GRet[-State1,], MARGIN = 2, FUN = sd)
}


# -----Test 3, Forecasting -----



