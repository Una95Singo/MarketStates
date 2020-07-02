# ICC replication script, implemeneted using the EM algorithm version
# una singo, SNGUNA003
# 9 June 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/ICC_EM.R")
source("R/simulate.R")

# libraries ---------------------
library(readxl)
library(MASS)
library(fossil)
library(tidyverse)
library(lubridate)
library(NetworkToolbox)
#library(huge)
library(xts)
#library(corrplot)
library(PerformanceAnalytics)
#library(timeSeries)
#library(xts)

# helper functions -----

naSums = function(x){sum(is.na(x))}

set.seed(5)

# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
Ret.dates = GRet[,1]
GRet = GRet[,-1]
GRet = diff(as.matrix(log(GRet)))
colnames(GRet)

# test variables

#hist(abs(diag((GRet - colMeans(GRet)) %*% (LoGo(GRet)) %*% t(GRet- colMeans(GRet)))), xlab = 'distance', main ='Distribution of Mahalanobis distance')

dim(GRet)
par(mfrow = c(1,1))
plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')


# ICC setup
gamma = 0
sparseMethod = 2
distanceFunction = 1 
K =2


ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 100)

FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches
# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)

#Visualise states
State1 = which(FinalPath == 1)
par(mfrow=c(2,2))
# state 1
plot(colMeans(GRet[State1,]), type ='h', ylab='Mu', col ='red', lwd=2, main='mean stock return', ylim = c(-0.01, 0.01))
plot(apply(GRet[State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Sigma', col ='red', lwd=2, main='standard dev')
#state 2
plot(colMeans(GRet[-State1,]), type ='h', ylab='Mu', col ='blue', lwd=2, ylim= c(-0.01, 0.01))
plot(apply(GRet[-State1,], MARGIN = 2, FUN = sd), type ='h', ylab='Sigma', col ='blue', lwd=2)

# Sharpe ratios
par(mfrow=c(1,1))

data <- GRet
dates <-seq(as.Date("2013-02-09"), length = 1258, by = "days")
Close <- xts(x = data, order.by = dates)


plot(c(SharpeRatio.annualized(Close[State1,])), type ='h', ylab='Mu', col ='red', lwd=2, main='Sharpe ratios', ylim = c(2, -2))
lines(c(SharpeRatio.annualized(Close[-State1,])), type ='h', ylab='Mu', col ='blue', lwd=2)

plot(x = colMeans(GRet), y=apply(GRet, MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch =19, col='red' , xlim = c(-0.02, 0.02), ylim = c(0, 0.03))
plot(x = colMeans(GRet[State1,]), y=apply(GRet[State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch =19, col='red' , xlim = c(-0.02, 0.02), ylim = c(0, 0.03))
points(x = colMeans(GRet[-State1,]), y=apply(GRet[-State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch = 22, col ='blue')





