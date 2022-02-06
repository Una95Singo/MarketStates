# ICC replication script, implemeneted using the EM algorithm version
# una singo, SNGUNA003
# 9 June 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/05_ICC.R")

# libraries ---------------------
library(readxl)
library(MASS)
library(fossil)
library(tidyverse)
library(lubridate)
#library(NetworkToolbox)
#library(huge)
library(xts)
#library(corrplot)
library(PerformanceAnalytics)
#library(timeSeries)
#library(xts)

# helper functions -----

naSums = function(x){sum(is.na(x))}
set.seed(1)

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
colnames(GRet)


# test variables

#hist(abs(diag((GRet - colMeans(GRet)) %*% (LoGo(GRet)) %*% t(GRet- colMeans(GRet)))), xlab = 'distance', main ='Distribution of Mahalanobis distance')

dim(GRet)
par(mfrow = c(1,1))
plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')
plot(1:1258, y =colMeans(t(GRet)) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')

par(mfrow=c(1,1))
# ICC setup
gamma = 17
sparseMethod = 2
distanceFunction = 1 
K =2


ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 50)

FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches

par(mfrow =c(1,2))
# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', gamma))

# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)


#Visualise states
State1 = which(FinalPath == 1)
par(mfrow=c(2,2))
# state 1
means1 = colMeans(GRet[State1,])
plot(means1, type ='h', ylab='Mu', col ='red', lwd=2, main='mean stock return', ylim = c(min(means1)*1.2, max(means1)*1.2))
sd1 = apply(GRet[State1,], MARGIN = 2, FUN = sd)
plot(sd1, type ='h', ylab='Sigma', col ='red', lwd=2, main='standard dev', ylim = c(min(sd1)*1.1, max(sd1)*1.1))
#state 2
means2 = colMeans(GRet[-State1,])
sd2 = apply(GRet[-State1,], MARGIN = 2, FUN = sd)
plot(means2, type ='h', ylab='Mu', col ='blue', lwd=2, ylim = c(min(means2)*1.2, max(means2)*1.2))
plot(sd2, type ='h', ylab='Sigma', col ='blue', lwd=2,  ylim = c(min(sd2)*1.1, max(sd2)*1.1) )




mean(means2-means1)
# Sharpe ratios
par(mfrow=c(1,1))

data <- GRet
dates <-seq(as.Date("2013-02-09"), length = 1258, by = "days")
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


# ---------
par(mfrow = c(1,1))
min_x = min(colMeans(GRet))
max_x = max(colMeans(GRet))
min_y = min(apply(GRet, MARGIN = 2, FUN = sd))
max_y = max(apply(GRet, MARGIN = 2, FUN = sd))

plot(x = colMeans(GRet[State1,]), y=apply(GRet[State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch =19, col='red' , xlim = c(1.5*min_x, 1.2*max_x), ylim = c(0.7*min_y, 1.2*max_y))
points(x = colMeans(GRet[-State1,]), y=apply(GRet[-State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch = 22, col ='blue')





