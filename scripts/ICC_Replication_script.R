# ICC replication script
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
library(matrixcalc)

#library(timeSeries)
#library(xts)

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


# test variables

#hist(abs(diag((GRet - colMeans(GRet)) %*% (LoGo(GRet)) %*% t(GRet- colMeans(GRet)))), xlab = 'distance', main ='Distribution of Mahalanobis distance')

dim(GRet)
par(mfrow = c(1,1))
plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')
plot(1:1258, y =colMeans(t(GRet)) , type='l', ylab = 'daily return', xlab = 'time', main ='SnP 500 sample return')

par(mfrow=c(1,1))
# ICC setup
gamma = 16
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

plot(main = 'Individual stock sample statistics',x = colMeans(GRet[State1,]), y=apply(GRet[State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch =19, col='red' , xlim = c(1.5*min_x, 1.2*max_x), ylim = c(0.7*min_y, 1.2*max_y))
points(x = colMeans(GRet[-State1,]), y=apply(GRet[-State1,], MARGIN = 2, FUN = sd), xlab='mu', ylab='sd', pch = 22, col ='blue')


# out of sample prediction ------------->


N = nrow (GRet)
p = ncol(GRet)
split = round(N*0.650)
train = GRet[1:split,]
test = GRet[-c(1:split),]

ISPath = FinalPath[1:split]
OOSPath = FinalPath[-c(1:split)]



# in sample training

gamma = 16
sparseMethod = 2
distanceFunction = 1 
K =2

IS.ICC.Output = ICC.cluster(returnsMatrix = t(train), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 30)

# calculate mu and precision
IS.Est = ICC.thetaEst(K=2, returns = t(train), stateSeq = IS.ICC.Output$OptimalPath, sparse = 2 )



# forecasting
d.squared = ICC.distance(t(GRet), IS.Est, K = 2, dist = 1)
p = nrow(IS.Est$precision[[1]])

Liklihood.state1 = 1/2 * (sum(2*log(diag(chol(IS.Est$precision[[1]])))) - d.squared[1,] - p*log(2*pi) )
Liklihood.state1
Liklihood.state2 =  1/2 * ( sum(2*log(diag(chol(IS.Est$precision[[2]])))) - d.squared[2,] - p*log(2*pi) )
Liklihood.state2

plot(Liklihood.state1, type ='l')
lines(Liklihood.state2, col = 'red', )
abline(v = split, col ='green')


roll.delta = 24 
Rt = c()

for( t in roll.delta:N){
  s = t-roll.delta+1
  Rt =c( Rt, sum(Liklihood.state1[s:t] - Liklihood.state2[s:t]) )
}



plot(Rt, type ='l', col = 'red')
abline(h = 0)


y = IS.ICC.Output$OptimalPath[roll.delta:split] - 1

logistic.data = cbind(Rt[1:(split-24+1)],y)
colnames(logistic.data) = c('R_t', 'Y_t')

fitted.model = glm(formula = Y_t~R_t, data = as.data.frame(logistic.data ), family = 'binomial')


IS.predicted = predict(fitted.model, data = as.data.frame(logistic.data), type ='response')

IS.predicted.classes =  ifelse(IS.predicted > 0.5, 1, 0) + 1

# in-sample proportion table
prop.table(table(IS.predicted.classes, y+1), 1)

confMatrix = table(IS.predicted.classes, y+1)
sum(diag(confMatrix))/sum(confMatrix)

# out-of-sample

OOS.Rt = Rt[-c(1:(split-24+1))]
length(OOS.Rt)
#OOS.predicted = predict(fitted.model, data = OOS.Rt, type ='response')
OOS.predicted = 1/(1+exp(-1*(fitted.model$coefficients[1] + OOS.Rt*fitted.model$coefficients[2])))
length(OOS.predicted)
OOS.predicted.classes =  ifelse(OOS.predicted > 0.50, 1, 0) + 1

prop.table(table(OOS.predicted.classes, OOSPath), 1)
confMatrix = table(OOS.predicted.classes, OOSPath)
sum(diag(confMatrix))/sum(confMatrix)

