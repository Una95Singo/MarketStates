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


# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period

# generate dataset
max.iterations = 5
DataSet = matrix (NA, nrow = (nrow(survivorStocks)-1), ncol =(100*max.iterations))
for (i in 1:max.iterations){
  smaller = sample(2:400, size = 100)
  GRet = survivorStocks[, c(smaller)]
  GRet = diff(as.matrix(log(GRet)))
  DataSet[, (i*100-99):(i*100)] = GRet
}



Gamma = seq(from = 0, to = 50, by = 2)

history.distance = matrix(NA, nrow = max.iterations, ncol = length(Gamma))
sparseMethod = 2
distanceFunction = 1 
K =2

for (i in 1:ncol(history.distance)){
  print(Gamma[i])
  for (j in 1:nrow(history.distance)){
    # ICC setup
    Gret = DataSet[, (j*100-99):(j*100)] 
    ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = Gamma[i], K = K, max.iters = 30)
    history.distance[j, i] = ICC.Output$OptimalViterbi$Final_Cost
  }
}


# 1. Open jpeg file
jpeg("images/GridSearchPlot.jpg", width = 350, height = 350)
#2. plot grid search results
plot(spline(y = colMeans(history.distance[,1:17]), x = Gamma[1:17]), type ='l', ylab = 'Average cost distance', xlab = 'Gamma', main = 'Grid search error')
# 3. Close the file
dev.off()


# gamma 24

smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
GRet = diff(as.matrix(log(GRet)))
ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = 24, K = K, max.iters = 30)
ICC.Output$OptimalViterbi$Final_Cost


FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches

# 1. Open jpeg file
jpeg("images/Gamma24.jpg", width = 350, height = 350)
#2. plot grid search results
# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', '24'), S0 = 100)
# 3. Close the file
dev.off()


ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = 5, K = K, max.iters = 30)
ICC.Output$OptimalViterbi$Final_Cost


FinalPath = ICC.Output$OptimalPath
Number.switches = sum(diff(FinalPath) != 0)
Number.switches

# 1. Open jpeg file
jpeg("images/Gamma5.jpg", width = 350, height = 350)
#2. plot grid search results
# visualise identified path
plotter.series(colMeans(t(GRet)),FinalPath
               , title = paste('Optimal market states, gamma:', '5'), S0 = 100)
# 3. Close the file
dev.off()
