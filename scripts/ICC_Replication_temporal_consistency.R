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
library(huge)
library(xts)
library(corrplot)
library(PerformanceAnalytics)
library(timeSeries)
library(xts)

# helper functions -----

naSums = function(x){sum(is.na(x))}


# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period

# history holders
history.switches = matrix(NA, nrow = 4, ncol = 100)
history.consitency = matrix(NA, nrow = 4, ncol = 100)



# Full precission Gamma Zero ----------------------------

# ICC setup
gamma = 0
sparseMethod = 1
distanceFunction = 1 
K = 2

for(j in 1:100){
  #sample data
  smaller = sample(2:400, size = 100)
  GRet = survivorStocks[, c(1,smaller)]
  Ret.dates = GRet[,1]
  GRet = GRet[,-1]
  GRet = diff(as.matrix(log(GRet)))
  
  ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 20)
  
  FinalPath = ICC.Output$OptimalPath
  Number.switches = sum(diff(FinalPath) != 0)
  
  history.switches[1, j] = Number.switches
  # visualise identified path
  plotter.series(colMeans(t(GRet)),FinalPath
                 , title = paste('Optimal market states, gamma:', gamma), S0 = 100)
  
}
  
history.switches[1,]



# Full precission Gamma Zero ----------------------------

# ICC setup
gamma = 0.0001/2
sparseMethod = 1
distanceFunction = 1 
K = 2

for(j in 1:100){
  #sample data
  print(j)
  smaller = sample(2:400, size = 100)
  GRet = survivorStocks[, c(1,smaller)]
  Ret.dates = GRet[,1]
  GRet = GRet[,-1]
  GRet = diff(as.matrix(log(GRet)))
  
  ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 20)
  
  FinalPath = ICC.Output$OptimalPath
  Number.switches = sum(diff(FinalPath) != 0)
  
  history.switches[2, j] = Number.switches
  #visualise identified path
  #plotter.series(colMeans(t(GRet)),FinalPath
  #               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)
  
}


median(history.switches[1, ])
median(history.switches[2,])

history.switches[2,]

