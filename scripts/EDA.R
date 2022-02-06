# Exploratory data analysis
# una singo, SNGUNA003
# 18 July 2021

# package import


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(depmixS4)
library(xts)
library(PerformanceAnalytics)
library(openxlsx)
library(Rcpp)
library(timeSeries)




#helper functions
naSums = function(x){sum(is.na(x))}

# clean and prep data -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocksVol = allStocks  %>% dplyr::select(date, volume, Name)
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
flatStocksVol = allStocksVol %>% spread(key = Name, value = volume, fill = NA ) 
head(flatStocks)
tail(flatStocks)
survivorStocks = flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
survirvorVolume = flatStocksVol  %>% select_if(apply(flatStocksVol,naSums, MARGIN = 2) == 0)
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
GRet = diff(as.matrix(log(GRet)))
GRet = timeSeries(GRet, as.Date(survivorStocks$date[-1]))

survivorGRet = diff(as.matrix(log(survivorStocks[,-1])))
survivorGRet = timeSeries(survivorGRet, as.Date(survivorStocks$date[-1]))


survirvorVolume = timeSeries(survirvorVolume[,-1],as.Date(survirvorVolume$date))

#plot data single stocks
plot(apply(GRet, MARGIN = 2, FUN = cumsum), plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Cummulative Returns",
     main = "Random Sample of S & P 500 stock returns (n=100)")

plot(GRet, plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Daily Returns",
     main = "Random Sample of S & P 500 stock returns (n=100)")

#plot the market
par(mfrow=c(2,1))
plot(apply(apply(GRet, MARGIN = 1, FUN = mean), MARGIN = 2, FUN = cumsum), plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Cummulative Returns",
     main = "Random Sample of S & P 500 stock returns (n=100)")

plot(apply(GRet, MARGIN = 1, FUN = mean), plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Daily Returns",
     main = "Random Sample S & P 500 stock returns (n=100)")

# full sample
par(mfrow=c(2,2))
plot(apply(apply(survivorGRet, MARGIN = 1, FUN = mean), MARGIN = 2, FUN = cumsum), plot.type = c("single"), format = "auto",
     at = pretty(survivorGRet),
     ylab = "Cummulative Returns",
     main = "S & P 500 stock returns (n=470)")

plot(apply(survivorGRet, MARGIN = 1, FUN = mean), plot.type = c("single"), format = "auto",
     at = pretty(survivorGRet),
     ylab = "Daily Returns",
     main = "S & P 500 stock returns (n=470)")

plot(apply(survirvorVolume, 1, sum)/100000, main ="Traded volume", ylab = "# of daily trades (Mn)")
hist(apply(survivorGRet, MARGIN = 1, FUN = mean), main="Daily returns historgram", xlab = "Daily log returns")
