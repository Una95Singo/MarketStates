# Simple script to test applicaton of Mahalanobis distance function to a matrix of observations
# una singo, SNGUNA003
# 24 June 2020

# source functions

# libraries
library('tidyverse')
library(readxl)
library(fossil)
#library(xtable)
library(NetworkToolbox)
#library(ape)
#library(network)


# read data and clean

# helper function
naSums = function(x){sum(is.na(x))}


# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:300, size = 100)
GRet = survivorStocks[,-1]
#Ret.dates = GRet[,1]
#GRet = GRet[,c(4, 5, 6, 7, 8, 11, 15, 20, 30)]
GRet = GRet[, smaller]
GRet = diff(as.matrix(log(GRet)))


# Test 1 Full Precition Matrix
mu_1 = colMeans(GRet)
pre_1 = solve(cov(GRet))

full.dist = mahalanobis(GRet, center = mu_1, cov = pre_1, inverted = T)
plot(full.dist, type = 'l', ylab = 'Mahalanobis Distance', xlab = 'Time', main = 'Sample Mahalanobis Distance')

# Test 2 Sparse precition matrix
mu_2  = colMeans(GRet)
pre_2 = LoGo(cov(GRet), na.data = 'none', partial = F)

LoGo.dist = mahalanobis(GRet, center = mu_2, cov = pre_2, inverted = T)
lines(LoGo.dist, type = 'l', col = 'red', lty = 4)

legend(1, 500, legend=c("Full precision", "Sparse precision"),
       col=c("black", "red"), lty=1:2, cex=0.8)



