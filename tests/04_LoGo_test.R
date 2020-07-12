# 12 July 2020
# LoGo test
# Una Singo 

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/ICC_EM.R")
source("R/simulate.R")
source("R/04_LoGo.R")

# libraries ---------------------
library(readxl)
library(xlsx)
library(MASS)
library(tidyverse)
library(lubridate)
library(xts)

# helper functions -----

naSums = function(x){sum(is.na(x))}

# create test data 
# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period



# test set 1 ------- N = 25
testData = list()
sizes = c(25, 50, 75, 100)
for (i in 1:4){
  smaller = sample(2:400, size = sizes[i])
  GRet = survivorStocks[, c(smaller)]
  GRet = diff(as.matrix(log(GRet)))
  covariance = cov(GRet)
  testData[[i]] = covariance
  write.xlsx(covariance, 'data/Test/LoGotest.xlsx', sheetName=paste("datatest",i), 
             col.names=F, row.names=F, append = T )
}

# R LoGo implementation test -------

R.LoGo.results = list()
for( i in 1:4){
  R.LoGo.results[[i]] = LoGo(abs(testData[[i]]))
  SparseM::image(R.LoGo.results[[i]], main =paste(i), col = rainbow(10))
}

Matlab.LoGo.results = list()
for ( i in 1:4){
  Matlab.LoGo.results[[i]] = as.matrix(read.xlsx(file ="data/Test/LoGoMatlabResults.xlsx", sheetIndex = i, header = F))
}


par(mfrow = c(4, 2))
image(Matlab.LoGo.results[[1]], main ='Aste Matlab implementation, N = 25', col = topo.colors(10))
image(R.LoGo.results[[1]], main = 'Our implementation, N = 25', col = topo.colors(10))

image(Matlab.LoGo.results[[2]], main ='Aste Matlab implementation, N = 50', col = topo.colors(20))
image(R.LoGo.results[[2]], main = 'Our implementation, N = 50', col = topo.colors(20))

image(Matlab.LoGo.results[[3]], main ='Aste Matlab implementation, N = 75', col = topo.colors(20))
image(R.LoGo.results[[3]], main = 'Our implementation, N = 75', col = topo.colors(20))

image(Matlab.LoGo.results[[4]], main ='Aste Matlab implementation, N = 100', col = topo.colors(50))
image(R.LoGo.results[[4]], main = 'Our implementation, N = 100', col = topo.colors(50))


# plot determinants and condition number as size increases


dets = c(NA, NA, NA, NA)
conds = c(NA, NA, NA, NA)

for ( i in 1:4){
  dets[i] =  det(R.LoGo.results[[i]])
  conds[i] = kappa(R.LoGo.results[[i]])
}
# note to self. determinants of 100>= stocks is infinity. Why?