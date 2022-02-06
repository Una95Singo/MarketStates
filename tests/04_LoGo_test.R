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
sizes = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)
d_mat = matrix(NA, nrow = 100, ncol = 11 )

for(i in 1:11){
  for (j in 1:100){
    smaller = sample(2:400, size = 100)
    GRet = survivorStocks[, c(smaller)]
    GRet = diff(as.matrix(log(GRet)))
    covariance = cov(GRet)
    precision = LoGo.solve((covariance))
    d_mat[j,i] = det(precision)
    #testData[[i]] = covariance
    #write.xlsx(covariance, 'data/Test/LoGotest2.xlsx', sheetName=paste("datatest",i), 
    #           col.names=F, row.names=F, append = T )
  }
}

colMeans(d_mat)
plot(y =colMeans(d_mat), x = sizes, log='y', ylab = 'Determinant', main = 'LoGo determinant estimates', xlab = 'Number of stocks', type ='l')
points(y =colMeans(d_mat), x = sizes, col ='blue')

barplot(height= colMeans(d_mat)[1:2])

det(testData[[1]])
J = solve(testData[[1]])
L = LoGo.solve(testData[[1]])
det( J )
det( L )


det(testData[[3]])
J = solve(testData[[3]])
L = LoGo.solve(abs(testData[[3]]))
det( J )
det( L )
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


par(mfrow = c(1, 2))
image(Matlab.LoGo.results[[1]], main ='Aste Matlab implementation, N = 25', col = topo.colors(10))
image(Matlab.LoGo.results[[1]], main = 'Our implementation, N = 25', col = topo.colors(10))

image(Matlab.LoGo.results[[2]], main ='Aste Matlab implementation, N = 50', col = topo.colors(20))
image(R.LoGo.results[[2]], main = 'Our implementation, N = 50', col = topo.colors(20))

image(Matlab.LoGo.results[[3]], main ='Aste Matlab implementation, N = 75', col = topo.colors(20))
image(R.LoGo.results[[3]], main = 'Our implementation, N = 75', col = topo.colors(20))

image(Matlab.LoGo.results[[4]], main ='Aste Matlab implementation, N = 100', col = topo.colors(50))
image(R.LoGo.results[[4]], main = 'Our implementation, N = 100', col = topo.colors(50))


#-----------------------------------------------------------------
#---------------------- Simulate and plot determinants---------------------
#----------------------------------------------------------------

sims = 100
stockSizes = seq(from = 10, to = 150,by = 10)
determinantsError = c()
determinantsFix = c()


for ( j in stockSizes){
  Ret = matrix(NA, nrow = 1000 , ncol = j)
  for(i in 1:j){
    Ret[, i] = rnorm(1000, mean = 0, 0.014) 
  }
  determinantsError = (c( determinantsError , det(LoGo.solve(cov(Ret)))))
  determinantsFix = c( determinantsFix , sum(2*log(diag(chol( LoGo.solve(cov(Ret)))) ) ))
}

log(determinantsError)
determinantsFix

plot(y =(determinantsFix), x=stockSizes, main = "Precision Matrix Log-Determinants for N Simulated Stocks",
     xlab = 'Number of Stocks', 
     ylab = "Log-Determinant",
     type = 'b', 
     lwd = 2,
      ylim = c(min(determinantsFix), max(determinantsFix)))

lines(y =log(determinantsError), x = stockSizes, col ='red', type ='b', lwd=2)

legend(x=20, y =1200, legend = c('Direct','Choleskey'), col =c('red', 'black'), lty = c(2,2) )


Ret = matrix(NA, nrow = N, ncol = D)

var = cov(Ret)
var.inv = solve(var)
dim(var)

plot(y = (100*cumsum(rowMeans(Ret))+100), x = 1:N, type='l')




gbm_loop <- function(nsim = 100, t = 25, mu = 0, sigma = 0.1, S0 = 100, dt = 1./365) {
  gbm <- matrix(ncol = nsim, nrow = t)
  for (simu in 1:nsim) {
    gbm[1, simu] <- S0
    for (day in 2:t) {
      epsilon <- rnorm(1)
      dt = 1 / 365
      gbm[day, simu] <- gbm[(day-1), simu] * exp((mu - sigma * sigma / 2) * dt + sigma * epsilon * sqrt(dt))
    }
  }
  return(gbm)
}

library(tidyverse)
nsim <- 150
t <- 100
mu <- 0
sigma <- 0.1
S0 <- 100
gbm <- gbm_loop(nsim, t, mu, sigma, S0)

plot(y=gbm[,4], x=1:100, type = 'l')

Ret = diff(log(gbm))
det(LoGo.solve(cov(Ret)))

library(matrixcalc)
is.positive.semi.definite(cov(GRet))
Ret
