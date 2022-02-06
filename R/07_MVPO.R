# MPT implementation
# una singo, SNGUNA003
# 18 July 2021

# package import


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
#source("R/stylisedFacts.R")
source("R/05_ICC.R")
source("R/NohAnzatsStateSimulation.R")
source('R/06_ASCP.R')
# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(markovchain)
library(depmixS4)
library(mclust)
library(xts)
library(PerformanceAnalytics)
library(openxlsx)
library(Rcpp)
library(timeSeries)
library(nloptr)
library(reticulate)



#setFinCenter(tsTAA) <- "Johannesburg"


#helper functions
naSums = function(x){sum(is.na(x))}

# clean and prep data -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
GRet = diff(as.matrix(log(GRet)))
GRet = timeSeries(GRet, as.Date(survivorStocks$date[-1]))

# timeseries object with dates
Mu = colMeans(GRet)*12
Sigma = cov(GRet)*12

#plot data
plot(apply(GRet, MARGIN = 2, FUN = cumsum), plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Returns",
     main = "Random Sample of SnP 500 stock returns (n=100)")

plot(GRet, plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Returns",
     main = "Random Sample of SnP 500 stock returns (n=100)")


plot(apply(GRet, MARGIN = 1, FUN = mean), plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Returns",
     main = "Random Sample of SnP average market return (n=100)")


States = replicate(nrow(GRet), 1)

par(mfrow=c(1,1))



#------------------------------------------------------------------------
#---------------------------------- Data setup ----------------------------
#------------------------------------------------------------------------

# ICC covariance comparison for the methodology -------------------------------------
ICC.Output = ICC.cluster(returnsMatrix =t(GRet), sparseMethod = 2, gamma = 10, K = 2, max.iters = 50)

States = ICC.Output$OptimalPath
plotter.series(colMeans((t(GRet))),States
               , title = paste('Identified Market States'))

par(mfrow=c(2,2))
plotter.series(colMeans((t(GRet))),States
               , title = paste('Identified Market States'))
image(cov(GRet[States==1,]), main = "Covariance mat: state 1", col = hcl.colors(150, rev = TRUE), useRaster = T)
image(cov(GRet[States==2,]), main = "Covariance: state 2",col = hcl.colors(150, rev = TRUE),useRaster = T)
image(cov(GRet), main = "Coveriance: no states",col = hcl.colors(150, rev = TRUE),useRaster = T)


par(mfrow=c(1,1))



# split data into training and test -------------------------------------

trainIndex = round(nrow(GRet)[1]*0.7)

train = GRet[c(1:trainIndex), ]
trainStates = States[c(1:trainIndex)]
test = GRet[-c(1:trainIndex), ]
testStates = States[-c(1:trainIndex)]



#------------------------------------------------------------------------
#---------------------------------- ICC ----------------------------
#------------------------------------------------------------------------
# compute train sample stats
mu1 = matrix(colMeans(train[trainStates==1,]),nrow=100, ncol=1)*12
mu2 = matrix(colMeans(train[trainStates==2,]), nrow=100, ncol=1)*12
sigma1 = cov(train[trainStates==1,])*12
sigma2 = cov(train[trainStates==2,]) *12

# single states
mu3 = matrix(colMeans(train),nrow=100, ncol=1)*12
sigma3 = cov(train) *12


# estimate state  weights -------------------------------------

# MVPO state 1
# timeseries object with dates
#mu1 = ICC.Output$ThetaEst$mu[[1]]*12
#sigma1 = ICC.Output$ThetaEst$precision[[1]]*12
Wts1 = MVOptimise(Mu = mu1, Sigma = sigma1)
#Wts1 = rep(0,100)
#Wts1[1] = 0.5
#Wts1[2] =0.5


# MVPO state 2
#mu2 = ICC.Output$ThetaEst$mu[[2]]*12
#sigma2 = ICC.Output$ThetaEst$precision[[2]]*12
Wts2 = MVOptimise(Mu = mu2, Sigma = sigma2)


euclidean(Wts1, Wts2)

#build weight matrix
LT_WtsA = matrix(NA, nrow = Dim(test)[1], ncol =Dim(test)[2])

for (i in (1:length(testStates))){
  if(testStates[i]==1){
    LT_WtsA[i, ] = Wts1
  } else if (testStates[i]==2){
    LT_WtsA[i, ] = Wts2
  }
}

LT_WtsA = timeSeries(LT_WtsA, rownames(test))

PortfolioReturnsA = Return.portfolio(
  test, 
  weights = LT_WtsA,
  geometric = T,
  wealth.index = F,
  value = 1000
)

PortfolioWealthA = Return.portfolio(
  test, 
  weights = LT_WtsA,
  geometric = F,
  wealth.index = T,
  value = 1000
)

# single state comparison

#mu3 = colMeans(GRet)*12
#sigma3 = cov(GRet)*12
Wts3 = MVOptimise(Mu = mu3, Sigma = sigma3)
#StatesSingle = rep(1, nrow(test))

LT_WtsB = matrix(NA, nrow = nrow(test), ncol =Dim(GRet)[2])

for (i in (1:nrow(test))){
    LT_WtsB[i, ] = Wts3
}

LT_WtsB = timeSeries(LT_WtsB, rownames(test))


euclidean(Wts1, Wts3)
euclidean(Wts2, Wts3)

PortfolioWealthB = Return.portfolio(
  test, 
  weights = LT_WtsB,
  geometric = F,
  wealth.index = T,
  value = 1000
)

PortfolioReturnsB = Return.portfolio(
  test, 
  weights = LT_WtsB,
  geometric = T,
  wealth.index = F,
  value = 1000
)



#------------------------------------------------------------------------
#---------------------------------- ASPC ----------------------------
#------------------------------------------------------------------------

py_run_file("R/Python/ASPCLionel.py")
it = 1
test = (GRet-rowMeans(GRet))
den = apply(GRet, 1, sd)
G_ = cor(t(test/matrix(den,1258,100)))

G_ = cor(t(GRet))
dim(G_)

#constrained
ASPC.Output =  agglo_spc(G = G_, cn = 2)


States = ASPC.Output
tail(States)
#States[330:340] = 1
#States[1000:1100] = 1
States
plotter.series(colMeans((t(GRet))),States
               , title = paste('Identified Market States'))


train = GRet[c(1:trainIndex), ]
trainStates = States[c(1:trainIndex)]
test = GRet[-c(1:trainIndex), ]
testStates = States[-c(1:trainIndex)]

mu4 = matrix(colMeans(train[trainStates==1,]),nrow=100, ncol=1)*12
mu5 = matrix(colMeans(train[trainStates==2,]), nrow=100, ncol=1)*12
sigma4 = cov(train[trainStates==1,])*12
sigma5 = cov(train[trainStates==2,]) *12


Wts4 = MVOptimise(Mu = mu4, Sigma = sigma4)

Wts5 = MVOptimise(Mu = mu5, Sigma = sigma5)


LT_WtsC = matrix(NA, nrow = Dim(test)[1], ncol =Dim(test)[2])

for (i in (1:length(testStates))){
  if(testStates[i]==1){
    LT_WtsC[i, ] = Wts4
  } else if (testStates[i]==2){
    LT_WtsC[i, ] = Wts5
  }
}

LT_WtsC = timeSeries(LT_WtsC, rownames(test))

PortfolioWealthC = Return.portfolio(
  test, 
  weights = LT_WtsC,
  geometric = F,
  wealth.index = T,
  value = 1000
)

PortfolioReturnsC = Return.portfolio(
  test, 
  weights = LT_WtsC,
  geometric = T,
  wealth.index = F,
  value = 1000
)


# unconstrained ASPC--------------

#constrained
ASPC.Output =  agglo_spc(G = G_, cn=5)
ASPC.Output2 =  agglo_spc(G = G_, cn=6)
ASPC.Output3 =  agglo_spc(G = G_, cn=7)


train = GRet[c(1:trainIndex), ]
trainStates = States[c(1:trainIndex)]
test = GRet[-c(1:trainIndex), ]
testStates = States[-c(1:trainIndex)]

mu6 = matrix(colMeans(train[trainStates==1,]),nrow=100, ncol=1)*12
mu7 = matrix(colMeans(train[trainStates==2,]), nrow=100, ncol=1)*12
mu8 = matrix(colMeans(train[trainStates==3,]),nrow=100, ncol=1)*12
mu9 = matrix(colMeans(train[trainStates==4,]), nrow=100, ncol=1)*12
mu10 = matrix(colMeans(train[trainStates==5,]),nrow=100, ncol=1)*12

sigma6 = cov(train[trainStates==1,])*12
sigma7 = cov(train[trainStates==2,]) *12
sigma8 = cov(train[trainStates==1,])*12
sigma9 = cov(train[trainStates==2,]) *12
sigma10 = cov(train[trainStates==1,])*12


Wts6 = MVOptimise(Mu = mu6, Sigma = sigma6)
Wts7 = MVOptimise(Mu = mu7, Sigma = sigma7)
Wts8 = MVOptimise(Mu = mu7, Sigma = sigma7)
Wts9 = MVOptimise(Mu = mu7, Sigma = sigma7)
Wts10 = MVOptimise(Mu = mu7, Sigma = sigma7)


LT_WtsD = matrix(NA, nrow = Dim(test)[1], ncol =Dim(test)[2])

for (i in (1:length(testStates))){
  if(testStates[i]==1){
    LT_WtsD[i, ] = Wts6
  } else if (testStates[i]==2){
    LT_WtsD[i, ] = Wts7
  }else if (testStates[i]==3){
    LT_WtsD[i, ] = Wts8
  }else if (testStates[i]==4){
    LT_WtsD[i, ] = Wts9
  }else if (testStates[i]==5){
    LT_WtsD[i, ] = Wts10
  }
}

LT_WtsD = timeSeries(LT_WtsD, rownames(test))


plotter.series(colMeans((t(GRet))),States
               , title = paste('Identified Market States'), S0=100)

PortfolioWealthD = Return.portfolio(
  test, 
  weights = LT_WtsD,
  geometric = F,
  wealth.index = T,
  value = 1000
)

PortfolioReturnsD = Return.portfolio(
  test, 
  weights = LT_WtsD,
  geometric = T,
  wealth.index = F,
  value = 1000
)


States = ASPC.Output
tail(States)
#States[330:340] = 1
#States[1000:1100] = 1
States
plotter.series(colMeans((t(GRet))),States
               , title = paste('Identified Market States'))




min_ =min(c(PortfolioWealthA, PortfolioWealthB, PortfolioWealthC ))
max_=max(c(PortfolioWealthA, PortfolioWealthB, PortfolioWealthC ))
# plot results
plot(PortfolioWealthA, plot.type = c("single"), format = "auto",
     at = pretty(GRet),
     ylab = "Returns",
     ylim = c(min_,max_),
     col = 'blue',
     main = "MVPO total wealth from investment strategies")

lines(PortfolioWealthB,
     ylab = "Returns",
     main = "Single state", col = 'black')

lines(PortfolioWealthC,
      ylab = "Returns",
      main = "Single state", col = 'red')

lines(PortfolioWealthD,
      ylab = "Returns",
      main = "Single state", col = 'green')

legend("topleft", inset=0.15, legend = c("ICC MVPO wealth", "Standard MVPO wealth", "ASPC MVPO wealth", "Unconstrained ASPC MVPO"), col=c('blue', 'black', 'red', 'green'), lty=1:1, cex=0.8)


SharpeRatio(R = PortfolioReturnsA)
SharpeRatio(R =PortfolioReturnsB)
SharpeRatio(R = PortfolioReturnsC)
SharpeRatio(R = PortfolioReturnsD)


mean(PortfolioReturnsA)/sd(PortfolioReturnsA)
mean(PortfolioReturnsB)/sd(PortfolioReturnsB)


#------------------------------------------------------------------------
#---------------------------------- ASPC ----------------------------
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#---------------------------------- FUNCTION ----------------------------
#------------------------------------------------------------------------

MVOptimise = function (Mu = NA, Sigma = NA){
  
  # annualise statistics
  Mu = 12 * Mu
  Sigma = 12 * Sigma
  
  # full invested
  A = cbind(1, diag(length(Mu)))
  b = 1
  
  # no short selling
  #A = cbind(1, diag(length(Mu)))
  #b = c(b , rep(0, length(Mu)))
  
  #initialisse weights
  #Wts = matrix (NA, 1)

  #Find The Sharpe Ratio maximising portfolio
  # Initial values as fully invested equally weighted portfolio
  Ones = seq(1,1, length.out = length (Mu))
  # equally weighted portfolio
  Wts0 <- Ones/length(Ones)
  # unit vector
  e <- rep(1, length ( Wts0)) # useful matrix ( ones )
  # initialise the weights
  Wts <- matrix ( NA ,1, length (Mu))
  # 8. Maximise the Sharpe Ratio
  fn0 <- function ( x ) { return ( -(x%*%Mu)/ sqrt ( x %*% Sigma %*% x ) )}
  # Fully Invested + Return Target
  heq0 <- function ( x ) { return ( x %*% e - 1)} # fully invested
  # Use SQP to solve for the tangency portfolio
  soln=slsqp( Wts0, fn = fn0, gr = NULL , # target returns
                  lower = rep(0, length ( Wts0)) , # no short - selling
                  upper = rep(1, length ( Wts0)) , # no leverage
                  heq = heq0, # fully invested constraint function
                  control = list (xtol_rel = 1e-8)) # SQP
  Wts = soln$par
  return(Wts)
}

