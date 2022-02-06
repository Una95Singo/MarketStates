# ICC implementation
# una singo, SNGUNA003
# 27 January 2020

# libraries ----------
source('R/01_viterbi.R')
source('R/04_LoGo.R')
library('Matrix')

#returnsMatrix = t(GRet)
#gamma = 16
#sparseMethod = 2
#distanceFunction = 1
#K =2

#returnsMatrix = t(GRet)
#gamma = 16
#sparseMethod = 2
#distanceFunction = 1 
#K =2
#set.seed(1)

ICC.cluster  = function(returnsMatrix, gamma = 1, sparseMethod = 1, distanceFunction = 1, K, max.iters = 50){
  #' ICC procedure
  #' @param returns matrix of reutrns stocks by time
  #' @param K number of states
  #' @param Time Time horizon of the dataset
  #' @param Sparse type of sparsifiation method
  #' @param dist which distance metric to use
  trace(print)
  Time = ncol(returnsMatrix)
  optimalViterbi = NA
  optimalStateNum = 0
  iters = 0
  ErrorRate = 0
  max.it = max.iters
  stop.crit = 0.00001
  ds = NA
  d.hist= c()
  #initialise cluster parameters
  states.est = ICC.shuffle(K = K, Time = Time)
  theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
  #print( det(theta.est$precision[[1]]))
  #print( det(theta.est$precision[[2]]))
  #det(LoGo.solve(abs(cov(t(returnsMatrix[, which(states.est == 1)])))))
  distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
  #plot(y = distance.est[1,], x = 1:Time, type='l', col = 'red', main = 'initialisation')
  #lines(y = distance.est[2,], x = 1:Time, col = 'blue')
  
  
  for (it in 1:max.it){
    #print("trap 2")
    # E - step // assin points to clusters using viterbi
    optimalViterbi = viterbi(D = t(distance.est), K=K, gamma=gamma)
    states.est = optimalViterbi$Final_Path
    
    # check if assingment is valid. if not start with random assingment
    if (length(unique(states.est)) != K || sum(diff(states.est)!=0)<3 ) {
      print(paste('Produced single state path at iteration', it))
      states.est = ICC.shuffle(K = K, Time = Time)
    }
    
    # M - step // update cluster parameters.
    T.error = try(
      {
        theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
      }
      , silent =T)
    
    if(class(T.error) == "try-error"){
      print(paste("Try error at iteration:", it, "Reshuffling states"))
      #print(T.error)
      ErrorRate = ErrorRate + 1
      states.est = ICC.shuffle(K = K, Time = Time)
      #break
    }
    
    distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
    d.hist =c(d.hist, -1*ICC.TotalDistance(distance.est, gamma = gamma, stateSeq = states.est))
    
    plot(y = distance.est[1,], x = 1:Time, type='l', col = 'red', main = paste('iteration: ', it))
    lines(y = distance.est[2,], x = 1:Time, col = 'blue')
    
    if (nrow(distance.est) == 3 ){ lines(y = distance.est[3,], x = 1:Time, col = 'black')}
    if (nrow(distance.est) == 4 ){ lines(y = distance.est[3,], x = 1:Time, col = 'yellow')}
    if (nrow(distance.est) == 5 ){ lines(y = distance.est[3,], x = 1:Time, col = 'pink')}
    
    #print(paste('Iteration', it))
    # check convergence
    #print(paste('Distance after iteration', it ,':', d.hist[length(d.hist)]))
    #if (it < 5 ){
    #  next
    #}
    #ds = diff(d.hist)
    #ds = ds[length(ds)] / d.hist[length(d.hist)-1]
    #if( abs(ds) < abs(0.05)){
    #  print('convergence')
    #  plot(y = d.hist, x = 1:it, type ='l', main = 'Distance history')
    #  break
    #}
  }
  #plot(y = d.hist, x = 1:it, type ='l', main = 'Distance history')
  return( list('OptimalViterbi' = optimalViterbi, 'OptimalPath' = optimalViterbi$Final_Path,  'ThetaEst' = theta.est, 'ER' = round(100*(ErrorRate/max.iters),2)) )
}

# ----------------- helper funcions -----------------
# shuffles observations into state sequences
ICC.shuffle = function(K=2, Time = 100){
  return(sample(1:K, Time, replace = T))
}

#K = K
#returns = returnsMatrix
#stateSeq = states.est
#sparse = sparseMethod

# evaluates sample statistics for each set of states
ICC.thetaEst = function(K=2, returns = matrix(runif(250), 5, 10), stateSeq = c(1,1,1,1,2,2,2,2,1,1), sparse = 1){
  thetaEsta.mu = list()
  thetaEsta.precision = list()
  for( state  in 1:K){
    temp.index = which(stateSeq==state) 
    temp.sample = returns[,temp.index]
    
    # compute mu vector
    thetaEsta.mu[[state]] = matrix(rowMeans(temp.sample), nrow = nrow(returns), ncol = 1, byrow = T)
    # compute precision
    thetaEsta.precision[[state]] = switch(sparse,
                                solve(abs( cov( t(temp.sample) ) ) ), # full precision
                                LoGo.solve(abs(cov( t( temp.sample ) ) ) ), # TMFG-LoGo
                                cov(t(temp.sample))) #Glasso - not currently used
  }
  return(list('mu' = thetaEsta.mu ,  'precision' = thetaEsta.precision))
}

#returns = returnsMatrix
#dist = distanceFunction
#K = K
#theta = theta.est

# evaluate disance matrix given returns, theta (estimate parameters) and state sequence
ICC.distance = function(returns, dist, theta, K){
  # compute D (distance matrix)
  Distance.matrix = matrix(NA, nrow = K, ncol = ncol(returns))
  for (state in 1:K){
    Distance.matrix[state, ]  = switch(dist,
                                      abs( mahalanobis(t(returns), center = c(theta$mu[[state]]), cov = theta$precision[[state]], inverted = T) ),
                                      1, # euclidian distance - not currently used
                                      1  # exponential family lilkihood - not currently used
    )
  }
  return(Distance.matrix)
}

# Not currentyl used.
#given a distance matrix and gamma, the function evaluates the total distance of the path
ICC.TotalDistance = function(distance, gamma, stateSeq){
  totalDist = distance[stateSeq[1], 1]
  Time = ncol(distance)
  for (obs in 2:Time){
    if( stateSeq[obs-1] == stateSeq[obs]){
      totalDist = totalDist + distance[stateSeq[obs], obs]
    }else{
      totalDist = totalDist + distance[stateSeq[obs], obs] + gamma
    }
  }
  return(totalDist)
}



#-----------------------------
#-----------------------------
#-----------------------------Testing
#-----------------------------


#naSums = function(x){sum(is.na(x))}

# SnP 500 data clean -------------------------------------
#allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
#allStocks = allStocks  %>% dplyr::select(date, close, Name)
#flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
#survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
#smaller = sample(2:400, size = 50)
#GRet = survivorStocks[, c(smaller)]
#GRet
#GRet = diff(as.matrix(log(GRet)))
#colnames(GRet)


# test variables

#hist(abs(diag((GRet - colMeans(GRet)) %*% (LoGo(GRet)) %*% t(GRet- colMeans(GRet)))), xlab = 'distance', main ='Distribution of Mahalanobis distance')

#dim(GRet)
#par(mfrow = c(1,1))
#plot(1:1258, y =(100*cumsum(colMeans(t(GRet))) + 100) , type='l', ylab = 'cummulative return', xlab = 'time', main ='SnP 500 sample return')
#plot(1:1258, y =colMeans(t(GRet)) , type='l', ylab = 'daily return', xlab = 'time', main ='SnP 500 sample return')

#par(mfrow=c(1,1))
# ICC setup
#gamma = 16
#sparseMethod = 2
#distanceFunction = 1 
#K =2


#ICC.Output = ICC.cluster(returnsMatrix = t(GRet), sparseMethod = sparseMethod, gamma = gamma, K = K, max.iters = 50)

#FinalPath = ICC.Output$OptimalPath
#Number.switches = sum(diff(FinalPath) != 0)
#Number.switches

#par(mfrow =c(1,2))
# visualise identified path
#plotter.series(colMeans(t(GRet)),FinalPath
#               , title = paste('Optimal market states, gamma:', gamma))

# visualise identified path
#plotter.series(colMeans(t(GRet)),FinalPath
#               , title = paste('Optimal market states, gamma:', gamma), S0 = 100)


#det(ICC.Output$ThetaEst$precision[[2]])
