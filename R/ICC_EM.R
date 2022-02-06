# ICC implementation
# una singo, SNGUNA003
# 27 January 2020

# libraries ----------
source('R/viterbi.R')
library('NetworkToolbox')
library('Matrix')

#returnsMatrix = t(GRet)
#gamma = 24
#sparseMethod = 2
#distanceFunction = 1
#K =2

#returnsMatrix = t(GRet)
#gamma = 0.05
#sparseMethod = 2
#distanceFunction = 1 
#K =2



#set.seed(1)

ICC.cluster  = function(returnsMatrix, gamma = 1, sparseMethod = 1, distanceFunction = 1, K, max.iters = 15){
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
  max.it = max.iters
  stop.crit = 0.00001
  ds = NA
  d.hist= c()
  #initial
  states.est = ICC.shuffle(K = K, Time = Time)
  for (it in 1:max.it){
    #print("trap 2")
    # E - step 
    theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
    distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
    par(mfrow = c(1,1))
    plot(y = distance.est[1,], x = 1:Time, type='l', col = 'red')
    lines(y = distance.est[2,], x = 1:Time, col = 'blue')
    # M - step
    optimalViterbi = viterbi(D = t(distance.est), K=K, gamma=gamma)
    d.hist = c(d.hist, optimalViterbi$Final_Cost)
    states.est = optimalViterbi$Final_Path
    
    if (length(unique(states.est)) != K || sum(diff(states.est)!=0)<3 ) {
      print(paste('Produced single state path at iteration', it))
      #print(states.est)
      break
    }
    print(paste('Distance after iteration', it ,':', d.hist[length(d.hist)]))
    if (it == 1 ){
      next
    }
    ds = diff(d.hist)
    ds = ds[length(ds)] / d.hist[length(d.hist)-1]
    if( abs(ds) < abs(stop.crit)){
      print('convergence')
      plot(y = d.hist, x = 1:it, type ='l', main = 'Distance history')
      break
    }
  }
  
  return( list('OptimalViterbi' = optimalViterbi, 'OptimalPath' = optimalViterbi$Final_Path,  'ThetaEst' = theta.est))
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
    #print(temp.sample)
    # compute mu vector
    thetaEsta.mu[[state]] = matrix(rowMeans(temp.sample), nrow = nrow(returns), ncol = 1, byrow = T)
    # compute precision
    thetaEsta.precision[[state]] = switch(sparse,
                                solve(cov(t(temp.sample))), # full precision
                                LoGo(cov(t(temp.sample)) , na.data = 'none', partial = F), #cov(t(temp.sample)), # Tringle
                                cov(t(temp.sample))) #Glasso 
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
                                      1, # euclidian distance
                                      1  # exponential family lilkihood
    )
  }
  return(Distance.matrix)
}
#given a distance matrix and gamma, the function evaluates the total distance of the path
ICC.TotalDistance = function(distance, gamma, stateSeq, Time){
  totalDist = distance[stateSeq[1], 1]
  for (obs in 2:Time){
    if( stateSeq[obs-1] == stateSeq[obs]){
      totalDist = totalDist + distance[stateSeq[obs], obs]
    }else{
      totalDist = totalDist + distance[stateSeq[obs], obs] + gamma
    }
  }
  return(totalDist)
}
