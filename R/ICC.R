# ICC implementation
# una singo, SNGUNA003
# 27 January 2020

# libraries ----------
source('R/viterbi.R')
library('NetworkToolbox')


#returnsMatrix = returns
#gamma = 0
#sparseMethod = 1
#distanceFunction = 1
#K =2

#returnsMatrix = t(GRet)
#gamma = Gamma[i]
#sparseMethod = 1 
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
  
  # continue until valid sequence returned
  while( optimalStateNum != K ){
    print(paste("Iteration: ", iters))
    print(paste("-------Segmentation------"))
    states.est = ICC.shuffle(K = K, Time = Time)
    theta.est = ICC.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
    distance.est = ICC.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
    #browse()
    print(paste("-------Cluster assignment------"))
    optimalViterbi = viterbi(D=t(distance.est), K=K, gamma=gamma)
    optimalStateNum = length(unique(optimalViterbi$Final_Path))
    iters = iters + 1
    print(optimalStateNum)
    if (iters == max.iters){
      print(paste("-------Unable to maximise number of states for given gamma------"))
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
                                cov(t(temp.sample)), # full precision
                                LoGo(t(temp.sample), partial = F),#cov(t(temp.sample)), # Tringle
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
                                      abs(diag(t(returns - c(theta$mu[[state]])) %*% theta$precision[[state]] %*% (returns - c(theta$mu[[state]])))), #mahalanobis distance
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
