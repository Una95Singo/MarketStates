# segmentation procedure ( modified EM algorithm )
# una singo, SNGUNA003
# 27 January 2020

# libraries ----------
source('R/viterbi.R')
source('R/plotter.R')
source('R/simulate.R')
library('NetworkToolbox')


#' segmentation procedure
#' @param returns matrix of reutrns stocks by time
#' @param K number of states
#' @param Time Time horizon of the dataset
#' @param Sparse type of sparsifiation method
#' @param dist which distance metric to use


set.seed(1)
# # test 1 step function - 2 states 30/50/20
#simulated = simulate.step_states(130, c(1,2,1), c(0.3,0.5,0.2))
#plotter.states(simulated, title = 'simulated states', thickness = 4)

#returns = twoStateData$stockMatrix

#seg.results = segmentation.procedure(returns = returns, Time = 130, gamma=0.001, iters =300, K = 2 )
#plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')


#minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

#minIndex

#plotter.states(seg.results$history$`iter: 31 states`, thickness = 4)

#seg.results$history$`iter: 31 states`
#simulated
#returnsMatrix = t(GRet)
#K = 2
#sparseMethod = 1
#distanceFunction = 1
#gamma = 0.00001

segmentation.procedure  = function(returnsMatrix, iters = 1, gamma = 1, sparseMethod = 1, distanceFunction = 1, K){
  trace(print)
  Time = ncol(returnsMatrix)
  totalDistance.history = c()
  viterbiHistory = list()
  totalDistance.historyDelta = c()
  optimalStates = NA
  
  # intialise --------
  print(paste("-------Iteration:",0, "------"))
  print(paste("-------Initialise step-------"))  
  # randomly assign obseravtions to states
  states.est = segmentation.shuffle(K = K, Time = Time)

  for (iter in 1:iters){
    print(paste("-------Iteration:", iter, "------")) 
    
    print(paste("-------Expectation step-------"))  
    theta.est = segmentation.thetaEst(K = K, returns = returnsMatrix, stateSeq = states.est, sparse = sparseMethod)
    distance.est = segmentation.distance(returns = returnsMatrix, dist = distanceFunction, theta = theta.est, K = K)
    
    print(paste("-------Maximisation step-------"))
    viterbiOptimal = viterbi(D=t(distance.est), K=K, gamma=gamma)
    # check if valid, if not valid reshuffle and move to next  
    if (length(unique(viterbiOptimal$Final_Path)) != K){
      print(paste("invalid states:", length(unique(viterbiOptimal$Final_Path))))
      states.est = segmentation.shuffle(K = K, Time = Time)
      next 
    }
    # check if it reduces overall distance
    distance.new = viterbiOptimal$Final_Cost
    if ( iter > 1){
      if(distance.new < distance.old){
        print("less than old")
        totalDistance.history = c(totalDistance.history, distance.new)
        final.theta = theta.est
        final.Viterbi = viterbiOptimal
        distance.old = distance.new
        states.est = viterbiOptimal$Final_Path
      }else{
        print("larger than old")
        states.est = segmentation.shuffle(K = K, Time = Time)
      }
      
    }else{
      totalDistance.history = c(totalDistance.history, distance.new)
      distance.old = distance.new
      states.est = viterbiOptimal$Final_Path
    }
  }
  print(totalDistance.history)
  return( list('history' = totalDistance.history, 'currentTheta' = final.theta, 'currentViterbi' = final.Viterbi))
}

# helper funcions -----------------
# shuffles observations into state sequences
segmentation.shuffle = function(K=2, Time = 100){
  return(sample(1:K, Time, replace = T))
}

# evaluates sample statistics for each set of states
segmentation.thetaEst = function(K=2, returns = matrix(runif(250), 5, 10), stateSeq = c(1,1,1,1,2,2,2,2,1,1), sparse = 1){
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

# evaluate disance matrix given returns, theta (estimate parameters) and state sequence
segmentation.distance = function(returns, dist, theta, K){
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
segmentation.TotalDistance = function(distance, gamma, stateSeq, Time){
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
