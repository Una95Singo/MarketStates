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

#returns =  rbind(simulated*0.2, simulated*0.3)

#seg.results = segmentation.procedure(returns = returns, Time = 130, gamma=0.001, iters =300, K = 2 )
#plot(y=seg.results$history$metric, x =1:length(seg.results$history$metric), main = 'history', type ='l', xlab = 'iterations', ylab  = 'dist')


#minIndex = which(seg.results$history$metric ==min(seg.results$history$metric))

#minIndex

#plotter.states(seg.results$history$`iter: 31 states`, thickness = 4)

#seg.results$history$`iter: 31 states`
#simulated

#Time = 1258 
#gamma=Gamma[10]
#iters =10
#K = 2
#sparse = 2
#dist = 1
#returns = t(GRet)


segmentation.procedure  = function(returns, Time, iters = 1, gamma = 1, sparse = 1, dist = 1, K){
  # instance variables ---------
  #-----------------------------
  seg.history = list() # stores history of computed distances, states and sample statistics
  seg.current = list() # list containing all current outputs
  seg.history[[paste('metric')]] = c()
  nStocks = nrow(returns)
  #---------------------
  # initialise ---------
  #---------------------
  # randomly assign each multivariate observation to a state
  seg.current[['states']] = sample(1:K, Time, replace = T) 
  # setup distance (D) matrix  
  seg.current[['distance']] = matrix(NA, nrow = K, ncol = Time)
  # compute sample estimates given staes
  for( state  in 1:K){
    temp.index = which(seg.current[['states']]==state) 
    temp.sample = returns[,temp.index]
    # compute mu vector
    seg.current[[paste('mu', state)]] = matrix(rowMeans(temp.sample), nrow = nStocks, ncol = 1, byrow = T)
    # compute precision
    seg.current[[paste('precision', state)]] = switch(sparse, 
                                                      cov(t(temp.sample)), # full precision
                                                      LoGo(cov(t(temp.sample)), partial = F),#cov(t(temp.sample)), # Tringle
                                                      cov(t(temp.sample))) #Glasso 
    # compute D (distance matrix)
    seg.current[['distance']][state,] = switch(dist,
                                               diag(t(returns - c(seg.current[[paste('mu', state)]])) %*% seg.current[[paste('precision', state)]] %*% (returns - c(seg.current[[paste('mu', state)]]))), #mahalanobis distance
                                               1, # euclidian distance
                                               1  # exponential family lilkihood
    )
    
  }
  
  # plotpoint
  plotter.states(seg.current$states, title = paste('iteration:', 0))
  iter_count = 0
  for ( iter in 1:iters){
    #---------------------
    # E-Step -------------
    #---------------------
    temp.viterbi = viterbi(D=t(seg.current$distance), K=K, gamma=gamma)
    temp.unique = length(unique(temp.viterbi$Final_Path))
    
    # plotpoint
    #plotter.states(temp.viterbi$Final_Path, title = paste('iteration:', iter))
    
    #---------------------
    # M-Step -------------
    #---------------------
    
    if (temp.unique != K ){ #resuffle
      # store history-----
      #seg.history[['metric']] = c(seg.history[['metric']],temp.viterbi$Final_Cost)
      #seg.history[[paste('iter:',iter,'states')]] = seg.current$states
      #for (state in 1:K){
      # compute mu vector
      #  seg.history[[paste('iter:', iter, 'mu', state)]] = seg.current[[paste('mu', state)]]
      #  seg.history[[paste('iter:', iter, 'precision', state)]] = seg.current[[paste('precision', state)]]
      
      #}
      
      # randomly assign each multivariate observation to a state
      seg.current[['states']] = sample(1:K, Time, replace = T) 
      # setup distance (D) matrix  
      seg.current[['distance']] = matrix(NA, nrow = K, ncol = Time)
      # compute sample estimates given staes
      for( state  in 1:K){
        temp.index = which(seg.current[['states']]==state) 
        temp.sample = returns[,temp.index]
        # compute mu vector
        seg.current[[paste('mu', state)]] = matrix(rowMeans(temp.sample), nrow = nStocks, ncol = 1, byrow = T)
        # compute precision
        seg.current[[paste('precision', state)]] = switch(sparse, 
                                                          cov(t(temp.sample)), # full precision
                                                          LoGo(cov(t(temp.sample)), partial = F),#cov(t(temp.sample)), # Tringle
                                                          cov(t(temp.sample))) #Glasso 
        # compute D (distance matrix)
        seg.current[['distance']][state,] = switch(dist,
                                                   diag(t(returns - c(seg.current[[paste('mu', state)]])) %*% seg.current[[paste('precision', state)]] %*% (returns - c(seg.current[[paste('mu', state)]]))), #mahalanobis distance
                                                   1, # euclidian distance
                                                   1  # exponential family lilkihood
        )
      }
    }
    else{
      iter_count = iter_count + 1
      plotter.states(temp.viterbi$Final_Path, title = paste('iteration:', iter_count))
      # store history-----
      seg.history[['metric']] = c(seg.history[['metric']],temp.viterbi$Final_Cost)
      seg.history[[paste('iter:',iter_count,'states')]] = temp.viterbi$Final_Path
      for (state in 1:K){
        # compute mu vector
        seg.history[[paste('iter:', iter_count, 'mu', state)]] = seg.current[[paste('mu', state)]]
        seg.history[[paste('iter:', iter_count, 'precision', state)]] = seg.current[[paste('precision', state)]]
        
      }
      
      # update sample estimates
      
      seg.current[['states']] = temp.viterbi$Final_Path
      # setup distance (D) matrix  
      seg.current[['distance']] = matrix(NA, nrow = K, ncol = Time)
      # compute sample estimates given staes
      for( state  in 1:K){
        temp.index = which(seg.current[['states']]==state) 
        temp.sample = returns[,temp.index]
        # compute mu vector
        seg.current[[paste('mu', state)]] = matrix(rowMeans(temp.sample), nrow = nStocks, ncol = 1, byrow = T)
        # compute precision
        seg.current[[paste('precision', state)]] = switch(sparse, 
                                                          cov(t(temp.sample)), # full precision
                                                          LoGo(cov(t(temp.sample)), partial = F),#cov(t(temp.sample)), # Tringle
                                                          cov(t(temp.sample))) #Glasso 
        # compute D (distance matrix)
        seg.current[['distance']][state,] = switch(dist,
                                                   diag(t(returns - c(seg.current[[paste('mu', state)]])) %*% seg.current[[paste('precision', state)]] %*% (returns - c(seg.current[[paste('mu', state)]]))), #mahalanobis distance
                                                   1, # euclidian distance
                                                   1  # exponential family lilkihood
        )
        
      }
      
      # plotpoint
      #plotter.states(seg.current$states, title = paste('iteration:', iter))
    }
  }
  return( list('history' = seg.history, 'current' = seg.current))
}


